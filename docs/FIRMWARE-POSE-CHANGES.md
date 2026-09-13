# Stewart pose-handling changes for RRF (LinearStewartKinematics)

Spec for firmware changes to the custom RepRapFirmware build at
`~/Documents/Projects/5 Axis Printer/Final_Printer/Firmwares/RRF/RRFBuild-3.6.3`.
Target file: `RepRapFirmware/src/Movement/Kinematics/LinearStewartKinematics.cpp`
(+ its header, and one `config.g` line). Written for an implementer with no
prior context on this machine.

## Machine background (read first)

- 6-carriage linear Stewart platform ("hexapod"). The user commands a 5-DOF
  nozzle pose via axes **X Y Z A B**: XYZ = nozzle tip in the *plate frame*
  (Z = height above the plate), **A = tilt theta (deg)**, **B = azimuth phi
  (deg)** of the build direction. C is a redundant sixth carriage axis that is
  never commanded and is deliberately exempt from M208 limiting (see the
  comment in `LimitPosition`, ~line 990).
- Kinematics is **segmented** (`M669 K15 ... S200 T0.2`): moves are chopped
  into ~5 ms segments; each segment target runs through the IK (`Transform`)
  which yields the six carriage heights. This is what keeps the pose - and the
  plate-to-nozzle relationship - controlled along a move. **Do not remove any
  axis from segmentation.**
- The pose parametrization is redundant. For direction
  `d(theta,phi) = (sin·cos, sin·sin, cos)` the complete same-pose set is
  exactly:
  1. `phi + 360k` (wrap ladder),
  2. `(-theta, phi + 180 + 360k)` (twin family),
  3. `theta = 0`: phi is fully degenerate (every phi is the same pose).
  There are no other cases; the changes below handle all three.
- The FK (forward transform, canonicalization around line ~455) picks between
  twin representations by continuity with `lastFkPose`, and the current
  machine config limits A to `M208 A0:45`. A previous incident: after a move
  abort, the FK reported the (-20, 0) twin of a commanded (20, 180) pose,
  which violated `M208 A0:45` and poisoned every subsequent move
  ("intermediate position outside machine limits" cascade). The owner has
  since added handling here - **do not rewrite the canonicalization; read what
  is there and keep it working**.

## Changes (in recommended order)

### 1. Allow signed tilt: `M208 A-45:45`

One line in `config.g` (`M208 A0:45` -> `M208 A-45:45`). Then, in the
firmware, audit everything that assumes theta >= 0:

- the FK validity filter (`FkMinValidTheta` / `FkMaxValidTheta`, used ~line
  470): must accept negative theta symmetrically;
- any clamp or assert on THETA_AXIS.

Negative tilt is physically valid: `(-t, phi)` is identical to `(t, phi+180)`.
The IK (`Transform`) already handles it through sin/cos. This change alone
makes the twin-cascade class structurally impossible.

### 2. Pose target normalization hook (the core change)

Add one function, called **once per move at target-set time, before
segmentation** (NOT per segment - find where the kinematics receives the
user-space target; on this build the owner has already modified
`src/GCodes/GCodes.cpp` and `Kinematics.cpp/.h`, so follow their pattern for
where kinematics-specific target processing can live). Given the current pose
`(a_c, b_c)` and commanded target `(a_t, b_t)`:

```
candidates = [ (a_t,  wrapShortest(b_t,  b_c)),          // as commanded
               (-a_t, wrapShortest(b_t + 180, b_c)) ]    // twin
```

where `wrapShortest(target, current)` returns
`current + NormaliseAngle(target - current)` (NormaliseAngle is at ~line 404,
maps to [-180, 180)). Selection:

- **Degenerate target**: if `|a_t| < 1 deg`, the target phi is meaningless:
  use `(sign-preserving a_t, b_c)` - i.e. B does not move at all. Skip the
  rest.
- **Degenerate current**: if `|a_c| < 1 deg`, B may be re-aimed for free -
  prefer the candidate with the smaller |A| motion and set its phi directly
  (B is physically a no-op at vertical, but still let it interpolate; it
  costs nothing).
- Otherwise compute the cost of each candidate as **actual carriage travel**:
  run `Transform()` on the current pose and on each candidate (at the move's
  commanded XYZ) and sum |delta carriage height| over the six towers. Pick
  the cheaper candidate **only if it wins by at least 2x**; on anything
  closer, keep the pose as commanded (prevents representation chatter on
  paths that oscillate near the tie boundary).

The chosen candidate replaces the user-space A/B target for the move, so
segmentation interpolates the short/cheap way.

Acceptance behavior:
- from (A20, B0), `G1 A20 B180` must **rock through vertical** (A sweeps
  +20 -> -20, B stays 0), not precess 180 deg of azimuth;
- from (A20, B179), `G1 A20 B-179` must move B by ~2 deg, not ~358;
- from (A20, B0), `G1 A15 B10` must behave exactly as commanded (direct);
- from (A0, B<anything>), `G1 A20 B90` must not sweep B while tilting -
  B goes (freely) to 90, tilt rises in that plane.

### 3. Exempt B from M208 limiting; fold it between moves

With wrapShortest the commanded B coordinate accumulates outside +/-180.
- In `LimitPosition` (~line 990): skip M208 clamping for PHI_AXIS, with a
  comment mirroring the existing C-axis exemption ("phi is periodic; every
  value is physically legal").
- After each move completes (or at a safe idle point), fold the stored user
  coordinate of B back into [-180, 180) - a position rewrite, not motion -
  so M114/DWC stay readable and the coordinate cannot grow without bound.
  The FK resync already reports folded values; with the hook in change 2
  operating relative to the *current* coordinate, a fold is transparent.

### 4. (Optional, quality-of-life) Pose axes excluded from feed pacing

RRF paces a move by total coordinate length, counting A/B degrees as mm, so
combined position+pose moves run the tip slower than F says. If addressed:
pace the move duration by XYZ length only and let A/B track along the
segments; clamp to per-carriage speed limits as usual (a huge pose delta on
a tiny XY move must still stretch the move rather than violate carriage
limits). This is independent of changes 1-3 and can be skipped initially -
G-code generators currently compensate by scaling F.

## Test plan (on the machine, no filament)

1. `M208 A-45:45` in config; restart; `G1 A-10 B0 F3000` from level - must
   tilt opposite to `A10 B0` and `M114` must report A-10 without error.
2. The four acceptance moves from change 2, watching the platform.
3. Wrap fold: `G1 A10 B170`, then `G1 A10 B-170` (expect ~20 deg short arc),
   then `M114` - B in [-180, 180).
4. Abort mid-move (M112 or pause) while tilted, `M114`, then command the same
   pose again - no cascade, no giant swing (this exercises the FK-resync +
   fold + relative-wrapShortest interplay; it is the scenario that previously
   bricked a print).
5. A full print of an existing known-good file (e.g. a TigerSlicer `-a20`
   file): behavior must be unchanged - those files use A in [0, 20] with
   isolated pose moves; the selection logic must leave them alone.

## Explicitly out of scope

- Do not modify the FK canonicalization logic (owner-maintained).
- Do not change segmentation parameters (M669 S/T).
- Do not touch the plate-mesh (M669 C coefficients) or M666 rail terms.
