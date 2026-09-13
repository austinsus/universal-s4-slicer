# S4 Pipeline

One command from an STL to ready-to-print **5-axis nonplanar G-code**:

```sh
python s4_pipeline.py model.stl
```

The pipeline deforms the model so overhangs become printable without support,
slices the deformed model with stock CuraEngine, maps the planar G-code back
through the inverse deformation into 4-axis "core R-theta" G-code, and finally
converts that into pose G-code (`X Y Z A B`) for a Stewart-platform printer —
or, via a small plug-in file, for **your** 5-axis kinematics.

The deformation stage is built on [S4_Slicer](https://github.com/jyjblrd/S4_Slicer)
by Joshua Bird ([video](https://www.youtube.com/watch?v=M51bMMVWbC8)) — this
repo turns that notebook into a headless end-to-end toolchain with mesh
repair, a printer config file, and a swappable kinematics back end.
License: GPL-3.0, like the original.

---

## Setup

1. **Get the code** — clone or download this repository.

2. **Python 3.10+ with the dependencies** (3.12 tested):

   ```sh
   python3 -m venv .venv
   .venv/bin/pip install -r requirements.txt
   ```

3. **CuraEngine** (used headlessly for stage 2):
   - macOS: install [UltiMaker Cura](https://ultimaker.com/software/ultimaker-cura/)
     — the bundled engine is found automatically inside the app.
   - Anything else: point the pipeline at the binary with
     `--cura-engine /path/to/CuraEngine` or `export CURA_ENGINE=...`
     (a `CuraEngine` on `PATH` is also picked up).

   The slicing profile lives in `cura_config.3mf` (a normal Cura project file:
   open it in Cura, tweak, re-save) with your machine-specific settings applied
   on top from `printer.json`.

4. **Describe your printer** — edit [`printer.json`](printer.json)
   (see [Configuring your printer](#configuring-your-printer)).

5. Run it:

   ```sh
   .venv/bin/python s4_pipeline.py "examples/... or your model.stl"
   ```

   Outputs land in `pipeline_out/<model name>/`; the file to print is
   `stewart_tilt-signed.gcode` (final pose G-code). `--help` lists every option.

---

## How it works

```
model.stl
   │  0. isotropic remesh (pymeshlab) - evenly sized triangles
   ▼
[1. deform]   tetrahedralize (TetGen; repair fallbacks if the STL is broken),
   │          optimize a per-tet rotation field so downward-facing surfaces
   │          rotate up to the printable overhang angle, then least-squares
   │          deform the mesh to realize that field  -> deformed.stl
   ▼
[2. slice]    headless CuraEngine slice of the DEFORMED model with the
   │          profile from cura_config.3mf + printer.json overrides
   │          -> sliced_planar.gcode
   ▼
[3. invmap]   map every planar G-code segment back through the inverse
   │          deformation: flat layers become curved, tilted layers on the
   │          real part -> s4_4axis.gcode  (G01 C.. X.. Z.. B.. E.. F..)
   ▼
[4. stewart]  kinematics conversion to machine pose G-code, clearance
   │          checking, feeds, fan, start/end G-code
   │          -> stewart_tilt-signed.gcode      (swap in your own machine
   ▼              here: see Custom kinematics)
[5. check]    optional slow full collision scan (--check)
```

Stage 1 details worth knowing:

- **Every input is isotropically remeshed first** (evenly sized, near
  equilateral triangles). This makes tetrahedralization dramatically more
  robust and keeps the tet count sane. `--no-remesh` skips it.
- **Broken STLs are handled.** If TetGen rejects or crashes on the surface
  (self-intersections, non-manifold edges — very common with sculpt/CAD
  exports), the pipeline walks a repair ladder automatically:
  direct → pymeshfix repair → signed-distance-field rebuild of a single
  watertight shell (inside/outside decided by a 6-direction ray-parity vote,
  robust against holes) → SDF + pymeshfix. You'll see which rung worked in
  the log.
- The deformation and the 4-axis output clamp layer tilt to the
  **rotation limit** (default `auto` = `max_tilt + max_clearance`), so the
  sliced file never demands steeper layers than the machine can print.

## Mesh resolution options

| option | default | effect |
|---|---|---|
| *(none)* | remesh at the STL's own triangle count | best fidelity; deform time grows with size |
| `--decimate N` | off | remesh to ~N triangles (≈8000 is plenty for most parts; use ~16000 for detailed models). This is the main speed/quality knob. |
| `--no-remesh` | remesh on | use the raw STL triangles; with `--decimate N` it falls back to plain quadric decimation |
| `--max-tet-volume MM3` | off | densify the tet mesh for low-poly models (a cube gives ~20 tets otherwise — too coarse for a smooth rotation field) |

Guideline: the deform stage cost scales with tet count; 15–35k tets
(≈8–16k triangles) deforms in 1–4 minutes and prints indistinguishably from
full resolution for most models.

---

## Configuring your printer

Everything machine-specific is in **[`printer.json`](printer.json)**:

```jsonc
{
  "machine": {              // defaults for the CLI options of the same name
    "mode": "tilt-signed",  // pose mode (see below)
    "max_tilt": 20.0,       // deg - see "The two key angles"
    "max_clearance": 25.0,  // deg
    "nozzle_offset": 42.0,  // mm, pivot-to-nozzle-tip arm of your head/platform
    "out_travel_hop": 2.0,  // mm hop added to travel moves
    "bed_safe_z": 1.0,      // mm below which the nozzle stays vertical
    "tilt_ramp": 0.5,       // mm band over which tilt fades back in above that
    "bed_phase": "off",     // "auto": keep whole print level until the last
                            //  bed-adjacent extrusion (models that only touch
                            //  the bed early); "off": per-move ramp only
    "speed_scale": 1.0,     // global print-speed multiplier
    "print_cap": 10.0,      // mm/s cap for print moves
    "temp": 200.0           // hotend temperature in the emitted header
  },
  "start_gcode": [ "...", "M109 S{temp:g}", "..." ],   // replaces the header
  "end_gcode":   [ "...", "G0 Z{lift_z:.1f}", "..." ], // replaces the footer
  "cura_overrides": {       // any Cura setting, applied on top of the .3mf
    "machine_nozzle_size": 0.4,
    "machine_width": 200, "machine_depth": 200, "machine_height": 150,
    "machine_center_is_zero": true
  }
}
```

- Any `machine` value can still be overridden per run on the command line
  (`--max-tilt 30` beats the file).
- A different file can be selected with `--printer other-printer.json` — keep
  one JSON per machine.
- **Start/end G-code**: plain lists of lines. `{temp}` is available in
  `start_gcode`, `{lift_z}` (highest printed Z + 10) in `end_gcode`; write
  literal `{`/`}` as `{{`/`}}`. One-off experiments can use
  `--start-gcode file.gcode` / `--end-gcode file.gcode` instead.

### The two key angles

- **`max_tilt`** — how far your platform/head can physically tilt. S4 layer
  tilts beyond this are clamped to it.
- **`max_clearance`** — the cone angle of free space around the nozzle tip:
  the largest angle between the nozzle axis and the local layer normal before
  the hotend body hits the printed part. Whatever tilt the machine could not
  reproduce (clamped by `max_tilt`, or levelled near the bed) is left as
  exactly this angle — the converter counts every extruding move that exceeds
  it and prints the worst offenders:

  ```
  worst nozzle-to-layer-normal angle: 24.9 deg (limit 25)
  ```

  A few violations in the first millimetre above the bed are usually
  acceptable (the near-bed levelling trade-off); violations mid-print mean
  you should lower `--rotation-limit`, raise `max_tilt`, or re-run with a
  shallower `--max-overhang`.
- **`rotation_limit`** *(CLI: `--rotation-limit`)* — clamp on layer tilt in
  the deformation itself. Default `auto` = `max_tilt + max_clearance`
  (never produce layers the machine can't print), a number for an explicit
  clamp, `off` for the unclamped original-notebook behaviour.

### Pose modes (`mode`)

- `tilt-signed` *(default)* — signed tilt with one continuous azimuth
  formula; needs shortest-path rotary handling in the firmware
  (see [docs/FIRMWARE-POSE-CHANGES.md](docs/FIRMWARE-POSE-CHANGES.md)).
- `tilt` — unsigned tilt, azimuth follows the bed angle; wraps at ±180°.
- `3axis` — vertical nozzle, nonplanar XYZ only (works on any printer that
  accepts plain `X Y Z E F` — a normal 3-axis machine with good Z travel).

---

## Custom kinematics

Stage 4 is machine-agnostic except for one small back end
(`StewartKinematics` in [`s4_to_stewart.py`](s4_to_stewart.py)). The
converter hands it, per move, the true nozzle-tip position `(x, y, z)`, the
local layer-normal tilt `b_deg` and bed angle `c_deg`; the back end returns
the pose your machine uses and formats the output line. Feeds, fan,
retraction, travel hops and clearance accounting stay in the converter.

To target your machine:

1. Copy [`kinematics_example.py`](kinematics_example.py) (documented
   template) and adapt `pose()` and `format_move()` — e.g. emit `B C` rotary
   words instead of `A B`, apply your rotary sign conventions, or resolve a
   different head geometry.
2. Run with `--kinematics my_machine.py` (works on `s4_pipeline.py` and
   standalone `s4_to_stewart.py`).

`start_lines()`/`end_lines()` can also be overridden in code, but prefer
`printer.json` for plain start/end G-code.

---

## Examples

[`examples/`](examples/) contains two complete slices produced by this exact
pipeline (deformed model + final pose G-code + the command used):

- **`pi-3mm/`** — a π sign, default settings end to end.
- **`christmas-tree/`** — heavy overhangs everywhere; also a demo of the
  mesh-repair ladder (the source STL is self-intersecting with 192 open
  edges) and of `--decimate 16000` for a detailed model.

## Repository layout

| file | role |
|---|---|
| `s4_pipeline.py` | the pipeline (stages 1–5), CLI |
| `s4_to_stewart.py` | stage 4: S4 4-axis → machine pose G-code (also standalone) |
| `check_nozzle_collisions.py` | stage 5: slow geometric collision scan |
| `kinematics_example.py` | template for your own machine back end |
| `printer.json` | your machine: limits, start/end G-code, Cura overrides |
| `cura_config.3mf` | Cura project holding the slicing profile |
| `docs/FIRMWARE-POSE-CHANGES.md` | firmware requirements for `tilt-signed` |

## Troubleshooting

- **`CuraEngine not found`** — set `--cura-engine` or `CURA_ENGINE` (see Setup).
- **`WARNING: N extruding moves exceed the max clearance angle`** — see
  [The two key angles](#the-two-key-angles).
- **`tetgen crashed` / repair-ladder messages** — informational; only a final
  `could not tetrahedralize this model` is fatal. If the SDF rung reports a
  volume very different from your model, inspect `deformed.stl` before printing.
- **Deform stage is slow** — pass `--decimate 8000` (see Mesh resolution).
- Every stage can be re-run in isolation: `--skip-deform` reuses a previous
  deformation, `--from-gcode` injects a hand-sliced file, `--stop-after`
  stops early. Artifacts are plain files in `pipeline_out/<name>/`.

## Credits & license

- Deformation approach and inverse map: [S4_Slicer](https://github.com/jyjblrd/S4_Slicer)
  © Joshua Bird, GPL-3.0.
- Pipeline, mesh repair, Cura integration, Stewart conversion: Austin Huang.
- License: [GPL-3.0](LICENSE).

```bibtex
@software{Bird_S4_Slicer,
  author = {Bird, Joshua}, license = {GPL-3.0},
  title = {{S4 Slicer}}, url = {https://github.com/jyjblrd/S4_Slicer}
}
```
