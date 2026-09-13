#!/usr/bin/env python3
"""Check Stewart-platform pose G-code (output of s4_to_stewart.py) for spots
where the nozzle/hotend intersects material that was already printed.

Model: the nozzle tip is at the commanded XYZ; the nozzle axis in the part
frame leans by A degrees toward azimuth B (the same convention s4_to_stewart
emits: B is the azimuth the nozzle leans toward; A=0 or no A word = vertical).
Everything the hotend occupies is modeled as a single clearance cone: the
volume above the plane through the tip perpendicular to the axis, rising at
--clearance degrees. A previously deposited point whose elevation angle above
that tip plane exceeds --clearance is a (potential) collision; how bad is
measured by the elevation excess and the penetration depth.

This matches s4_to_stewart's --max-clearance parameter: the report tells you
what clearance angle your hotend actually needs for a given file.
"""
import argparse
import math
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

WORD = re.compile(r"([A-Z])(-?\d+\.?\d*)")


def parse_moves(path):
    """Yield (line_no, x, y, z, a_deg, b_deg, extruding) for every move with
    a full position. A/B persist modally from the last move that set them."""
    x = y = z = None
    a = b = 0.0
    moves = []
    with open(path) as fh:
        for ln, raw in enumerate(fh):
            line = raw.strip()
            if not line.startswith(("G0 ", "G1 ", "G0\t", "G1\t")):
                continue
            words = dict((k, float(v)) for k, v in WORD.findall(line.split(";")[0]))
            x = words.get("X", x)
            y = words.get("Y", y)
            z = words.get("Z", z)
            a = words.get("A", a)
            b = words.get("B", b)
            if x is None or y is None or z is None:
                continue
            e = words.get("E", 0.0)
            moves.append((ln, x, y, z, a, b, e > 0))
    return moves


def axis_from_pose(a_deg, b_deg):
    a = math.radians(a_deg)
    p = math.radians(b_deg)
    return np.array([math.sin(a) * math.cos(p),
                     math.sin(a) * math.sin(p),
                     math.cos(a)])


def main():
    ap = argparse.ArgumentParser(
        description="Find nozzle-vs-printed-material collisions in Stewart pose G-code.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("gcode")
    ap.add_argument("-c", "--clearance", type=float, default=25.0, metavar="DEG",
                    help="hotend clearance angle: material must stay below this "
                         "elevation above the nozzle tip plane")
    ap.add_argument("--max-range", type=float, default=25.0, metavar="MM",
                    help="ignore deposited material farther than this from the tip "
                         "(set to the tip-to-effector-plate distance)")
    ap.add_argument("--min-dist", type=float, default=2.0, metavar="MM",
                    help="ignore material closer than this to the tip (the bead "
                         "being laid down and its immediate neighborhood)")
    ap.add_argument("--lag", type=int, default=150, metavar="MOVES",
                    help="ignore material deposited within this many moves before "
                         "the current one (the still-active path)")
    ap.add_argument("--sample", type=float, default=1.0, metavar="MM",
                    help="spacing of deposited-material sample points along extrusions")
    ap.add_argument("--check-every", type=int, default=1, metavar="N",
                    help="check every Nth move (speedup for huge files)")
    ap.add_argument("--top", type=int, default=12, help="collision clusters to print")
    args = ap.parse_args()

    moves = parse_moves(args.gcode)
    print(f"{len(moves)} positioned moves parsed from {args.gcode}")

    # --- deposited material samples (chronological) ---
    dep_pts, dep_move, dep_line = [], [], []
    prev = None
    n_phantom = 0
    for mi, (ln, x, y, z, a, b, ext) in enumerate(moves):
        p = np.array([x, y, z])
        if ext and prev is not None:
            seg = p - prev
            L = float(np.linalg.norm(seg))
            # prime/transition artifacts (e.g. the S4 file's first E move)
            # interpolate as steep diagonals through mid-air; no real bead
            # climbs anywhere near 60 deg
            if L > 1e-9 and abs(seg[2]) / L > math.sin(math.radians(60)):
                n_phantom += 1
                prev = p
                continue
            n = max(1, int(L / args.sample))
            for k in range(1, n + 1):
                dep_pts.append(prev + seg * (k / n))
                dep_move.append(mi)
                dep_line.append(ln)
        prev = p
    dep_pts = np.asarray(dep_pts)
    dep_move = np.asarray(dep_move)
    dep_line = np.asarray(dep_line)
    print(f"{len(dep_pts)} deposited sample points "
          f"(every {args.sample} mm along extrusions; "
          f"{n_phantom} steep prime/transition segments excluded)")

    tree = cKDTree(dep_pts)
    tan_clear = math.tan(math.radians(args.clearance))

    hits = []          # (excess_deg, penetration_mm, move_idx, dep_idx)
    worst_elev_all = 0.0
    checked = 0
    for mi in range(0, len(moves), args.check_every):
        ln, x, y, z, a, b, ext = moves[mi]
        t = np.array([x, y, z])
        u = axis_from_pose(a, b)
        idx = tree.query_ball_point(t, args.max_range)
        if not idx:
            continue
        idx = np.asarray(idx)
        idx = idx[dep_move[idx] < mi - args.lag]
        if len(idx) == 0:
            continue
        checked += 1
        d = dep_pts[idx] - t
        dist = np.linalg.norm(d, axis=1)
        keep = dist > args.min_dist
        if not keep.any():
            continue
        idx, d, dist = idx[keep], d[keep], dist[keep]
        axial = d @ u
        radial = np.linalg.norm(d - np.outer(axial, u), axis=1)
        elev = np.degrees(np.arctan2(axial, radial))
        wi = int(np.argmax(elev))
        if elev[wi] > worst_elev_all:
            worst_elev_all = float(elev[wi])
        viol = elev > args.clearance
        if viol.any():
            # penetration: distance from point to the cone surface (approx)
            pen = dist[viol] * np.sin(np.radians(elev[viol] - args.clearance))
            for j, dj in zip(np.nonzero(viol)[0], pen):
                hits.append((float(elev[j] - args.clearance), float(dj),
                             mi, int(idx[j])))

    print(f"\nclearance angle tested: {args.clearance} deg  "
          f"(range {args.min_dist}-{args.max_range} mm from tip, lag {args.lag} moves)")
    print(f"worst material elevation above tip plane anywhere: {worst_elev_all:.1f} deg")
    if not hits:
        print("no collisions found.")
        return

    # --- cluster by tip position (5 mm voxels), keep worst hit per cluster ---
    clusters = {}
    n_travel = 0
    for exc, pen, mi, di in hits:
        ln, x, y, z, a, b, ext = moves[mi]
        if not ext:
            n_travel += 1
        key = (round(x / 5), round(y / 5), round(z / 5))
        cur = clusters.get(key)
        if cur is None or pen > cur[1]:
            clusters[key] = (exc, pen, mi, di)
        else:
            clusters[key] = cur
    print(f"{len(hits)} violating (move, material-point) pairs "
          f"({n_travel} on travel moves) in {len(clusters)} 5mm regions\n")

    print(f"top {min(args.top, len(clusters))} collision regions "
          f"(worst penetration first):")
    ranked = sorted(clusters.values(), key=lambda h: -h[1])[:args.top]
    for exc, pen, mi, di in ranked:
        ln, x, y, z, a, b, ext = moves[mi]
        r = math.hypot(x, y)
        dz = dep_pts[di]
        print(f"  line {ln:>7}  tip({x:7.2f},{y:7.2f},{z:6.2f}) r={r:5.1f}  "
              f"A={a:5.1f} B={b:7.1f} {'print ' if ext else 'TRAVEL'} "
              f"pen {pen:5.2f} mm (+{exc:4.1f} deg over) vs material from "
              f"line {dep_line[di]} at z={dz[2]:.2f}")

    # --- correlation: violation rate by tip radius and by commanded tilt ---
    def bin_stats(getter, edges, label):
        tot = {}
        bad = {}
        vmoves = set(mi for _, _, mi, _ in hits)
        for mi in range(0, len(moves), args.check_every):
            v = getter(moves[mi])
            for lo, hi in edges:
                if lo <= v < hi:
                    tot[(lo, hi)] = tot.get((lo, hi), 0) + 1
                    if mi in vmoves:
                        bad[(lo, hi)] = bad.get((lo, hi), 0) + 1
        print(f"\ncolliding-move rate by {label}:")
        for k in edges:
            t = tot.get(k, 0)
            if t:
                b_ = bad.get(k, 0)
                print(f"  {k[0]:>5g}-{k[1]:<5g}: {b_:>6}/{t:<7} ({100*b_/t:5.1f}%)")

    bin_stats(lambda m: math.hypot(m[1], m[2]),
              [(0, 5), (5, 10), (10, 15), (15, 20), (20, 30)], "tip radius (mm)")
    bin_stats(lambda m: abs(m[4]),
              [(0, 5), (5, 10), (10, 15), (15, 20.01)], "commanded tilt |A| (deg)")

    req = args.clearance + max(h[0] for h in hits)
    print(f"\nthis file needs a hotend clearance angle of at least {req:.1f} deg "
          f"(worst elevation seen {worst_elev_all:.1f} deg)")


if __name__ == "__main__":
    sys.exit(main())
