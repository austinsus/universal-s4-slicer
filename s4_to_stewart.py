#!/usr/bin/env python3
"""Convert S4_Slicer (Joshua Bird core R-theta 4-axis) G-code to Stewart-platform
pose G-code (X Y Z A B, absolute plate-frame coordinates, G94 feeds).

S4 output format per move:  G00/G01 C<deg,unwrapped> X<radial> Z<height> B<tilt deg> [E<rel>] F<...>
  - F is inverse-time (G93, 1/min) except in G94 spans (travel/retract lines, mm/min)
  - The slicer pre-compensated for the nozzle-pivot offset; travels also carry a
    baked hop. Both are undone here to recover the true tip position.

Modes:
  3axis (default) : A0 B0 everywhere - nonplanar XYZ with a vertical nozzle.
  tilt            : A = clamp(|B_s4|, 0, max_tilt), azimuth follows the bed angle.
                    WARNING: the azimuth wraps at +/-180; if the firmware does not
                    treat B as a continuous rotary, each wrap commands a 360 swing.
  tilt-signed     : A = clamp(B_s4, -max_tilt, max_tilt) with its sign; phi is the
                    single continuous formula C+180 (no flips at tilt-sign changes).
                    REQUIRES the firmware changes in FIRMWARE-POSE-CHANGES.md
                    (M208 A-45:45, shortest-path B). Do NOT print on stock firmware.

Clearance check:
  The S4 B word is the tilt of the local layer normal. Whatever tilt the pose
  cannot reproduce (|B_s4| - |A_commanded|) is left as an angle between the
  nozzle axis and the layer normal. If that residual exceeds --max-clearance,
  the hotend body can hit the already-printed part. Policy per --clearance-policy:
    warn  (default) : count + report the worst offenders, emit the move anyway
    mark            : same as warn, plus a ";CLEARANCE" comment on each move
    skip            : drop the move's pose down to the clearance limit is NOT
                      possible (the path itself is baked in), so the move is
                      emitted with A at max tilt but flagged; use warn/mark and
                      re-slice with a smaller S4 rotation range instead
    error           : abort the conversion at the first violation
"""
import argparse
import importlib.util
import json
import math
import os
import re
import sys

WORD = re.compile(r"([A-Z])(-?\d+\.?\d*)")

# Machine start/end G-code. Override per printer in printer.json
# ("start_gcode" / "end_gcode" lists) or with --start-gcode/--end-gcode files.
# Available placeholders: {temp} in start lines, {lift_z} in end lines
# (write literal braces as {{ }}).
DEFAULT_START_GCODE = [
    "G21 ; mm",
    "G90 ; absolute coordinates",
    "M83 ; relative extrusion",
    "G28 ; home",
    "M109 S{temp:g} ; heat and wait",
    "G1 E10 F300 ; prime (as in the source)",
    "G0 Z20 A0 B0 F6000 ; start pose",
]
DEFAULT_END_GCODE = [
    "M107 ; fan off",
    "M104 S0 ; heater off",
    "G0 Z{lift_z:.1f} F6000 ; lift",
    "G0 A0 B0 F6000 ; level",
    "M84",
]


def norm180(a):
    a = math.fmod(a + 180.0, 360.0)
    if a < 0:
        a += 360.0
    return a - 180.0


class StewartKinematics:
    """Default pose backend: Stewart-platform "X Y Z A B" G-code.

    This class is the machine-specific half of the converter; convert() owns
    everything machine-agnostic (parsing the S4 file, undoing the baked
    pivot/hop compensation to recover the true tip position, feeds, fan,
    clearance accounting). To target a different 5-axis machine, copy
    kinematics_example.py, adapt these four methods, and pass the file with
    --kinematics (or s4_pipeline.py --kinematics):

      start_lines(cfg)                          header G-code lines
      pose(cfg, cmd, x, y, z, b_deg, c_deg,
           in_bed_phase, prev_pose)             -> (a_out, phi, capped, levelled)
           map the S4 layer tilt b_deg / bed angle c_deg at tip position
           (x, y, z) to your machine's tilt+azimuth (or whatever your axes
           mean); capped/levelled are stats flags
      format_move(cfg, cmd, x, y, z, a_out,
                  phi, e, f_cmd)                -> one G-code line
      end_lines(cfg, z_max)                     footer G-code lines
    """

    def start_lines(self, cfg):
        lines = cfg.start_gcode if cfg.start_gcode else DEFAULT_START_GCODE
        return [ln.format(temp=cfg.temp) for ln in lines]

    def end_lines(self, cfg, z_max):
        lines = cfg.end_gcode if cfg.end_gcode else DEFAULT_END_GCODE
        return [ln.format(lift_z=(z_max or 0) + 10) for ln in lines]

    def pose(self, cfg, cmd, x, y, z, b_deg, c_deg, in_bed_phase, prev_pose):
        mode = cfg.mode
        capped = levelled = False
        if in_bed_phase and mode in ("tilt", "tilt-signed"):
            a_out, phi = 0.0, prev_pose[1]
        elif mode == "tilt-signed":
            # signed tilt, one continuous azimuth formula: (b, C+180) is
            # the same pose as (-b, C); sign changes rock through vertical
            a_out = max(-cfg.max_tilt, min(b_deg, cfg.max_tilt))
            if abs(b_deg) > cfg.max_tilt:
                capped = True
            a_allowed = cfg.max_tilt * min(max((z - cfg.bed_safe_z) / cfg.tilt_ramp, 0.0), 1.0)
            if abs(a_out) > a_allowed:
                a_out = a_allowed if a_out > 0 else -a_allowed
                levelled = True
            if abs(a_out) < cfg.tilt_epsilon:
                a_out, phi = 0.0, prev_pose[1]
            else:
                # Empirically verified on the machine: the correct emission
                # is B = C + 180. (The contract's "machine B = model phi -
                # 180" combined with the geometric derivation gave B = C,
                # which tilted the wrong way - the model-frame assumption
                # and the S4 frame differ by 180 somewhere in the chain.)
                phi = norm180(c_deg + 180.0)
        elif mode == "tilt" and abs(b_deg) >= cfg.tilt_epsilon:
            a_out = min(abs(b_deg), cfg.max_tilt)
            if abs(b_deg) > cfg.max_tilt:
                capped = True
            # near the bed the nozzle stays vertical; the path (layer
            # orientation) is untouched, only the pose is levelled
            a_allowed = cfg.max_tilt * min(max((z - cfg.bed_safe_z) / cfg.tilt_ramp, 0.0), 1.0)
            if a_out > a_allowed:
                a_out = a_allowed
                levelled = True
            if a_out < cfg.tilt_epsilon:
                a_out, phi = 0.0, prev_pose[1]
            else:
                phi = norm180(c_deg + (180.0 if b_deg > 0 else 0.0))
        else:
            a_out, phi = 0.0, prev_pose[1]
        return a_out, phi, capped, levelled

    def format_move(self, cfg, cmd, x, y, z, a_out, phi, e, f_cmd):
        parts = [cmd, f"X{x:.3f}", f"Y{y:.3f}", f"Z{z:.3f}"]
        if cfg.mode in ("tilt", "tilt-signed"):
            parts += [f"A{a_out:.2f}", f"B{phi:.2f}"]
        if e is not None:
            parts.append(f"E{e:.4f}")
        parts.append(f"F{f_cmd:.0f}")
        return " ".join(parts)


def load_kinematics(path):
    """Import a python file that defines get_kinematics() -> backend object
    (see StewartKinematics for the required methods)."""
    spec = importlib.util.spec_from_file_location("custom_kinematics", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not hasattr(mod, "get_kinematics"):
        sys.exit(f"--kinematics {path}: file must define get_kinematics()")
    return mod.get_kinematics()


def load_printer_config(path, required=False):
    """Read printer.json: {"machine": {...option defaults...},
    "start_gcode": [...], "end_gcode": [...], "cura_overrides": {...}}."""
    if path is None:
        return {}
    if not os.path.exists(path):
        if required:
            sys.exit(f"printer config not found: {path}")
        return {}
    with open(path) as fh:
        cfg = json.load(fh)
    return {k: v for k, v in cfg.items() if not k.startswith("_")}


def apply_machine_defaults(parser, machine, strict=False):
    """Use printer.json's "machine" section as argparse defaults (CLI flags
    still win). strict=True rejects unknown keys (the pipeline, which knows
    every machine option, runs strict so typos don't pass silently); the
    standalone converter only warns, since some keys belong to the pipeline."""
    dests = {a.dest for a in parser._actions}
    unknown = sorted(set(machine) - dests)
    if unknown and strict:
        sys.exit(f"printer config: unknown machine option(s): "
                 f"{', '.join(unknown)}")
    for key in unknown:
        print(f"note: machine option {key!r} is not one of this tool's "
              f"options - ignored here", file=sys.stderr)
    parser.set_defaults(**{k: v for k, v in machine.items() if k in dests})


def build_parser():
    p = argparse.ArgumentParser(
        description="Convert S4_Slicer 4-axis G-code to Stewart-platform pose G-code.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("src", help="S4_Slicer output .gcode file")
    p.add_argument("-o", "--output", default=None,
                   help="output path (default: <src>-stewart-<mode>[...].gcode)")
    p.add_argument("-m", "--mode", choices=("3axis", "tilt", "tilt-signed"),
                   default="3axis", help="pose mode")

    g = p.add_argument_group("main pose parameters")
    g.add_argument("-a", "--max-tilt", type=float, default=20.0, metavar="DEG",
                   help="max nozzle tilt A; S4 tilts beyond this are clamped")
    g.add_argument("-c", "--max-clearance", type=float, default=25.0, metavar="DEG",
                   help="max nozzle clearance angle: largest allowed angle between "
                        "the nozzle axis and the local layer normal (the tilt the "
                        "pose could not reproduce). Violations mean the hotend can "
                        "collide with the printed part.")
    g.add_argument("--clearance-policy", choices=("warn", "mark", "error"),
                   default="warn", help="what to do when a move violates --max-clearance")

    g = p.add_argument_group("geometry (match your S4 slice settings)")
    g.add_argument("--nozzle-offset", type=float, default=42.0, metavar="MM",
                   help="S4 pivot-to-tip arm baked into the source coordinates")
    g.add_argument("--travel-hop", type=float, default=1.0, metavar="MM",
                   help="hop baked into S4 travel (G0) moves")
    g.add_argument("--out-travel-hop", type=float, default=1.0, metavar="MM",
                   help="hop re-added to emitted travel moves. The original "
                        "TigerSlicer converter removed the baked hop and never "
                        "re-added one, so travels dragged at surface height; "
                        "0 reproduces that behavior")

    g = p.add_argument_group("near-bed behavior")
    g.add_argument("--bed-safe-z", type=float, default=1.0, metavar="MM",
                   help="below this tip height the nozzle stays vertical "
                        "(tilting near the bed brings the platform edge into the nozzle body)")
    g.add_argument("--tilt-ramp", type=float, default=1.0, metavar="MM",
                   help="height band above --bed-safe-z over which tilt fades back in")
    g.add_argument("--tilt-epsilon", type=float, default=2.0, metavar="DEG",
                   help="below this tilt A stays 0 (avoids azimuth flapping)")
    g.add_argument("--bed-phase", choices=("auto", "off"), default="auto",
                   help="auto: keep the whole print level until the last "
                        "source line that extrudes below --bed-safe-z (good "
                        "when only the first layers touch the bed); off: rely "
                        "on the per-move z ramp only (use when nonplanar "
                        "layers dip to the bed deep into the print, where "
                        "auto would keep a third of the file vertical)")

    g = p.add_argument_group("feeds")
    g.add_argument("--speed-scale", type=float, default=1.0,
                   help="global speed multiplier (print moves only)")
    g.add_argument("--print-cap", type=float, default=None, metavar="MM/S",
                   help="cap print-move speed at this many mm/s")
    g.add_argument("--feed-cap", type=float, default=10200.0, metavar="MM/MIN",
                   help="absolute feed ceiling for every emitted move")
    g.add_argument("--feed-fallback", type=float, default=6000.0, metavar="MM/MIN",
                   help="feed used when a move has no usable length/feed")

    g = p.add_argument_group("start/end & fan")
    g.add_argument("--temp", type=float, default=200.0, metavar="C",
                   help="M109 hotend temperature in the emitted header")
    g.add_argument("--no-fan", action="store_true",
                   help="do not emit automatic M106/M107 fan toggles")
    g.add_argument("--fan-speed", type=int, default=255, metavar="0-255",
                   help="fan PWM used when the fan is toggled on")

    g = p.add_argument_group("printer config & customization")
    g.add_argument("--printer", default=None, metavar="JSON",
                   help="printer.json whose 'machine' section provides "
                        "defaults for the options above and whose "
                        "'start_gcode'/'end_gcode' lists replace the "
                        "built-in header/footer")
    g.add_argument("--start-gcode", default=None, metavar="FILE",
                   help="file of G-code lines replacing the default start "
                        "sequence ({temp} placeholder available); wins over "
                        "printer.json")
    g.add_argument("--end-gcode", default=None, metavar="FILE",
                   help="file of G-code lines replacing the default end "
                        "sequence ({lift_z} placeholder available); wins "
                        "over printer.json")
    g.add_argument("--kinematics", default=None, metavar="PY",
                   help="python file defining get_kinematics() -> object "
                        "with the StewartKinematics interface, to emit "
                        "G-code for a different 5-axis machine")
    return p


class ClearanceError(RuntimeError):
    pass


def convert(cfg, kin=None):
    if kin is None:
        kin = (load_kinematics(cfg.kinematics) if cfg.kinematics
               else StewartKinematics())
    src, mode = cfg.src, cfg.mode
    arm0 = cfg.nozzle_offset

    # First pass: find the last source line that extrudes below bed_safe_z.
    # Until that point the whole print stays at A0 ("first layers are level"),
    # because in a nonplanar file bed-adjacent paths are interleaved with
    # higher ones - a per-move height rule tilts far too early.
    bed_phase_end = -1
    line_no = -1
    if cfg.bed_phase == "auto":
        with open(src) as fh:
            for raw in fh:
                line_no += 1
                if not raw.startswith(('G01', 'G1', 'G0', 'G00')):
                    continue
                m = re.search(r"\bB(-?[\d.]+)", raw)
                e = re.search(r"\bE(0?\.[\d]+|[1-9][\d.]*)", raw)
                zc = re.search(r"\bZ(-?[\d.]+)", raw)
                bx = re.search(r"\bX(-?[\d.]+)", raw)
                if e and zc and bx and m is not None:
                    b_r = math.radians(float(m.group(1)))
                    z_tip = float(zc.group(1)) - (math.cos(b_r) - 1.0) * arm0
                    if z_tip < cfg.bed_safe_z:
                        bed_phase_end = line_no

    g93 = False
    src_line_no = -1
    prev = None           # previous tip position (x, y, z)
    prev_pose = (0.0, 0.0)
    stats = {"moves": 0, "capped": 0, "wraps": 0, "r_max": 0.0,
             "z_min": None, "z_max": None, "skipped": 0, "bed_levelled": 0,
             "fan_toggles": 0, "clearance_violations": 0, "worst_clearance": 0.0}
    worst = []  # (residual_deg, src_line_no, z) worst clearance offenders
    out = [
        f"; converted from S4_Slicer output: {src}",
        f"; mode={mode}  nozzle_offset={cfg.nozzle_offset}  max_tilt={cfg.max_tilt}"
        f"  max_clearance={cfg.max_clearance}",
        "; NOTE: no heater commands in the source - preheat before starting",
    ] + kin.start_lines(cfg)
    fan_on = False
    with open(src) as fh:
        for raw in fh:
            src_line_no += 1
            line = raw.strip()
            if not line:
                continue
            if line.startswith("G93"):
                g93 = True
                continue
            if line.startswith("G94"):
                g93 = False
                continue
            if line.startswith(("G28", "M83", "G90", "G21", "G1 E10")):
                continue  # already in our header
            m = re.match(r"^(G0[01]?)\b", line)
            if not m or " C" not in f" {line}":
                # passthrough of anything that is not an S4 move line
                if line.startswith(("M", "G")):
                    out.append(line)
                continue
            cmd = "G0" if m.group(1) in ("G0", "G00") else "G1"
            words = dict((a, float(v)) for a, v in WORD.findall(line))
            if "C" not in words or "X" not in words or "Z" not in words:
                continue
            c_deg = words["C"]
            b_deg = words.get("B", 0.0)
            hop = cfg.travel_hop if cmd == "G0" else 0.0
            arm = arm0 + hop
            b_rad = math.radians(b_deg)
            # undo the slicer's pivot-offset compensation -> true tip position
            r_tip = words["X"] + math.sin(b_rad) * arm
            z_tip = words["Z"] - (math.cos(b_rad) - 1.0) * arm - hop
            c_rad = math.radians(c_deg)
            x = r_tip * math.cos(c_rad)
            y = r_tip * math.sin(c_rad)
            z = z_tip
            if cmd == "G0":
                z += cfg.out_travel_hop
            if z < -0.05:
                stats["skipped"] += 1
                continue
            in_bed_phase = src_line_no <= bed_phase_end
            a_out, phi, capped, levelled = kin.pose(
                cfg, cmd, x, y, z, b_deg, c_deg, in_bed_phase, prev_pose)
            if capped:
                stats["capped"] += 1
            if levelled:
                stats["bed_levelled"] += 1
            if abs(phi - prev_pose[1]) > 180.0:
                stats["wraps"] += 1

            # clearance check: the tilt the pose could not reproduce is the
            # angle left between the nozzle axis and the local layer normal
            clearance_flag = False
            if cmd == "G1" and "E" in words and words["E"] > 0:
                residual = abs(b_deg) - abs(a_out)
                if residual > stats["worst_clearance"]:
                    stats["worst_clearance"] = residual
                if residual > cfg.max_clearance:
                    stats["clearance_violations"] += 1
                    clearance_flag = True
                    if len(worst) < 5 or residual > worst[-1][0]:
                        worst.append((residual, src_line_no, z))
                        worst.sort(reverse=True)
                        del worst[5:]
                    if cfg.clearance_policy == "error":
                        raise ClearanceError(
                            f"source line {src_line_no}: layer normal tilted "
                            f"{abs(b_deg):.1f} deg but nozzle only reaches "
                            f"{abs(a_out):.1f} deg -> {residual:.1f} deg residual "
                            f"exceeds max clearance {cfg.max_clearance:g} deg. "
                            f"Raise --max-tilt/--max-clearance or re-slice with a "
                            f"smaller S4 rotation range.")

            # feed
            dist = 0.0
            if prev is not None:
                dist = math.dist(prev, (x, y, z))
            f_in = words.get("F")
            if g93 and f_in:
                f_base = dist * f_in if dist > 1e-9 else cfg.feed_fallback
            else:
                f_base = f_in if f_in else cfg.feed_fallback
            if cmd == "G1" and "E" in words:
                f_base *= cfg.speed_scale
                if cfg.print_cap is not None:
                    f_base = min(f_base, cfg.print_cap * 60.0)
            # NO pose-length feed compensation: the firmware paces moves by
            # sqrt(XYZ^2 + great-circle-direction^2) itself (see
            # FIRMWARE-POSE-CHANGES-IMPLEMENTED.md); scaling here would
            # double-count and run moves too fast.
            f_cmd = min(max(f_base, 60.0), cfg.feed_cap)
            move_line = kin.format_move(cfg, cmd, x, y, z, a_out, phi,
                                        words.get("E"), f_cmd)
            if clearance_flag and cfg.clearance_policy == "mark":
                move_line += " ;CLEARANCE"
            # Nonplanar paths revisit bed heights late in the file, so the fan
            # follows the nozzle: off near the bed for adhesion, on elsewhere.
            if not cfg.no_fan and cmd == "G1" and words.get("E", 0) > 0:
                if not fan_on and z >= cfg.bed_safe_z + 0.2:
                    out.append(f"M106 S{cfg.fan_speed}")
                    fan_on = True
                    stats["fan_toggles"] += 1
                elif fan_on and z < cfg.bed_safe_z - 0.2:
                    out.append("M107 ; near the bed")
                    fan_on = False
                    stats["fan_toggles"] += 1
            out.append(move_line)
            stats["moves"] += 1
            stats["r_max"] = max(stats["r_max"], math.hypot(x, y))
            stats["z_min"] = z if stats["z_min"] is None else min(stats["z_min"], z)
            stats["z_max"] = z if stats["z_max"] is None else max(stats["z_max"], z)
            prev = (x, y, z)
            prev_pose = (a_out, phi)
    out += kin.end_lines(cfg, stats["z_max"])
    with open(cfg.output, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return stats, worst


def main(argv=None):
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--printer", default=None)
    pre_cfg, _ = pre.parse_known_args(argv)
    parser = build_parser()
    printer = load_printer_config(pre_cfg.printer, required=True)
    apply_machine_defaults(parser, printer.get("machine", {}))
    cfg = parser.parse_args(argv)
    # start/end gcode: --start-gcode file > printer.json list > built-in
    for which in ("start_gcode", "end_gcode"):
        path = getattr(cfg, which)
        if path:
            with open(path) as fh:
                setattr(cfg, which, [ln.rstrip("\n") for ln in fh])
        else:
            setattr(cfg, which, printer.get(which))
    if cfg.output is None:
        suffix = f"-stewart-{cfg.mode}"
        if cfg.max_tilt != 20.0 and cfg.mode != "3axis":
            suffix += f"-a{cfg.max_tilt:g}"
        if cfg.speed_scale != 1.0:
            suffix += f"-x{cfg.speed_scale:g}"
        if cfg.print_cap:
            suffix += f"-cap{cfg.print_cap:g}mms"
        cfg.output = cfg.src.rsplit(".", 1)[0] + suffix + ".gcode"
    try:
        s, worst = convert(cfg)
    except ClearanceError as e:
        sys.exit(f"CLEARANCE ERROR: {e}")
    print(f"wrote {cfg.output}")
    print(f"moves: {s['moves']}  skipped(z<0): {s['skipped']}  tilt-capped: {s['capped']}  "
          f"azimuth wraps: {s['wraps']}  near-bed levelled: {s['bed_levelled']}  "
          f"fan toggles: {s['fan_toggles']}")
    print(f"max radius: {s['r_max']:.1f} mm   Z range: [{s['z_min']:.2f}, {s['z_max']:.2f}]")
    print(f"worst nozzle-to-layer-normal angle: {s['worst_clearance']:.1f} deg "
          f"(limit {cfg.max_clearance:g})")
    if s["clearance_violations"]:
        print(f"WARNING: {s['clearance_violations']} extruding moves exceed the "
              f"max clearance angle - the hotend body may hit the part.")
        for residual, ln, z in worst:
            print(f"  source line {ln}: residual {residual:.1f} deg at Z {z:.2f}")
        print("  Fix: raise --max-tilt, raise --max-clearance if your hotend "
              "geometry allows, or re-slice in S4 with a smaller rotation range "
              "(MAX_ROTATION / MIN_ROTATION).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
