#!/usr/bin/env python3
"""Template for a custom kinematics backend.

s4_to_stewart.py's convert() handles everything machine-agnostic: parsing the
S4 4-axis file, recovering the true nozzle-tip position (x, y, z) plus the
local layer tilt b_deg and bed angle c_deg per move, feeds, fan, retraction
and the clearance accounting. Your backend only decides (1) what pose your
machine uses to realize that tilt and (2) how a move line is written.

Use it with:
    python s4_pipeline.py model.stl --kinematics my_machine.py
    python s4_to_stewart.py file.gcode --kinematics my_machine.py

This example emits the same X Y Z A B pose G-code as the built-in Stewart
backend but shows every override point. Subclass StewartKinematics and only
override what differs, or write an independent class with the same four
methods.
"""
from s4_to_stewart import StewartKinematics, norm180


class MyMachineKinematics(StewartKinematics):
    # --- header/footer -----------------------------------------------------
    # Normally you customize these in printer.json (start_gcode/end_gcode)
    # instead; override in code only if they must depend on logic.
    # def start_lines(self, cfg): ...
    # def end_lines(self, cfg, z_max): ...

    # --- pose mapping -------------------------------------------------------
    def pose(self, cfg, cmd, x, y, z, b_deg, c_deg, in_bed_phase, prev_pose):
        """Map the S4 layer orientation to your machine's rotary axes.

        In:  tip position (x, y, z) [mm], layer-normal tilt b_deg [deg,
             signed], bed angle c_deg [deg, unwrapped], in_bed_phase (True
             while the whole print must stay level), prev_pose = the
             (a, phi) you returned for the previous move.
        Out: (a, phi, capped, levelled)
             a    = tilt your machine will actually use [deg]
             phi  = azimuth of that tilt [deg, -180..180]
             capped   = True if you clamped |b_deg| to the machine limit
             levelled = True if you reduced the tilt near the bed
        convert() uses |b_deg| - |a| as the residual clearance angle.
        """
        # Example: identical to the Stewart tilt-signed mapping.
        return super().pose(cfg, cmd, x, y, z, b_deg, c_deg,
                            in_bed_phase, prev_pose)

    # --- move formatting ----------------------------------------------------
    def format_move(self, cfg, cmd, x, y, z, a_out, phi, e, f_cmd):
        """One output G-code line for this move. e is None for non-extruding
        moves. Rename/reorder words, convert to your controller's rotary
        convention, or split into multiple lines (return one string with
        embedded newlines) here."""
        parts = [cmd, f"X{x:.3f}", f"Y{y:.3f}", f"Z{z:.3f}",
                 f"A{a_out:.2f}", f"B{phi:.2f}"]
        if e is not None:
            parts.append(f"E{e:.4f}")
        parts.append(f"F{f_cmd:.0f}")
        return " ".join(parts)


def get_kinematics():
    return MyMachineKinematics()
