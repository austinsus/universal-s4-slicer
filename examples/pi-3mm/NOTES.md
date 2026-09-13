# pi 3mm — default settings

    python s4_pipeline.py "pi 3mm.stl"

Everything default: isotropic remesh at the model's own resolution (870
triangles), rotation limit auto (45°), tilt-signed pose mode with the
printer.json in the repo root.

Result: 90,977 moves, worst nozzle-to-layer-normal angle 11.4° (limit 25°),
zero clearance violations.
