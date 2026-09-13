# christmas tree — detailed model with a broken source mesh

    python s4_pipeline.py "christmass tree.stl" --decimate 16000

The source STL (334k triangles) is self-intersecting with 192 open edges, so
the deform stage's repair ladder ends up on its last rung (SDF rebuild +
pymeshfix — watch for "tetrahedralized via ..." in the log). --decimate 16000
keeps more detail than the ~8000 default guideline; the result is 33,278 tets.

Result: 245k moves, worst nozzle-to-layer-normal angle 24.9° (limit 25°),
zero clearance violations — heavy overhangs printed by tilting up to 45°.
