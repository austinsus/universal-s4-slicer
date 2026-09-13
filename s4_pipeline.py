#!/usr/bin/env python3
"""One-command S4 -> Stewart pipeline.

    s4_pipeline.py model.stl [options]

Stages (artifacts land in --outdir, default pipeline_out/<name>/):
  1. deform   : S4_Slicer mesh deformation (tetrahedralize, rotation field,
                least-squares deformation).  Ported from S4_Slicer/main.ipynb.
  2. slice    : headless CuraEngine slice of the deformed mesh, using the
                profile stack inside S4_Slicer/cura_config.3mf with the
                frontend-style '=expression' settings resolved here.
  3. invmap   : map the planar G-code back through the inverse deformation
                to 4-axis S4 G-code (C/X/Z/B).  Ported from main.ipynb.
  4. stewart  : s4_to_stewart.py conversion to Stewart pose G-code.
  5. check    : check_nozzle_collisions.py scan of the result.

Key physical parameters:
  --max-tilt      A the platform can reach (deg)
  --max-clearance angle the hotend tolerates between nozzle axis and layer
                  normal (deg)
  --rotation-limit clamp on the deformation rotation field AND the S4 B
                  output. Default 'auto' = max_tilt + max_clearance, so the
                  sliced file never demands layers steeper than the machine
                  can physically print. Pass degrees for an explicit clamp,
                  or 'off' for unclamped deformation (original notebook).
  --decimate      target triangle count for the mesh that gets tetrahedralized.
                  OFF by default - the STL's full resolution is kept. Pass a
                  count (~8000) to speed up the deform stage.
  --no-remesh     every STL is first isotropically remeshed (pymeshlab) to
                  evenly sized triangles at the --decimate count (or its own
                  count); this flag skips that and uses the raw STL faces.
"""
import os
import sys

# re-exec under the bundled venv so any `python3 s4_pipeline.py` works
_VENV = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".venv")
_VENV_PY = os.path.join(_VENV, "bin", "python")
if sys.prefix != _VENV and os.path.exists(_VENV_PY):
    os.execv(_VENV_PY, [_VENV_PY] + sys.argv)

import argparse
import ast
import configparser
import json
import math
import pickle
import re
import shutil
import subprocess
import tempfile
import time
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s4_to_stewart  # noqa: E402

# prefer 5.8 (matches the engine cura_config.3mf was made with); 5.4 fallback
_ENGINE_CANDIDATES = [
    "/Applications/UltiMaker Cura 5.8.app/Contents/Resources/CuraEngine",
    "/Applications/UltiMaker Cura.app/Contents/MacOS/CuraEngine",
]
DEFAULT_CURA_ENGINE = (os.environ.get("CURA_ENGINE")
                       or next((p for p in _ENGINE_CANDIDATES
                                if os.path.exists(p)), None)
                       or shutil.which("CuraEngine")
                       or _ENGINE_CANDIDATES[-1])
DEFAULT_CURA_CONFIG = os.path.join(HERE, "cura_config.3mf")

UP = np.array([0, 0, 1])


# --------------------------------------------------------------------------
# Stage 1: mesh deformation (port of main.ipynb cells 2-12)
# --------------------------------------------------------------------------

_TETGEN_CHILD = r"""
import sys, numpy as np, tetgen
inp, out, maxvol = sys.argv[1:4]
d = np.load(inp)
tg = tetgen.TetGen(d["v"], d["f"])
if maxvol:
    tg.tetrahedralize(quality=True, fixedvolume=True, maxvolume=float(maxvol))
else:
    # -Y (no Steiner points on the boundary) keeps the input surface exactly and
    # avoids tetgen exploding a sliver-y surface into 100k+ tets; fall back to
    # the default switches if this surface can't be recovered without them
    try:
        tg.tetrahedralize(nobisect=True)
    except RuntimeError as e:
        print(f"tetgen -Y failed ({e}); retrying with default switches",
              file=sys.stderr)
        tg = tetgen.TetGen(d["v"], d["f"])
        tg.tetrahedralize()
tg.grid.save(out)
"""


def isotropic_remesh(verts, tris, target_tris, iterations=10):
    """pymeshlab isotropic explicit remeshing to ~target_tris evenly sized,
    near-equilateral triangles. The target edge length comes from the surface
    area: N equilateral triangles of side L cover N * sqrt(3)/4 * L^2."""
    import pymeshlab as ml
    ms = ml.MeshSet()
    ms.add_mesh(ml.Mesh(vertex_matrix=np.asarray(verts, dtype=float),
                        face_matrix=np.asarray(tris, dtype=np.int32)))
    area = ms.get_geometric_measures()["surface_area"]
    edge = float(np.sqrt(4.0 * area / (np.sqrt(3.0) * max(int(target_tris), 4))))
    ms.apply_filter("meshing_isotropic_explicit_remeshing",
                    iterations=iterations, targetlen=ml.PureValue(edge),
                    adaptive=False)
    out = ms.current_mesh()
    return out.vertex_matrix().copy(), out.face_matrix().copy()


def load_and_center_mesh(stl_path, max_tet_volume=None, decimate=None,
                         remesh=True):
    import open3d as o3d
    import tetgen

    mesh = o3d.io.read_triangle_mesh(stl_path)
    if len(mesh.vertices) == 0:
        sys.exit(f"could not read mesh from {stl_path}")
    mesh.remove_duplicated_vertices()
    orig_verts = np.asarray(mesh.vertices).copy()
    orig_tris = np.asarray(mesh.triangles).copy()

    # remesh target: --decimate if given, else keep the STL's own resolution
    remesh_target = decimate if decimate is not None else len(orig_tris)

    if remesh:
        n0 = len(mesh.triangles)
        rv, rf = isotropic_remesh(orig_verts, orig_tris, remesh_target)
        mesh = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(rv),
                                         o3d.utility.Vector3iVector(rf))
        mesh.remove_duplicated_vertices()
        mesh.remove_degenerate_triangles()
        print(f"[deform] remeshed {n0} -> {len(mesh.triangles)} triangles "
              f"(target {remesh_target})")
    # the remesh already lands on the target count; plain quadric decimation
    # is the no-remesh path (or a safety net if the remesh overshot badly)
    if decimate is not None and len(mesh.triangles) > 1.25 * decimate:
        n0 = len(mesh.triangles)
        mesh = mesh.simplify_quadric_decimation(target_number_of_triangles=decimate)
        mesh.remove_duplicated_vertices()
        mesh.remove_degenerate_triangles()
        print(f"[deform] decimated {n0} -> {len(mesh.triangles)} triangles")
    verts = np.asarray(mesh.vertices)
    tris = np.asarray(mesh.triangles)

    def tetra(v, f):
        # TetGen is C++ and SEGFAULTS (rather than raising) on self-intersecting
        # / non-manifold input, which would kill the whole pipeline before the
        # repair fallbacks below get a chance. Run it in a child process so a
        # crash turns into a RuntimeError. The child's cwd is a temp dir so the
        # _skipped.node/.face debris tetgen dumps on failure lands there.
        import pyvista as pv
        with tempfile.TemporaryDirectory(prefix="s4_tetgen_") as td:
            inp = os.path.join(td, "surf.npz")
            out = os.path.join(td, "tet.vtu")
            np.savez(inp, v=np.asarray(v, dtype=float), f=np.asarray(f))
            r = subprocess.run(
                [sys.executable, "-c", _TETGEN_CHILD, inp, out,
                 "" if max_tet_volume is None else repr(float(max_tet_volume))],
                cwd=td, capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(out):
                how = (f"signal {-r.returncode}" if r.returncode < 0
                       else f"exit {r.returncode}")
                tail = (r.stderr or r.stdout).strip().splitlines()[-1:]
                raise RuntimeError(f"tetgen crashed ({how})"
                                   + (f": {tail[0]}" if tail else ""))
            return pv.read(out)

    def meshfix(v, f):
        import pymeshfix
        fixer = pymeshfix.MeshFix(np.asarray(v, dtype=float), np.asarray(f))
        fixer.repair(joincomp=True, remove_smallest_components=False)
        print(f"[deform] meshfix: {len(f)} -> {len(fixer.faces)} triangles")
        if len(fixer.faces) < 0.3 * len(f):
            raise RuntimeError("meshfix discarded most of the mesh")
        return fixer.points, fixer.faces

    def sdf_remesh():
        # rebuild a single watertight surface from a signed-distance field of
        # the ORIGINAL mesh (sign via ray parity - robust for multi-shell /
        # self-intersecting sculpts), then decimate for tetgen
        import pyvista as pv
        scene = o3d.t.geometry.RaycastingScene()
        scene.add_triangles(o3d.t.geometry.TriangleMesh.from_legacy(
            o3d.geometry.TriangleMesh(
                o3d.utility.Vector3dVector(orig_verts),
                o3d.utility.Vector3iVector(orig_tris))))
        lo = orig_verts.min(axis=0) - 2.0
        hi = orig_verts.max(axis=0) + 2.0
        spacing = max(float((hi - lo).max()) / 160.0, 0.4)
        dims = (np.ceil((hi - lo) / spacing).astype(int) + 1).tolist()
        grid = pv.ImageData(dimensions=dims, spacing=(spacing,) * 3, origin=lo)
        pts = np.asarray(grid.points, dtype=np.float32)
        # sign by majority parity vote over 6 axis-aligned rays: a hole or
        # non-manifold patch in the source casts a "shadow streak" of wrongly
        # signed samples along a single ray direction (rod/debris artifacts in
        # the rebuilt surface); the other directions outvote it
        votes = np.zeros(len(pts), dtype=np.int8)
        for d in np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                           [0, -1, 0], [0, 0, 1], [0, 0, -1]], np.float32):
            rays = o3d.core.Tensor(np.hstack(
                [pts, np.broadcast_to(d, pts.shape)]).astype(np.float32))
            votes += (scene.count_intersections(rays).numpy() % 2).astype(np.int8)
        dist = scene.compute_distance(o3d.core.Tensor(pts)).numpy()
        grid["d"] = np.where(votes >= 4, -dist, dist)
        def tri_only(pd):
            # clean() turns degenerate triangles into line/vertex cells, which
            # decimate() rejects and tetgen can't use - keep polygons only
            return pv.PolyData(pd.points, faces=pd.faces).triangulate()

        surf = tri_only(grid.contour([0.0], scalars="d").clean())
        target = decimate or 8000
        sv, sf = np.asarray(surf.points), surf.faces.reshape(-1, 4)[:, 1:]
        if remesh:
            # even triangles keep tetgen from filling slivers with Steiner tets
            sv, sf = isotropic_remesh(sv, sf, target)
        elif surf.n_cells > target:
            surf = tri_only(surf.decimate(1.0 - target / surf.n_cells).clean())
            sv, sf = np.asarray(surf.points), surf.faces.reshape(-1, 4)[:, 1:]
        print(f"[deform] SDF remesh at {spacing:.2f} mm voxels -> "
              f"{len(sf)} triangles")
        return sv, sf

    tet = None
    remeshed = None
    attempts = [
        ("direct", lambda: (verts, tris)),
        ("meshfix repair", lambda: meshfix(verts, tris)),
        ("SDF remesh", sdf_remesh),
        ("SDF remesh + meshfix", lambda: meshfix(*remeshed)
            if remeshed is not None else meshfix(*sdf_remesh())),
    ]
    for label, get in attempts:
        try:
            got = get()
            if label == "SDF remesh":
                remeshed = got
            tet = tetra(*got)
            if label != "direct":
                print(f"[deform] tetrahedralized via {label}")
            break
        except RuntimeError as e:
            print(f"[deform] {label} failed: {e}")
    if tet is None:
        sys.exit("[deform] could not tetrahedralize this model")
    if tet.n_cells < 500:
        print(f"[deform] WARNING: only {tet.n_cells} tets - the rotation "
              f"field will be very coarse and tilt will bleed into bed "
              f"cells. Use --max-tet-volume (mm^3) to densify low-poly "
              f"meshes (a few thousand tets is a good target).")
    x_min, x_max, y_min, y_max, z_min, z_max = tet.bounds
    tet.points -= np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])
    return tet


def build_topology(tet):
    import networkx as nx

    cell_neighbour_dict = {nt: {c: [] for c in range(tet.number_of_cells)}
                           for nt in ["point", "edge", "face"]}
    for neighbour_type in ["point", "edge", "face"]:
        cell_neighbours = []
        for cell_index in range(tet.number_of_cells):
            for neighbour in tet.cell_neighbors(cell_index, f"{neighbour_type}s"):
                if neighbour > cell_index:
                    cell_neighbours.append((cell_index, neighbour))
        for c1, c2 in np.array(cell_neighbours):
            cell_neighbour_dict[neighbour_type][c1].append(c2)
            cell_neighbour_dict[neighbour_type][c2].append(c1)
        tet.field_data[f"cell_{neighbour_type}_neighbours"] = np.array(cell_neighbours)

    graph = nx.Graph()
    cell_centers = tet.cell_centers().points
    for edge in tet.field_data["cell_point_neighbours"]:
        distance = np.linalg.norm(cell_centers[edge[0]] - cell_centers[edge[1]])
        graph.add_weighted_edges_from([(edge[0], edge[1], distance)])
    return cell_neighbour_dict, graph


def encode_object(obj):
    import base64
    return base64.b64encode(pickle.dumps(obj)).decode("utf-8")


def decode_object(s):
    import base64
    return pickle.loads(base64.b64decode(s))


def update_tet_attributes(tet, topo_graph):
    import networkx as nx

    surface_mesh = tet.extract_surface()
    cell_to_face = decode_object(tet.field_data["cell_to_face"])

    cells = tet.cells.reshape(-1, 5)[:, 1:]
    tet.add_field_data(cells, "cells")
    tet.add_field_data(tet.points, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]
    tet.add_field_data(faces, "faces")
    tet.add_field_data(surface_mesh.points, "face_vertices")

    tet.cell_data["face_normal"] = np.full((tet.number_of_cells, 3), np.nan)
    surface_normals = surface_mesh.face_normals
    for cell_index, face_indices in cell_to_face.items():
        fn = surface_normals[face_indices]
        tet.cell_data["face_normal"][cell_index] = fn[np.argmin(fn[:, 2])]
    tet.cell_data["face_normal"] /= np.linalg.norm(
        tet.cell_data["face_normal"], axis=1)[:, None]

    tet.cell_data["face_center"] = np.full((tet.number_of_cells, 3), np.nan)
    surface_centers = surface_mesh.cell_centers().points
    for cell_index, face_indices in cell_to_face.items():
        fc = surface_centers[face_indices]
        tet.cell_data["face_center"][cell_index] = fc[np.argmin(fc[:, 2])]

    tet.cell_data["cell_center"] = tet.cell_centers().points

    bottom_threshold = np.nanmin(tet.cell_data["face_center"][:, 2]) + 0.3
    bottom_mask = tet.cell_data["face_center"][:, 2] < bottom_threshold
    tet.cell_data["is_bottom"] = bottom_mask
    bottom_cells = np.where(bottom_mask)[0]

    face_normals = tet.cell_data["face_normal"].copy()
    face_normals[bottom_mask] = np.nan
    tet.cell_data["overhang_angle"] = np.arccos(np.dot(face_normals, UP))

    overhang_direction = face_normals[:, :2].copy()
    overhang_direction /= np.linalg.norm(overhang_direction, axis=1)[:, None]
    tet.cell_data["overhang_direction"] = overhang_direction

    IN_AIR_THRESHOLD = 1
    tet.cell_data["in_air"] = np.full(tet.number_of_cells, False)
    _, paths_to_bottom = nx.multi_source_dijkstra(topo_graph, set(bottom_cells))
    tet.cell_data["path_to_bottom"] = np.full(
        (tet.number_of_cells, max(len(x) for x in paths_to_bottom.values())), -1)
    for cell_index, path in paths_to_bottom.items():
        tet.cell_data["path_to_bottom"][cell_index, :len(path)] = path
    for cell_index in range(tet.number_of_cells):
        path = paths_to_bottom.get(cell_index, [])
        if len(path) > 1:
            heights = tet.cell_data["cell_center"][path, 2]
            if np.any(heights > tet.cell_data["cell_center"][cell_index, 2]
                      + IN_AIR_THRESHOLD):
                tet.cell_data["in_air"][cell_index] = True
    return tet


def calculate_tet_attributes(tet, topo_graph):
    surface_mesh = tet.extract_surface()

    cells = tet.cells.reshape(-1, 5)[:, 1:]
    tet.add_field_data(cells, "cells")
    tet.add_field_data(tet.points, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]
    tet.add_field_data(faces, "faces")
    face_vertices = surface_mesh.points
    tet.add_field_data(face_vertices, "face_vertices")

    cell_to_face = {}
    face_to_cell = {i: [] for i in range(len(faces))}
    cell_to_face_vertices = {}
    for cvi, cv in enumerate(tet.field_data["cell_vertices"].reshape(-1, 3)):
        fvi = np.where((face_vertices == cv).all(axis=1))[0]
        if len(fvi) == 1:
            cell_to_face_vertices[cvi] = fvi[0]

    for cell_index, cell in enumerate(tet.field_data["cells"]):
        fvis = [cell_to_face_vertices[v] for v in cell if v in cell_to_face_vertices]
        if len(fvis) >= 3:
            extracted = surface_mesh.extract_points(fvis, adjacent_cells=False)
            if extracted.number_of_cells >= 1:
                cell_to_face[cell_index] = list(extracted.cell_data["vtkOriginalCellIds"])
                for fi in extracted.cell_data["vtkOriginalCellIds"]:
                    face_to_cell[fi].append(cell_index)

    tet.add_field_data(encode_object(cell_to_face), "cell_to_face")
    tet.add_field_data(encode_object(face_to_cell), "face_to_cell")
    tet.cell_data["has_face"] = np.zeros(tet.number_of_cells)
    for cell_index in cell_to_face:
        tet.cell_data["has_face"][cell_index] = 1

    tet = update_tet_attributes(tet, topo_graph)
    bottom_mask = tet.cell_data["is_bottom"]
    bottom_cells = np.where(bottom_mask)[0]
    tet.cell_data["overhang_angle"][bottom_cells] = np.nan
    return tet, bottom_mask, bottom_cells


def plane_fit(points):
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1))
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    return ctr, svd(np.dot(x, x.T))[0][:, -1]


def calculate_path_length_to_base_gradient(tet, topo, graph, bottom_cells,
                                           max_overhang, smoothing, set_zero):
    import networkx as nx

    gradient = np.zeros(tet.number_of_cells)
    dist_to_bottom = np.full(tet.number_of_cells, np.nan)
    distances, paths = nx.multi_source_dijkstra(graph, set(bottom_cells))
    closest_bottom = np.zeros(tet.number_of_cells, dtype=int)
    for ci in range(tet.number_of_cells):
        normal = tet.cell_data["face_normal"][ci]
        is_overhang = np.arccos(np.dot(normal, UP)) > np.deg2rad(90 + max_overhang)
        if is_overhang and ci not in bottom_cells:
            closest_bottom[ci] = paths[ci][0]
            dist_to_bottom[ci] = distances[ci]
    tet.cell_data["cell_distance_to_bottom"] = dist_to_bottom

    for ci in range(tet.number_of_cells):
        if np.isnan(dist_to_bottom[ci]):
            continue
        local = np.hstack((topo["edge"][ci], ci))
        lengths = np.array([dist_to_bottom[c] for c in local])
        local = np.array(local)[~np.isnan(lengths)]
        lengths = lengths[~np.isnan(lengths)]
        if len(lengths) < 3:
            target = tet.cell_data["cell_center"][closest_bottom[ci], :2]
            direction = target - tet.cell_data["cell_center"][ci, :2]
            direction /= np.linalg.norm(direction)
            center = tet.cell_data["cell_center"][ci, :2].copy()
            center /= np.linalg.norm(center)
            d = np.dot(center, direction)
            g = d / np.abs(d) if not np.isnan(d) and d != 0 else 0
            gradient[ci] = 0 if np.isnan(g) else g
        else:
            pts = np.hstack((tet.cell_data["cell_center"][local, :2],
                             lengths[:, None]))
            _, plane_normal = plane_fit(pts.T)
            cdir = tet.cell_data["cell_center"][ci, :2]
            cdir = cdir / np.linalg.norm(cdir)
            g = np.dot(cdir, plane_normal[:2])
            if np.isnan(g):
                neigh = gradient[local][~np.isnan(gradient[local])]
                g = np.mean(neigh) if len(neigh) else 0
                if np.isnan(g):
                    g = 0
            gradient[ci] = g

    # NOTE: preserved notebook behavior - the loop recomputes from the
    # unsmoothed field every pass, so any smoothing>0 is one effective pass
    if smoothing != 0:
        for _ in range(smoothing):
            smoothed = np.zeros(tet.number_of_cells)
            for ci in range(tet.number_of_cells):
                if gradient[ci] != 0:
                    neighbours = topo["point"][ci]
                    local = neighbours.copy()
                    for n in neighbours:
                        local.extend(topo["point"][n])
                    local = np.array(list(set(local)))
                    local = local[gradient[local] != 0]
                    smoothed[ci] = np.mean(gradient[local])
        gradient = smoothed

    if not set_zero:
        gradient[gradient == 0] = np.nan
    tet.cell_data["path_length_to_base_gradient"] = gradient
    return gradient


def calculate_rotation_matrices(tet, rotation_field):
    from scipy.spatial.transform import Rotation as R
    cc2 = tet.cell_data["cell_center"][:, :2]
    tangential = np.cross(np.array([0, 0, 1]),
                          np.column_stack([cc2, np.zeros(len(cc2))]))
    tangential /= np.linalg.norm(tangential, axis=1)[:, None]
    tangential[np.isnan(tangential).any(axis=1)] = [1, 0, 0]
    return R.from_rotvec(rotation_field[:, None] * tangential).as_matrix()


def optimize_rotations(tet, topo, graph, bottom_cells, p):
    from scipy.optimize import least_squares
    from scipy.sparse import csr_matrix

    initial = np.abs(np.deg2rad(90 + p.max_overhang) - tet.cell_data["overhang_angle"])
    gradient = calculate_path_length_to_base_gradient(
        tet, topo, graph, bottom_cells, p.max_overhang,
        p.field_smoothing, p.set_initial_rotation_to_zero)
    if p.steep_overhang_compensation:
        in_air = tet.cell_data["in_air"]
        initial[in_air] += 2 * (np.deg2rad(180) - tet.cell_data["overhang_angle"][in_air])
    initial *= gradient
    initial = np.clip(initial * p.rotation_multiplier,
                      -np.deg2rad(360), np.deg2rad(360))
    rot_limit = np.deg2rad(p.rotation_limit)
    initial = np.clip(initial, -rot_limit, rot_limit)
    tet.cell_data["initial_rotation_field"] = initial

    valid = np.where(~np.isnan(initial))[0]
    n_valid = len(valid)
    neighbours = tet.field_data["cell_face_neighbours"]

    def objective(field):
        diffs = field[neighbours[:, 0]] - field[neighbours[:, 1]]
        return np.concatenate((p.neighbour_loss_weight * diffs**2,
                               (field[valid] - initial[valid])**2))

    n_nb = len(neighbours)
    jac_rows = np.concatenate([np.arange(n_nb), np.arange(n_nb),
                               n_nb + np.arange(n_valid)])
    jac_cols = np.concatenate([neighbours[:, 0], neighbours[:, 1], valid])
    jac_shape = (n_nb + n_valid, tet.number_of_cells)

    def jacobian(field):
        diffs = field[neighbours[:, 0]] - field[neighbours[:, 1]]
        data = np.concatenate([2 * p.neighbour_loss_weight * diffs,
                               -2 * p.neighbour_loss_weight * diffs,
                               2 * (field[valid] - initial[valid])])
        return csr_matrix((data.astype(np.float32), (jac_rows, jac_cols)),
                          shape=jac_shape)

    def sparsity():
        return csr_matrix((np.ones(len(jac_rows), dtype=np.int8),
                           (jac_rows, jac_cols)), shape=jac_shape)

    result = least_squares(objective, np.zeros(tet.number_of_cells),
                           jac=jacobian, max_nfev=p.rot_iterations,
                           jac_sparsity=sparsity(), verbose=1,
                           method="trf", ftol=1e-6)
    # the optimizer can overshoot the initial-field clip; enforce the limit
    return np.clip(result.x, -rot_limit, rot_limit)


def calculate_deformation(tet, rotation_field, iterations):
    from scipy.optimize import least_squares
    from scipy.sparse import csr_matrix

    N = np.eye(4) - 1 / 4 * np.ones((4, 4))
    params = tet.points.copy().flatten()
    rotation_matrices = calculate_rotation_matrices(tet, rotation_field)
    old_vertices = tet.field_data["cell_vertices"][tet.field_data["cells"]]
    old_transformed = np.einsum("ijk,ikl->ijl", rotation_matrices,
                                (N @ old_vertices).transpose(0, 2, 1))

    cells = tet.field_data["cells"]
    cell_idx = np.repeat(np.arange(tet.number_of_cells), cells.shape[1])
    vert_idx = np.ravel(cells)

    def objective(params):
        new_vertices = params.reshape(-1, 3)
        new_transformed = (N @ new_vertices[cells]).transpose(0, 2, 1)
        return np.linalg.norm(new_transformed - old_transformed, axis=(1, 2))**2

    jac_rows = np.tile(cell_idx, 3)
    jac_cols = np.concatenate([vert_idx * 3 + dim for dim in range(3)])
    jac_shape = (tet.number_of_cells, len(params))

    def jacobian(params):
        new_vertices = params.reshape(-1, 3)
        diff = ((N @ new_vertices[cells]).transpose(0, 2, 1)
                - old_transformed).transpose(0, 2, 1)
        data = np.concatenate([2 * diff[:, :, dim].ravel() for dim in range(3)])
        return csr_matrix((data.astype(np.float32), (jac_rows, jac_cols)),
                          shape=jac_shape)

    def sparsity():
        return csr_matrix((np.ones(len(jac_rows), dtype=np.int8),
                           (jac_rows, jac_cols)), shape=jac_shape)

    result = least_squares(objective, params, max_nfev=iterations, verbose=1,
                           jac=jacobian, jac_sparsity=sparsity(),
                           method="trf", x_scale="jac")
    return result.x.reshape(-1, 3)


def deform_stage(cfg, paths):
    import pyvista as pv

    t0 = time.time()
    print(f"[deform] loading + tetrahedralizing {cfg.model}")
    input_tet = load_and_center_mesh(cfg.model, cfg.max_tet_volume, cfg.decimate,
                                     remesh=not cfg.no_remesh)
    print(f"[deform] {input_tet.n_cells} tets; building topology")
    topo, graph = build_topology(input_tet)
    input_tet, _, bottom_cells = calculate_tet_attributes(input_tet, graph)

    undeformed = input_tet.copy()
    deformed = None
    for pass_no in range(cfg.passes):
        print(f"[deform] pass {pass_no + 1}/{cfg.passes}: rotation field "
              f"(limit +/-{cfg.rotation_limit:g} deg)")
        field = optimize_rotations(undeformed, topo, graph, bottom_cells, cfg)
        print(f"[deform] pass {pass_no + 1}: deformation "
              f"({cfg.deform_iterations} max iterations)")
        new_vertices = calculate_deformation(undeformed, field, cfg.deform_iterations)
        deformed = pv.UnstructuredGrid(
            undeformed.cells,
            np.full(undeformed.number_of_cells, pv.CellType.TETRA), new_vertices)
        for key in undeformed.field_data.keys():
            deformed.field_data[key] = undeformed.field_data[key]
        for key in undeformed.cell_data.keys():
            deformed.cell_data[key] = undeformed.cell_data[key]
        deformed = update_tet_attributes(deformed, graph)
        undeformed = deformed.copy()

    x_min, x_max, y_min, y_max, z_min, z_max = deformed.bounds
    deformed.points -= np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])
    # refresh the copies the inverse map reads for barycentric tests - they
    # went stale when the points were recentered (the notebook re-ran
    # calculate_tet_attributes after unpickling for the same reason)
    deformed.field_data["cell_vertices"] = deformed.points
    deformed.cell_data["cell_center"] = deformed.cell_centers().points

    deformed.extract_surface().save(paths["stl"])
    with open(paths["bundle"], "wb") as fh:
        pickle.dump({"input_tet": input_tet, "deformed_tet": deformed}, fh)
    print(f"[deform] done in {time.time() - t0:.0f}s -> {paths['stl']}")


# --------------------------------------------------------------------------
# Stage 2: headless CuraEngine slice with resolved profile settings
# --------------------------------------------------------------------------

class CuraSettingResolver:
    """Resolve the setting stack inside a Cura project 3mf the way the Cura
    frontend does: definition defaults + '=expression' values, overlaid by
    definition_changes < quality < quality_changes < user profiles."""

    PRECEDENCE = {"definition_changes": 0, "quality": 1, "intent": 2,
                  "quality_changes": 3, "user": 4}

    def __init__(self, threemf_path):
        self.defaults = {}     # name -> default_value
        self.exprs = {}        # name -> '=...' expression from definitions
        self.overrides = {}    # name -> literal or '=expr' from the stack
        self.cache = {}
        self.resolving = set()
        self.failed = []
        self._load(threemf_path)

    def _walk_settings(self, node):
        for name, spec in node.items():
            if not isinstance(spec, dict):
                continue
            if "default_value" in spec:
                self.defaults[name] = spec["default_value"]
            if "value" in spec:
                self.exprs[name] = spec["value"]
            if "children" in spec:
                self._walk_settings(spec["children"])

    def _load(self, threemf_path):
        tmp = tempfile.mkdtemp(prefix="cura3mf_")
        with zipfile.ZipFile(threemf_path) as zf:
            zf.extractall(tmp)
        cura_dir = os.path.join(tmp, "Cura")

        for def_name in ("custom.def.json", "custom_extruder_1.def.json"):
            path = os.path.join(cura_dir, def_name)
            if os.path.exists(path):
                d = json.load(open(path))
                self._walk_settings(d.get("settings", {}))
                for name, spec in d.get("overrides", {}).items():
                    if "default_value" in spec:
                        self.defaults[name] = spec["default_value"]
                    if "value" in spec:
                        self.exprs[name] = spec["value"]

        layers = []
        for fn in os.listdir(cura_dir):
            if not fn.endswith(".inst.cfg"):
                continue
            cp = configparser.ConfigParser(interpolation=None)
            cp.read(os.path.join(cura_dir, fn))
            ctype = cp.get("metadata", "type", fallback="")
            if ctype in self.PRECEDENCE and cp.has_section("values"):
                layers.append((self.PRECEDENCE[ctype], fn, dict(cp["values"])))
        for _, _, values in sorted(layers, key=lambda x: x[0]):
            self.overrides.update(values)
        shutil.rmtree(tmp, ignore_errors=True)

    @staticmethod
    def _parse_literal(text):
        try:
            return ast.literal_eval(text)
        except (ValueError, SyntaxError):
            lowered = text.strip().lower()
            if lowered in ("true", "false"):
                return lowered == "true"
            return text

    def set_override(self, name, value):
        self.overrides[name] = value
        self.cache.pop(name, None)

    def resolve(self, name):
        if name in self.cache:
            return self.cache[name]
        if name in self.resolving:  # cycle: fall back to default
            return self.defaults.get(name)
        self.resolving.add(name)
        try:
            if name in self.overrides:
                raw = self.overrides[name]
                if isinstance(raw, str) and raw.startswith("="):
                    value = self._eval(name, raw[1:])
                else:
                    value = self._parse_literal(raw) if isinstance(raw, str) else raw
            elif name in self.exprs:
                raw = self.exprs[name]
                if isinstance(raw, str) and raw.startswith("="):
                    raw = raw[1:]
                value = (self._eval(name, raw)
                         if isinstance(raw, str) else raw)
            elif name in self.defaults:
                value = self.defaults[name]
            else:
                value = None
        finally:
            self.resolving.discard(name)
        self.cache[name] = value
        return value

    def _eval(self, name, expr):
        resolver = self

        class Names(dict):
            def __missing__(self, key):
                v = resolver.resolve(key)
                if v is None:
                    raise KeyError(key)
                return v

        env = {
            "__builtins__": {"min": min, "max": max, "abs": abs, "round": round,
                             "int": int, "float": float, "bool": bool, "str": str,
                             "len": len, "sum": sum, "any": any, "all": all,
                             "map": map, "sorted": sorted, "True": True,
                             "False": False, "None": None},
            "math": math,
            "extruderValue": lambda n, key: resolver.resolve(key),
            "extruderValues": lambda key: [resolver.resolve(key)],
            "resolveOrValue": lambda key: resolver.resolve(key),
            "defaultExtruderPosition": lambda: "0",
            "anyExtruderWithMaterial": lambda *a: "0",
            "anyExtruderNrWithOrDefault": lambda *a: 0,
            "valueFromContainer": lambda key, idx: resolver.resolve(key),
            "extruderPosition": lambda *a: 0,
        }
        try:
            return eval(expr, env, Names())  # noqa: S307 - trusted profile
        except Exception:
            self.failed.append(name)
            return self.defaults.get(name)

    def all_resolved(self):
        names = set(self.defaults) | set(self.exprs) | set(self.overrides)
        out = {}
        for name in sorted(names):
            v = self.resolve(name)
            if v is not None:
                out[name] = v
        return out


def format_setting(value):
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, (list, dict)):
        return json.dumps(value)
    return str(value)


def slice_stage(cfg, paths):
    t0 = time.time()
    resolver = CuraSettingResolver(cfg.cura_config)

    # kill the start-gcode prime: with relative extrusion + center-is-zero it
    # extrudes a blob in mid-air that maps into the print space (this was the
    # phantom z~19 blob in the original prints)
    if not cfg.keep_start_gcode:
        resolver.set_override("machine_start_gcode", "G28 ; home")
    overrides = {k: v for k, v in
                 cfg.printer_config.get("cura_overrides", {}).items()
                 if not k.startswith("_")}
    for key, value in overrides.items():
        resolver.set_override(key, value)
    if overrides:
        print(f"[slice] applied {len(overrides)} Cura overrides "
              f"from {cfg.printer}")
    if cfg.layer_height is not None:
        resolver.set_override("layer_height", str(cfg.layer_height))
        resolver.set_override("layer_height_0", str(cfg.layer_height))
    if cfg.retraction_amount is not None:
        resolver.set_override("retraction_amount", str(cfg.retraction_amount))
    if cfg.retraction_speed is not None:
        resolver.set_override("retraction_speed", str(cfg.retraction_speed))
    if cfg.nozzle_size is not None:
        resolver.set_override("machine_nozzle_size", str(cfg.nozzle_size))
    for kv in cfg.cura_set or []:
        key, _, value = kv.partition("=")
        if not value:
            sys.exit(f"--cura-set needs key=value, got {kv!r}")
        resolver.set_override(key.strip(), value.strip())

    settings = resolver.all_resolved()
    if resolver.failed:
        uniq = sorted(set(resolver.failed))
        print(f"[slice] note: {len(uniq)} settings fell back to defaults "
              f"(unresolvable expressions): {', '.join(uniq[:8])}"
              + (" ..." if len(uniq) > 8 else ""))

    # sanity: the inverse mapper depends on these
    if not settings.get("relative_extrusion"):
        sys.exit("profile must use relative_extrusion=True for the inverse map")
    retraction = float(settings.get("retraction_amount", 1.0))

    defs = os.path.join(os.path.dirname(cfg.cura_engine),
                        "share", "cura", "resources", "definitions",
                        "fdmprinter.def.json")
    args = [cfg.cura_engine, "slice", "-j", defs]
    for key, value in settings.items():
        args += ["-s", f"{key}={format_setting(value)}"]
    args += ["-e0"]
    for key, value in settings.items():
        args += ["-s", f"{key}={format_setting(value)}"]
    args += ["-l", paths["stl"], "-o", paths["sliced"]]

    print(f"[slice] CuraEngine with {len(settings)} resolved settings "
          f"(layer {settings.get('layer_height')}, retraction {retraction}, "
          f"line width {settings.get('line_width')})")
    proc = subprocess.run(args, capture_output=True, text=True)
    errors = [ln for ln in proc.stderr.splitlines() if "[error]" in ln]
    if errors:
        print(f"[slice] CuraEngine reported {len(errors)} errors, e.g.:")
        for ln in errors[:5]:
            print("   ", ln)
    if not os.path.exists(paths["sliced"]) or os.path.getsize(paths["sliced"]) == 0:
        print(proc.stderr[-3000:])
        sys.exit("[slice] CuraEngine produced no output")
    size = os.path.getsize(paths["sliced"])
    print(f"[slice] done in {time.time() - t0:.0f}s -> "
          f"{paths['sliced']} ({size // 1024} KB)")
    return retraction


# --------------------------------------------------------------------------
# Stage 3: inverse map planar gcode -> S4 4-axis gcode (cells 15-18)
# --------------------------------------------------------------------------

def tetrahedron_volume(p1, p2, p3, p4):
    return np.abs(np.linalg.det(np.vstack([p2 - p1, p3 - p1, p4 - p1]))) / 6


def calc_barycentric(a, b, c, d, point):
    total = tetrahedron_volume(a, b, c, d)
    if total == 0:
        raise ValueError("zero-volume tetrahedron")
    return np.array([tetrahedron_volume(point, b, c, d),
                     tetrahedron_volume(point, a, c, d),
                     tetrahedron_volume(point, a, b, d),
                     tetrahedron_volume(point, a, b, c)]) / total


def project_onto_plane(x_axis, y_axis, points):
    return np.array([np.sum(x_axis * points, axis=1),
                     np.sum(y_axis * points, axis=1)]).T


def invmap_stage(cfg, paths, retraction_length):
    from pygcode import Line

    t0 = time.time()
    with open(paths["bundle"], "rb") as fh:
        bundle = pickle.load(fh)
    input_tet = bundle["input_tet"]
    deformed_tet = bundle["deformed_tet"]

    max_rot = (cfg.s4_max_rotation if cfg.s4_max_rotation is not None
               else cfg.rotation_limit)
    min_rot = (cfg.s4_min_rotation if cfg.s4_min_rotation is not None
               else -cfg.rotation_limit)
    print(f"[invmap] S4 B clamp [{min_rot:g}, {max_rot:g}] deg, "
          f"seg {cfg.seg_size} mm, retraction {retraction_length} mm")

    vertex_transformations = deformed_tet.points - input_tet.points

    cells_arr = np.asarray(deformed_tet.field_data["cells"])
    num_cells_per_vertex = np.zeros(input_tet.number_of_points)
    np.add.at(num_cells_per_vertex, cells_arr.ravel(), 1)

    # per-cell Kabsch rotation in the radial plane, batched
    new_v = deformed_tet.field_data["cell_vertices"][cells_arr] \
        - deformed_tet.cell_data["cell_center"][:, None, :]
    old_v = np.asarray(input_tet.field_data["cell_vertices"])[cells_arr] \
        - input_tet.cell_data["cell_center"][:, None, :]
    px = np.zeros((len(cells_arr), 3))
    cc2 = input_tet.cell_data["cell_center"][:, :2]
    nrm = np.linalg.norm(cc2, axis=1)
    ok = nrm > 0
    px[ok, :2] = cc2[ok] / nrm[ok, None]
    px[~ok] = [1.0, 0.0, 0.0]
    new_p = np.stack([np.einsum("cij,cj->ci", new_v, px), new_v[:, :, 2]], axis=-1)
    old_p = np.stack([np.einsum("cij,cj->ci", old_v, px), old_v[:, :, 2]], axis=-1)
    cov = np.einsum("cik,cil->ckl", new_p, old_p)
    U, _, Vt = np.linalg.svd(cov)
    Rm = U @ Vt
    cell_rotations = -np.arccos(np.clip(Rm[:, 0, 0], -1.0, 1.0))
    cell_rotations = np.where(Rm[:, 1, 0] < 0, -cell_rotations, cell_rotations)
    cell_rotations = np.clip(cell_rotations,
                             np.deg2rad(min_rot), np.deg2rad(max_rot))
    vertex_rotations = np.zeros(deformed_tet.number_of_points)
    np.add.at(vertex_rotations, cells_arr,
              cell_rotations[:, None] / num_cells_per_vertex[cells_arr])

    def batched_volume(p):  # p: (N, 4, 3)
        m = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]],
                     axis=1)
        return np.abs(np.linalg.det(m)) / 6

    z_squish = (batched_volume(np.asarray(input_tet.field_data["cell_vertices"])[cells_arr])
                / batched_volume(deformed_tet.field_data["cell_vertices"][cells_arr]))

    # read + segment the planar gcode
    pos = np.array([0., 0., 20.])
    feed = 5000.0
    gcode_points = []
    with open(paths["sliced"]) as fh:
        for line_text in fh:
            line = Line(line_text)
            if not line.block.gcodes:
                continue
            for gcode in sorted(line.block.gcodes):
                if gcode.word not in ("G00", "G01"):
                    continue
                prev_pos = pos.copy()
                if gcode.X is not None:
                    pos[0] = gcode.X
                if gcode.Y is not None:
                    pos[1] = gcode.Y
                if gcode.Z is not None:
                    pos[2] = gcode.Z
                for word in line.block.words:
                    if word.letter == "F":
                        feed = word.value
                extrusion = None
                for param in line.block.modal_params:
                    if param.letter == "E":
                        extrusion = param.value
                delta = pos - prev_pos
                distance = np.linalg.norm(delta)
                if distance > 0:
                    n_seg = int(-(-distance // cfg.seg_size))
                    seg_len = distance / n_seg
                    t_move = (1 / feed) * seg_len
                    itf = None if t_move == 0 else 1 / t_move
                    for i in range(n_seg):
                        gcode_points.append({
                            "position": prev_pos + delta * (i + 1) / n_seg,
                            "command": gcode.word,
                            "extrusion": extrusion / n_seg if extrusion is not None else None,
                            "inv_time_feed": itf, "feed": feed})
                else:
                    gcode_points.append({
                        "position": pos.copy(), "command": gcode.word,
                        "extrusion": extrusion, "inv_time_feed": None,
                        "feed": feed})

    print(f"[invmap] {len(gcode_points)} segmented points; locating cells")
    positions = [p["position"] for p in gcode_points]
    containing = deformed_tet.find_containing_cell(positions)
    closest = deformed_tet.find_closest_cell(positions)

    # batched barycentric interpolation for every point (same math as the
    # notebook's per-point calc_barycentric, incl. the sum<=1.01 gate)
    pos_arr = np.asarray(positions)
    cont_arr = np.asarray(containing)
    is_print = np.array([p["command"] == "G01" for p in gcode_points])
    eff = cont_arr.copy()
    fallback = (eff == -1) & is_print
    eff[fallback] = np.asarray(closest)[fallback]
    pt_valid = eff != -1
    vidx_all = cells_arr[np.clip(eff, 0, None)]
    cv = deformed_tet.field_data["cell_vertices"][vidx_all]  # (P,4,3)
    with np.errstate(divide="ignore", invalid="ignore"):
        total = batched_volume(cv)
        bary = np.stack([
            batched_volume(np.stack([pos_arr, cv[:, 1], cv[:, 2], cv[:, 3]], axis=1)),
            batched_volume(np.stack([pos_arr, cv[:, 0], cv[:, 2], cv[:, 3]], axis=1)),
            batched_volume(np.stack([pos_arr, cv[:, 0], cv[:, 1], cv[:, 3]], axis=1)),
            batched_volume(np.stack([pos_arr, cv[:, 0], cv[:, 1], cv[:, 2]], axis=1)),
        ], axis=1) / total[:, None]
    with np.errstate(invalid="ignore"):
        pt_valid &= ~(np.nansum(bary, axis=1) > 1.01) & ~np.isnan(bary).any(axis=1)
    transform_all = np.einsum("pij,pi->pj", vertex_transformations[vidx_all], bary)
    new_pos_all = pos_arr - transform_all
    rot_all = np.einsum("pi,pi->p", vertex_rotations[vidx_all], bary)

    ROTATION_ALPHA = 0.2
    ROTATION_MAX_DELTA = np.deg2rad(1)
    MAX_EXTRUSION_MULT = 10
    new_points = []
    prev_new_position = None
    prev_rotation = 0
    prev_command = "G00"
    prev_travelling = False
    travelling = False
    travelling_over_air = False
    highest_printed = 0
    lost = 0

    for gi, (point, containing_ci) in enumerate(zip(gcode_points, containing)):
        position = point["position"]
        command = point["command"]
        extrusion = point["extrusion"]
        dont_smooth = False
        if pt_valid[gi]:
            new_position, rotation = new_pos_all[gi].copy(), float(rot_all[gi])
        else:
            new_position, rotation = None, None
        if new_position is None:
            if command == "G01":
                lost += 1
                continue
            elif command == "G00" and not travelling_over_air \
                    and prev_new_position is not None:
                new_position = np.array([prev_new_position[0],
                                         prev_new_position[1], highest_printed])
                rotation = max(min(prev_rotation, np.deg2rad(45)), np.deg2rad(-45))
                dont_smooth = True
                travelling_over_air = True
            else:
                continue
        else:
            if travelling_over_air:
                new_position[2] = highest_printed
                rotation = max(min(rotation, np.deg2rad(45)), np.deg2rad(-45))
                dont_smooth = True
            travelling_over_air = False

        mult = 1
        if extrusion is not None and extrusion != retraction_length \
                and extrusion != -retraction_length:
            mult *= z_squish[containing_ci]
            extrusion = extrusion * min(mult, MAX_EXTRUSION_MULT)
        elif extrusion == -retraction_length:
            travelling = True
        elif extrusion == retraction_length:
            travelling = False
        if prev_rotation is not None and not dont_smooth:
            rotation = ROTATION_ALPHA * rotation + (1 - ROTATION_ALPHA) * prev_rotation

        if prev_rotation is not None and prev_new_position is not None \
                and np.abs(rotation - prev_rotation) > ROTATION_MAX_DELTA:
            delta_rot = rotation - prev_rotation
            n_interp = int(np.abs(delta_rot) / ROTATION_MAX_DELTA) + 1
            delta_pos = new_position - prev_new_position
            for i in range(n_interp):
                new_points.append({
                    "position": prev_new_position + delta_pos * (i + 1) / n_interp,
                    "rotation": prev_rotation + delta_rot * (i + 1) / n_interp,
                    "command": prev_command,
                    "extrusion": extrusion / n_interp if extrusion is not None else None,
                    "inv_time_feed": point["inv_time_feed"] * n_interp
                        if point["inv_time_feed"] is not None else None,
                    "travelling": prev_travelling})
        else:
            new_points.append({
                "position": new_position, "rotation": rotation,
                "command": command, "extrusion": extrusion,
                "inv_time_feed": point["inv_time_feed"],
                "travelling": travelling})

        prev_rotation = rotation
        prev_new_position = new_position.copy()
        prev_travelling = travelling
        prev_command = command
        if command == "G01" and extrusion is not None and extrusion > 0 \
                and (highest_printed != 0 or new_position[2] < 1):
            highest_printed = max(highest_printed, new_position[2])

    print(f"[invmap] lost {lost} print points outside the mesh")

    # emit S4 gcode (identical format to the notebook)
    prev_theta = 0.0
    theta_accum = 0.0
    with open(paths["s4"], "w") as fh:
        fh.write("G94 ; mm/min feed  \n")
        fh.write("G28 ; home \n")
        fh.write("M83 ; relative extrusion \n")
        fh.write("G1 E10 ; prime extruder \n")
        fh.write("G94 ; mm/min feed \n")
        fh.write("G90 ; absolute positioning \n")
        fh.write("G0 C0 X0 Z20 B0 ; go to start \n")
        fh.write("G93 ; inverse time feed \n")
        for point in new_points:
            position = point["position"]
            rotation = point["rotation"]
            if np.all(np.isnan(position)) or position[2] < 0:
                continue
            z_hop = 1 if point["travelling"] else 0
            r = np.linalg.norm(position[:2])
            theta = np.arctan2(position[1], position[0])
            z = position[2]
            r += -np.sin(rotation) * (cfg.nozzle_offset + z_hop)
            z += (np.cos(rotation) - 1) * (cfg.nozzle_offset + z_hop) + z_hop
            delta_theta = theta - prev_theta
            if delta_theta > np.pi:
                delta_theta -= 2 * np.pi
            if delta_theta < -np.pi:
                delta_theta += 2 * np.pi
            theta_accum += delta_theta
            string = (f"{point['command']} C{np.rad2deg(theta_accum):.5f} "
                      f"X{r:.5f} Z{z:.5f} B{np.rad2deg(rotation):.5f}")
            if point["extrusion"] is not None:
                string += f" E{point['extrusion']:.4f}"
            if point["inv_time_feed"] is not None:
                string += f" F{point['inv_time_feed']:.4f}"
                fh.write(string + "\n")
            else:
                fh.write("G94\n" + string + " F20000\n" + "G93\n")
            prev_theta = theta
    print(f"[invmap] done in {time.time() - t0:.0f}s -> {paths['s4']}")


# --------------------------------------------------------------------------
# Stages 4+5: stewart conversion + collision check
# --------------------------------------------------------------------------

def stewart_stage(cfg, paths):
    argv = [paths["s4"], "-o", paths["stewart"], "-m", cfg.mode,
            "-a", str(cfg.max_tilt), "-c", str(cfg.max_clearance),
            "--nozzle-offset", str(cfg.nozzle_offset),
            "--travel-hop", "1.0",
            "--out-travel-hop", str(cfg.out_travel_hop),
            "--bed-safe-z", str(cfg.bed_safe_z),
            "--tilt-ramp", str(cfg.tilt_ramp),
            "--bed-phase", cfg.bed_phase,
            "--speed-scale", str(cfg.speed_scale),
            "--temp", str(cfg.temp)]
    if cfg.print_cap is not None:
        argv += ["--print-cap", str(cfg.print_cap)]
    if os.path.exists(cfg.printer):
        argv += ["--printer", cfg.printer]
    if cfg.kinematics:
        argv += ["--kinematics", cfg.kinematics]
    if cfg.start_gcode:
        argv += ["--start-gcode", cfg.start_gcode]
    if cfg.end_gcode:
        argv += ["--end-gcode", cfg.end_gcode]
    print(f"[stewart] s4_to_stewart {' '.join(argv)}")
    s4_to_stewart.main(argv)


def check_stage(cfg, paths):
    args = [sys.executable, os.path.join(HERE, "check_nozzle_collisions.py"),
            paths["stewart"], "-c", str(cfg.max_clearance),
            "--check-every", str(cfg.check_every)]
    print(f"[check] scanning {paths['stewart']} (this takes a few minutes)")
    proc = subprocess.run(args, capture_output=True, text=True)
    report = proc.stdout + proc.stderr
    with open(paths["report"], "w") as fh:
        fh.write(report)
    print(report)
    print(f"[check] report saved to {paths['report']}")


# --------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="S4 deform -> CuraEngine -> inverse map -> Stewart pose "
                    "gcode, in one command.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("model", nargs="?", default=None, help="input STL")
    p.add_argument("--model", dest="model_opt", default=None,
                   help="input STL (same as the positional argument)")
    p.add_argument("-n", "--name", default=None, help="job name (default: STL stem)")
    p.add_argument("--outdir", default=None,
                   help="output directory (default: pipeline_out/<name>)")

    g = p.add_argument_group("machine limits")
    g.add_argument("-a", "--max-tilt", type=float, default=20.0, metavar="DEG")
    g.add_argument("-c", "--max-clearance", type=float, default=25.0, metavar="DEG")
    g.add_argument("--rotation-limit", default="auto", metavar="DEG|auto|off",
                   help="clamp on layer tilt in deformation AND S4 B output. "
                        "Default 'auto' = max_tilt + max_clearance; give "
                        "degrees for an explicit clamp, or 'off' for "
                        "unclamped deformation (like the original notebook)")
    g.add_argument("--no-rotation-limit", action="store_true",
                   help=argparse.SUPPRESS)  # now the default; kept for old cmds
    g.add_argument("-m", "--mode", choices=("3axis", "tilt", "tilt-signed"),
                   default="tilt-signed")

    g = p.add_argument_group("deformation (stage 1)")
    g.add_argument("--max-overhang", type=float, default=30.0, metavar="DEG")
    g.add_argument("--max-tet-volume", type=float, default=None, metavar="MM3",
                   help="densify the tet mesh to this max tet volume; needed "
                        "for low-poly STLs (target a few thousand tets)")
    g.add_argument("--no-remesh", action="store_true",
                   help="skip the pymeshlab isotropic remesh that runs first "
                        "on every STL (evenly sized triangles -> well-shaped "
                        "tets). Target triangle count = --decimate, or the "
                        "STL's own count when --decimate is not given.")
    g.add_argument("--decimate", type=int, default=None, metavar="TRIS",
                   help="target triangle count for the remesh (or plain "
                        "quadric decimation with --no-remesh) before "
                        "tetrahedralizing (the deform stage cost grows with "
                        "tet count; ~8000 is plenty). Off by default: the "
                        "STL's full resolution is kept.")
    g.add_argument("--passes", type=int, default=1,
                   help="deformation passes (notebook 'run another iteration')")
    g.add_argument("--rot-iterations", type=int, default=100)
    g.add_argument("--deform-iterations", type=int, default=1000)
    g.add_argument("--neighbour-loss-weight", type=float, default=20.0)
    g.add_argument("--rotation-multiplier", type=float, default=2.0)
    g.add_argument("--field-smoothing", type=int, default=30)
    g.add_argument("--no-steep-overhang-compensation", dest="steep_overhang_compensation",
                   action="store_false")
    g.add_argument("--set-initial-rotation-to-zero", action="store_true")

    g = p.add_argument_group("slicing (stage 2)")
    g.add_argument("--cura-engine", default=DEFAULT_CURA_ENGINE,
                   help="CuraEngine binary (default: $CURA_ENGINE, then the "
                        "macOS UltiMaker Cura app bundle, then PATH)")
    g.add_argument("--cura-config", default=DEFAULT_CURA_CONFIG,
                   help="Cura project 3mf holding the profile stack")
    g.add_argument("--layer-height", type=float, default=None,
                   help="override profile layer height")
    g.add_argument("--retraction-amount", type=float, default=4.0, metavar="MM",
                   help="retraction (bowden default 4; the inverse map is "
                        "kept in sync automatically)")
    g.add_argument("--retraction-speed", type=float, default=None, metavar="MM/S")
    g.add_argument("--nozzle-size", type=float, default=None, metavar="MM")
    g.add_argument("--cura-set", action="append", metavar="KEY=VALUE",
                   help="arbitrary Cura setting override (repeatable), e.g. "
                        "--cura-set speed_print=50 --cura-set infill_sparse_density=15")
    g.add_argument("--keep-start-gcode", action="store_true",
                   help="keep the profile's start gcode incl. the E3 prime "
                        "blob (default: replaced with plain G28)")

    g = p.add_argument_group("inverse map (stage 3)")
    g.add_argument("--seg-size", type=float, default=0.6, metavar="MM")
    g.add_argument("--nozzle-offset", type=float, default=42.0, metavar="MM")
    g.add_argument("--s4-max-rotation", type=float, default=None, metavar="DEG",
                   help="override the +B clamp (default: rotation-limit)")
    g.add_argument("--s4-min-rotation", type=float, default=None, metavar="DEG",
                   help="override the -B clamp (default: -rotation-limit)")

    g = p.add_argument_group("printer config & customization")
    g.add_argument("--printer", default=os.path.join(HERE, "printer.json"),
                   metavar="JSON",
                   help="printer config: 'machine' section = defaults for "
                        "max tilt / clearance / nozzle offset / temp / ... "
                        "(CLI flags still win), 'start_gcode'/'end_gcode' "
                        "lists = custom start/end G-code, 'cura_overrides' "
                        "= Cura settings for your machine")
    g.add_argument("--kinematics", default=None, metavar="PY",
                   help="python file defining get_kinematics() to emit "
                        "G-code for a non-Stewart 5-axis machine (see "
                        "kinematics_example.py)")
    g.add_argument("--start-gcode", default=None, metavar="FILE",
                   help="file of start G-code lines ({temp} placeholder); "
                        "wins over printer.json")
    g.add_argument("--end-gcode", default=None, metavar="FILE",
                   help="file of end G-code lines ({lift_z} placeholder); "
                        "wins over printer.json")

    g = p.add_argument_group("stewart conversion (stage 4)")
    g.add_argument("--out-travel-hop", type=float, default=2.0, metavar="MM")
    g.add_argument("--bed-safe-z", type=float, default=1.0, metavar="MM",
                   help="below this tip height the nozzle stays vertical")
    g.add_argument("--tilt-ramp", type=float, default=0.5, metavar="MM",
                   help="band above --bed-safe-z over which tilt fades back in")
    g.add_argument("--bed-phase", choices=("auto", "off"), default="off",
                   help="off: per-move z ramp only; auto: whole print level "
                        "until the last bed-adjacent extrusion (original "
                        "behavior, for models that only touch the bed early)")
    g.add_argument("--speed-scale", type=float, default=1.0)
    g.add_argument("--print-cap", type=float, default=10.0, metavar="MM/S")
    g.add_argument("--temp", type=float, default=200.0, metavar="C")

    g = p.add_argument_group("collision check (stage 5)")
    g.add_argument("--check", action="store_true",
                   help="run the collision scan (off by default - it dominates "
                        "run time; run it once on a final candidate before "
                        "printing). The converter's residual stats always print.")
    g.add_argument("--check-every", type=int, default=4)

    g = p.add_argument_group("flow control")
    g.add_argument("--from-gcode", default=None, metavar="GCODE",
                   help="skip stages 1-2: use this hand-sliced gcode of the "
                        "deformed STL (needs the bundle.pkl of a previous run "
                        "in --outdir)")
    g.add_argument("--skip-deform", action="store_true",
                   help="reuse deformed.stl + bundle.pkl already in --outdir")
    g.add_argument("--stop-after", choices=("deform", "slice", "invmap", "stewart"),
                   default=None)
    return p


def main(argv=None):
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--printer", default=os.path.join(HERE, "printer.json"))
    pre_cfg, _ = pre.parse_known_args(argv)
    parser = build_parser()
    printer_config = s4_to_stewart.load_printer_config(
        pre_cfg.printer, required=pre_cfg.printer
                                  != os.path.join(HERE, "printer.json"))
    machine = dict(printer_config.get("machine", {}))
    # options the stewart converter knows but the pipeline does not
    for key in ("tilt_epsilon", "clearance_policy", "feed_cap",
                "feed_fallback", "travel_hop", "no_fan", "fan_speed"):
        machine.pop(key, None)
    s4_to_stewart.apply_machine_defaults(parser, machine, strict=True)
    cfg = parser.parse_args(argv)
    cfg.printer_config = printer_config
    cfg.model = cfg.model or cfg.model_opt
    if cfg.model is None and sys.stdin.isatty():
        while cfg.model is None:
            try:
                raw = input("model STL path (drag & drop into the terminal works): ").strip()
            except (EOFError, KeyboardInterrupt):
                sys.exit("\nno input model given")
            if not raw:
                sys.exit("no input model given")
            raw = raw.strip("'\"")
            # a path dragged into the terminal arrives backslash-escaped
            for cand in (raw, re.sub(r"\\(.)", r"\1", raw)):
                cand = os.path.expanduser(cand)
                if os.path.exists(cand):
                    cfg.model = cand
                    break
            else:
                print(f"  not found: {raw}")
    if cfg.model is None:
        sys.exit("no input STL given (positional argument or --model)")
    if cfg.name is None:
        cfg.name = os.path.splitext(os.path.basename(cfg.model))[0]
    if (cfg.no_rotation_limit
            or str(cfg.rotation_limit).lower() in ("off", "none", "unlimited")):
        cfg.rotation_limit = 3600.0  # unclamped
    elif str(cfg.rotation_limit).lower() == "auto":
        cfg.rotation_limit = cfg.max_tilt + cfg.max_clearance
    else:
        try:
            cfg.rotation_limit = float(cfg.rotation_limit)
        except ValueError:
            sys.exit(f"--rotation-limit: expected degrees or 'auto', "
                     f"got {cfg.rotation_limit!r}")
    if cfg.outdir is None:
        cfg.outdir = os.path.join(HERE, "pipeline_out", cfg.name)
    os.makedirs(cfg.outdir, exist_ok=True)

    paths = {
        "stl": os.path.join(cfg.outdir, "deformed.stl"),
        "bundle": os.path.join(cfg.outdir, "bundle.pkl"),
        "sliced": os.path.join(cfg.outdir, "sliced_planar.gcode"),
        "s4": os.path.join(cfg.outdir, "s4_4axis.gcode"),
        "stewart": os.path.join(cfg.outdir, f"stewart_{cfg.mode}.gcode"),
        "report": os.path.join(cfg.outdir, "collision_report.txt"),
    }

    print(f"=== s4_pipeline: {cfg.name} -> {cfg.outdir}")
    print(f"=== max_tilt {cfg.max_tilt}  max_clearance {cfg.max_clearance}  "
          f"rotation_limit {cfg.rotation_limit}  mode {cfg.mode}")

    if cfg.from_gcode:
        if not os.path.exists(paths["bundle"]):
            sys.exit(f"--from-gcode needs {paths['bundle']} from a previous "
                     f"deform run of this model")
        shutil.copyfile(cfg.from_gcode, paths["sliced"])
        retraction = cfg.retraction_amount if cfg.retraction_amount else 1.0
    else:
        if cfg.skip_deform:
            if not (os.path.exists(paths["stl"]) and os.path.exists(paths["bundle"])):
                sys.exit("--skip-deform: no deformed.stl/bundle.pkl in outdir")
            print("[deform] skipped (reusing existing artifacts)")
        else:
            deform_stage(cfg, paths)
        if cfg.stop_after == "deform":
            return 0
        retraction = slice_stage(cfg, paths)
        if cfg.stop_after == "slice":
            return 0

    invmap_stage(cfg, paths, retraction)
    if cfg.stop_after == "invmap":
        return 0
    stewart_stage(cfg, paths)
    if cfg.stop_after == "stewart" or not cfg.check:
        return 0
    check_stage(cfg, paths)
    return 0


if __name__ == "__main__":
    sys.exit(main())
