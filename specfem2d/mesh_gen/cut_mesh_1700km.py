#!/usr/bin/env python3
"""
Cut the full-Earth SPECFEM2D mesh to a 1700 km radius domain,
matching the DFDM computational domain.

Creates new mesh files with:
- Elements with centroid radius < R_CUT
- Absorbing boundary at the outer edge (~R_CUT)
- No free surface
- Renumbered nodes and elements
"""

import numpy as np
from collections import defaultdict
import os
import sys

# ============================================================
# Configuration
# ============================================================
R_CUT = 1700.0e3  # meters (1700 km)

MESH_DIR = os.path.join(os.path.dirname(__file__), "MESH")

# Input files (full Earth)
NODES_IN = os.path.join(MESH_DIR, "Nodes_AK135F_NO_MUD")
MESH_IN = os.path.join(MESH_DIR, "Mesh_AK135F_NO_MUD")
MATERIAL_IN = os.path.join(MESH_DIR, "Material_AK135F_NO_MUD")

# Output files (1700 km domain)
SUFFIX = "_1700km"
NODES_OUT = os.path.join(MESH_DIR, f"Nodes{SUFFIX}")
MESH_OUT = os.path.join(MESH_DIR, f"Mesh{SUFFIX}")
MATERIAL_OUT = os.path.join(MESH_DIR, f"Material{SUFFIX}")
SURF_ABS_OUT = os.path.join(MESH_DIR, f"Surf_abs{SUFFIX}")
SURF_FREE_OUT = os.path.join(MESH_DIR, f"Surf_free{SUFFIX}")

# SPECFEM2D 9-node element edge definitions (0-indexed node positions)
# Node ordering:  3 - 6 - 2
#                 |       |
#                 7   8   5
#                 |       |
#                 0 - 4 - 1
EDGES = {
    1: (0, 1),  # bottom: corner nodes 0, 1
    2: (1, 2),  # right:  corner nodes 1, 2
    3: (2, 3),  # top:    corner nodes 2, 3
    4: (3, 0),  # left:   corner nodes 3, 0
}


def load_nodes(filepath):
    """Load node coordinates. Returns (num_nodes, coords[num_nodes, 2])."""
    with open(filepath) as f:
        num_nodes = int(f.readline().strip())
    coords = np.loadtxt(filepath, skiprows=1)
    assert coords.shape == (num_nodes, 2), f"Expected {num_nodes} nodes, got {coords.shape[0]}"
    print(f"  Loaded {num_nodes} nodes")
    return num_nodes, coords


def load_mesh(filepath):
    """Load element connectivity. Returns (num_elem, connectivity[num_elem, 9])."""
    with open(filepath) as f:
        num_elem = int(f.readline().strip())
    conn = np.loadtxt(filepath, skiprows=1, dtype=int)
    assert conn.shape == (num_elem, 9), f"Expected {num_elem}x9, got {conn.shape}"
    print(f"  Loaded {num_elem} elements (9-node quads)")
    return num_elem, conn


def load_materials(filepath):
    """Load material IDs for each element."""
    mat = np.loadtxt(filepath, dtype=int)
    print(f"  Loaded {len(mat)} material assignments")
    return mat


def find_boundary_edges(kept_conn, kept_indices):
    """
    Find edges on the outer boundary of the kept elements.
    An edge is on the boundary if it appears in only one element.

    Returns list of (new_elem_id_1indexed, edge_type, node1, node2).
    """
    # Build edge->element mapping using corner nodes
    edge_to_elems = defaultdict(list)

    for new_idx, conn_row in enumerate(kept_conn):
        for edge_type, (i, j) in EDGES.items():
            n1, n2 = conn_row[i], conn_row[j]
            edge_key = tuple(sorted((n1, n2)))
            edge_to_elems[edge_key].append((new_idx, edge_type, n1, n2))

    # Boundary edges: appear in exactly one element
    boundary_edges = []
    for edge_key, elem_list in edge_to_elems.items():
        if len(elem_list) == 1:
            new_idx, edge_type, n1, n2 = elem_list[0]
            boundary_edges.append((new_idx + 1, edge_type, n1, n2))  # 1-indexed

    return boundary_edges


def main():
    print("=" * 60)
    print(f"Cutting mesh to R < {R_CUT/1e3:.0f} km")
    print("=" * 60)

    # ---- Load input data ----
    print("\nLoading input mesh...")
    num_nodes, coords = load_nodes(NODES_IN)
    num_elem, conn = load_mesh(MESH_IN)
    mat = load_materials(MATERIAL_IN)

    # ---- Select elements within R_CUT ----
    print(f"\nSelecting elements with centroid r < {R_CUT/1e3:.0f} km...")
    centroids = np.zeros((num_elem, 2))
    for i in range(num_elem):
        node_ids = conn[i] - 1  # 0-indexed
        centroids[i] = coords[node_ids].mean(axis=0)

    centroid_r = np.sqrt(centroids[:, 0]**2 + centroids[:, 1]**2)
    keep_mask = centroid_r < R_CUT
    kept_indices = np.where(keep_mask)[0]

    print(f"  Kept {len(kept_indices)} of {num_elem} elements "
          f"({100*len(kept_indices)/num_elem:.1f}%)")

    kept_mat = mat[kept_indices]
    for mid in np.unique(kept_mat):
        print(f"    Material {mid}: {np.sum(kept_mat == mid)} elements")

    # ---- Renumber nodes ----
    print("\nRenumbering nodes...")
    kept_conn = conn[kept_indices]  # still using old node IDs
    used_old_nodes = np.unique(kept_conn.flatten())
    num_new_nodes = len(used_old_nodes)

    # Create mapping: old_node_id (1-indexed) -> new_node_id (1-indexed)
    old_to_new = {}
    for new_id_0, old_id in enumerate(used_old_nodes):
        old_to_new[old_id] = new_id_0 + 1  # 1-indexed

    # Apply remapping to connectivity
    new_conn = np.zeros_like(kept_conn)
    for i in range(len(kept_indices)):
        for j in range(9):
            new_conn[i, j] = old_to_new[kept_conn[i, j]]

    # Extract coordinates for kept nodes (in new order)
    new_coords = coords[used_old_nodes - 1]  # old IDs are 1-indexed

    print(f"  {num_new_nodes} nodes retained (was {num_nodes})")

    # Verify domain extent
    new_r = np.sqrt(new_coords[:, 0]**2 + new_coords[:, 1]**2)
    print(f"  Radius range: [{new_r.min()/1e3:.1f}, {new_r.max()/1e3:.1f}] km")

    # ---- Find absorbing boundary edges ----
    print("\nFinding absorbing boundary edges...")
    boundary_edges = find_boundary_edges(new_conn, kept_indices)

    # Separate outer boundary (r > some threshold) from potential inner edges
    abs_edges = []
    for elem_id, edge_type, n1, n2 in boundary_edges:
        # Check if both nodes are at large radius (outer boundary)
        r1 = np.sqrt(new_coords[n1 - 1, 0]**2 + new_coords[n1 - 1, 1]**2)
        r2 = np.sqrt(new_coords[n2 - 1, 0]**2 + new_coords[n2 - 1, 1]**2)
        avg_r = (r1 + r2) / 2
        if avg_r > 0.5 * R_CUT:  # outer boundary
            abs_edges.append((elem_id, edge_type, n1, n2))

    print(f"  Total boundary edges: {len(boundary_edges)}")
    print(f"  Absorbing boundary edges (outer): {len(abs_edges)}")

    if boundary_edges:
        inner_count = len(boundary_edges) - len(abs_edges)
        if inner_count > 0:
            print(f"  WARNING: {inner_count} inner boundary edges found (central cube?)")

    # ---- Write output files ----
    print("\nWriting output files...")

    # Nodes file
    with open(NODES_OUT, 'w') as f:
        f.write(f"      {num_new_nodes}\n")
        for x, z in new_coords:
            f.write(f"  {x:.16E}  {z:.16E}\n")
    print(f"  {NODES_OUT}")

    # Mesh connectivity file
    num_new_elem = len(kept_indices)
    with open(MESH_OUT, 'w') as f:
        f.write(f"      {num_new_elem}\n")
        for row in new_conn:
            f.write("  " + "  ".join(f"{n:10d}" for n in row) + "\n")
    print(f"  {MESH_OUT}")

    # Material file
    with open(MATERIAL_OUT, 'w') as f:
        for m in kept_mat:
            f.write(f"  {m}\n")
    print(f"  {MATERIAL_OUT}")

    # Absorbing surface file
    with open(SURF_ABS_OUT, 'w') as f:
        f.write(f"        {len(abs_edges)}\n")
        for elem_id, edge_type, n1, n2 in abs_edges:
            f.write(f"       {elem_id:5d}  {edge_type}       {n1:10d}      {n2:10d}\n")
    print(f"  {SURF_ABS_OUT}")

    # Free surface file (empty - no free surface at 1700 km)
    with open(SURF_FREE_OUT, 'w') as f:
        f.write("           0\n")
    print(f"  {SURF_FREE_OUT}")

    # ---- Summary ----
    print("\n" + "=" * 60)
    print("Summary:")
    print(f"  Domain radius:    {new_r.max()/1e3:.1f} km")
    print(f"  Elements:         {num_new_elem} (was {num_elem})")
    print(f"  Nodes:            {num_new_nodes} (was {num_nodes})")
    print(f"  Absorbing edges:  {len(abs_edges)}")
    print(f"  Free surface:     0 (none)")
    print(f"\n  Output prefix: *{SUFFIX}")
    print("=" * 60)

    # ---- Verify source/receiver locations are inside ----
    print("\nVerifying source/receiver locations...")
    locations = {
        "Source (elem 3)": (-775207.84, 0.0),
        "REC01": (0.0, -775207.84),
        "REC02": (-1312983.0, 0.0),
        "REC03": (-775207.84, 0.0),
    }
    for name, (x, z) in locations.items():
        r = np.sqrt(x**2 + z**2)
        inside = "OK" if r < new_r.max() else "OUTSIDE!"
        print(f"  {name}: r = {r/1e3:.1f} km  [{inside}]")


if __name__ == "__main__":
    main()
