#!/usr/bin/env python3
"""
Compute DFDM element center coordinates for SPECFEM2D station placement.

Reads DFDM grid files (grid_x_N, grid_z_N) and config.toml to determine the
physical coordinates of source, receivers, and virtual receivers. Adapts
automatically when precision parameters (ppw, order) change the grid.

Usage:
  # Interactive: print coordinates
  python3 compute_dfdm_coordinates.py --dfdm-data-dir /path/to/sample_out

  # Write STATIONS file and coordinates JSON for pipeline automation
  python3 compute_dfdm_coordinates.py \
      --dfdm-data-dir /path/to/sample_out \
      --config /path/to/config.toml \
      --stations-out /path/to/STATIONS \
      --coords-json-out /path/to/dfdm_coords.json

COORDINATE METHODS
------------------
1. Bilinear interpolate at (0.5, 0.5):
   DFDM's bilinear_interpolate (src/utils.cpp:280-306) maps parametric
   coordinates to physical space. This is where receivers record data.

2. Nearest grid point to mean:
   For elements NOT in receiver_elements, data is extracted from wavefield
   snapshots at the grid point closest to the element center.
"""

import numpy as np
import math
import json
import argparse
import os
import sys

# Defaults (overridden by CLI args)
DEFAULT_DFDM_DATA_DIR = "/Users/yuantian/Desktop/Berkeley_work/DFDM_2D_1.0_bench/build/sample_out"
DEFAULT_SOURCE_ELEM = 3
DEFAULT_RECEIVER_ELEMS = [1, 2, 3]


def parse_config_toml(config_path):
    """Parse a simple TOML config file to extract source_element and receiver_elements.

    Uses basic string parsing (no toml library needed).
    """
    source_elem = DEFAULT_SOURCE_ELEM
    receiver_elems = list(DEFAULT_RECEIVER_ELEMS)

    if not os.path.exists(config_path):
        print(f"  WARNING: Config file not found: {config_path}, using defaults")
        return source_elem, set(receiver_elems)

    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("source_element"):
                source_elem = int(line.split("=")[1].strip())
            elif line.startswith("receiver_elements"):
                # Parse [1, 2, 3] format
                val = line.split("=")[1].strip()
                val = val.strip("[]")
                receiver_elems = [int(x.strip()) for x in val.split(",")]

    return source_elem, set(receiver_elems)


def bilinear_interpolate(grid, Nx, Nz, xlocal, zlocal):
    """Replicate DFDM's bilinear_interpolate exactly (from utils.cpp).

    Maps xlocal in [0,1] to fractional index (Nx-1)*xlocal, then
    bilinearly interpolates between the 4 surrounding grid points.
    """
    ix = (Nx - 1) * xlocal
    iz = (Nz - 1) * zlocal
    ix_f = max(0, min(int(math.floor(ix)), Nx - 1))
    ix_c = max(0, min(int(math.ceil(ix)), Nx - 1))
    iz_f = max(0, min(int(math.floor(iz)), Nz - 1))
    iz_c = max(0, min(int(math.ceil(iz)), Nz - 1))
    rx = ix - ix_f
    rz = iz - iz_f
    return (grid[ix_f, iz_f] * (1 - rx) * (1 - rz) +
            grid[ix_c, iz_f] * rx * (1 - rz) +
            grid[ix_f, iz_c] * (1 - rx) * rz +
            grid[ix_c, iz_c] * rx * rz)


def nearest_grid_point_to_mean(gx, gz):
    """Find the grid point nearest to (gx.mean(), gz.mean()).

    This is how virtual receivers are extracted from snapshots.
    Returns (grid_point_x, grid_point_z, row_index, col_index).
    """
    ax, az = gx.mean(), gz.mean()
    dist = (gx - ax)**2 + (gz - az)**2
    ci, cj = np.unravel_index(np.argmin(dist), gx.shape)
    return gx[ci, cj], gz[ci, cj], ci, cj


def find_virtual_receivers_on_grid(data_dir, source_elem, source_x, source_z,
                                    target_distances_km=(50, 100, 200, 300),
                                    margin_frac=0.1):
    """Pick grid points within the source element at target distances from source.

    Selects exact DFDM grid points so virtual receiver extraction is precise.
    Prefers z > 0 direction (angular) at roughly constant depth.

    Args:
        data_dir: DFDM output directory containing grid files
        source_elem: source element ID
        source_x, source_z: source coordinates (bilinear center)
        target_distances_km: desired distances from source in km
        margin_frac: fraction of grid extent to avoid at boundaries

    Returns:
        list of (rec_name, x, z, row_idx, col_idx) tuples
    """
    gx = np.genfromtxt(os.path.join(data_dir, f"grid_x_{source_elem}"), delimiter=',')
    gz = np.genfromtxt(os.path.join(data_dir, f"grid_z_{source_elem}"), delimiter=',')
    Nx, Nz = gx.shape

    # Compute distance from source for all grid points
    dist_from_source = np.sqrt((gx - source_x)**2 + (gz - source_z)**2)

    # Create margin mask to avoid boundary grid points
    margin_r = int(max(1, Nx * margin_frac))
    margin_c = int(max(1, Nz * margin_frac))
    interior_mask = np.zeros_like(gx, dtype=bool)
    interior_mask[margin_r:Nx-margin_r, margin_c:Nz-margin_c] = True

    # Prefer z > 0 direction (positive angular direction)
    prefer_mask = interior_mask & (gz > source_z)

    virtual_receivers = []
    rec_idx = 10  # start at REC10

    for target_km in target_distances_km:
        target_m = target_km * 1e3

        # Find grid point closest to target distance, preferring z > 0
        cost = np.abs(dist_from_source - target_m)
        cost[~prefer_mask] = 1e20  # penalize boundary and z < 0 points

        best_idx = np.unravel_index(np.argmin(cost), gx.shape)
        best_x = float(gx[best_idx])
        best_z = float(gz[best_idx])
        actual_dist = float(dist_from_source[best_idx])

        rec_name = f"REC{rec_idx}"
        virtual_receivers.append((rec_name, best_x, best_z, int(best_idx[0]), int(best_idx[1])))
        print(f"    {rec_name}: grid[{best_idx[0]},{best_idx[1]}] = "
              f"({best_x:.1f}, {best_z:.1f}) m, "
              f"dist={actual_dist/1e3:.1f} km (target {target_km} km)")
        rec_idx += 1

    return virtual_receivers


def count_elements(data_dir):
    """Count the number of elements by looking for grid_x_N files."""
    n = 0
    while os.path.exists(os.path.join(data_dir, f"grid_x_{n}")):
        n += 1
    return n


def compute_all_coordinates(data_dir, source_elem, receiver_elems):
    """Compute coordinates for all elements.

    Returns dict: eid -> {x, z, method, grid_shape}
    """
    n_elem = count_elements(data_dir)
    if n_elem == 0:
        print(f"ERROR: No grid files found in {data_dir}")
        sys.exit(1)

    all_coords = {}

    for eid in range(n_elem):
        gx = np.genfromtxt(os.path.join(data_dir, f"grid_x_{eid}"), delimiter=',')
        gz = np.genfromtxt(os.path.join(data_dir, f"grid_z_{eid}"), delimiter=',')
        Nx, Nz = gx.shape

        # bilinear_interpolate at (0.5, 0.5) — where receivers record
        bil_x = bilinear_interpolate(gx, Nx, Nz, 0.5, 0.5)
        bil_z = bilinear_interpolate(gz, Nx, Nz, 0.5, 0.5)

        # nearest grid point to mean — where virtual receivers are
        ngp_x, ngp_z, ci, cj = nearest_grid_point_to_mean(gx, gz)

        if eid in receiver_elems:
            use_x, use_z = bil_x, bil_z
            method = "bilinear"
        else:
            use_x, use_z = ngp_x, ngp_z
            method = "nearest_grid_point"

        all_coords[eid] = {
            'x': float(use_x), 'z': float(use_z),
            'bil_x': float(bil_x), 'bil_z': float(bil_z),
            'ngp_x': float(ngp_x), 'ngp_z': float(ngp_z),
            'method': method,
            'grid_shape': [int(Nx), int(Nz)],
        }

    return all_coords


def get_source_bilinear_coordinates(all_coords, source_elem):
    """Get source coordinates using DFDM's actual method (always bilinear)."""
    src = all_coords[source_elem]
    return float(src['bil_x']), float(src['bil_z'])


def write_stations_file(all_coords, source_elem, stations_path, virtual_receivers=None):
    """Write SPECFEM2D STATIONS file with coordinates matching DFDM grid.

    Includes both element-center receivers (REC01-09) and virtual receivers
    (REC10+) so SPECFEM records at the same locations as DFDM virtual extraction.
    """
    # Element-to-SPECFEM receiver mapping
    # REC01-03 = high-res receivers (elem 1,2,3)
    # REC04-09 = virtual receivers (elem 0,4,5,6,7,8)
    elem_order = [1, 2, 3, 0, 4, 5, 6, 7, 8]
    rec_names = ["REC01", "REC02", "REC03", "REC04", "REC05",
                 "REC06", "REC07", "REC08", "REC09"]

    lines = []
    for rec_name, eid in zip(rec_names, elem_order):
        if eid not in all_coords:
            continue
        c = all_coords[eid]
        lines.append(f"{rec_name:<10s} SY  {c['x']:16.4f}  {c['z']:16.4f}  0.0  0.0")

    # Add virtual receivers (REC10+) at exact grid-point positions
    if virtual_receivers:
        for rec_name, vx, vz, _, _ in virtual_receivers:
            lines.append(f"{rec_name:<10s} SY  {vx:16.4f}  {vz:16.4f}  0.0  0.0")

    with open(stations_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    n_virt = len(virtual_receivers) if virtual_receivers else 0
    print(f"  Wrote STATIONS file: {stations_path} ({len(lines)} receivers: 9 element + {n_virt} virtual)")


def write_coords_json(all_coords, source_elem, json_path, virtual_receivers=None):
    """Write coordinates JSON for compare_snapshots.py to consume."""
    # Build the structure that compare_snapshots.py expects
    elem_order = [1, 2, 3, 0, 4, 5, 6, 7, 8]
    rec_names = ["REC01", "REC02", "REC03", "REC04", "REC05",
                 "REC06", "REC07", "REC08", "REC09"]

    source_x, source_z = get_source_bilinear_coordinates(all_coords, source_elem)

    output = {
        'source_element': source_elem,
        'source_x': source_x,
        'source_z': source_z,
        'source_method': 'bilinear',
        'receivers': {},
        'all_receivers': {},
    }

    # SPECFEM receiver mapping
    for rec_name, eid in zip(rec_names, elem_order):
        if eid not in all_coords:
            continue
        c = all_coords[eid]
        output['receivers'][rec_name] = {
            'elem_id': eid, 'x': c['x'], 'z': c['z'], 'method': c['method'],
        }

    # All elements (E0-E8) for compare_snapshots.py
    for eid, c in sorted(all_coords.items()):
        output['all_receivers'][f"E{eid}"] = {
            'elem_id': eid, 'x': c['x'], 'z': c['z'],
            'method': c['method'], 'grid_shape': c['grid_shape'],
        }

    # Add virtual receivers at exact DFDM grid points within source element
    if virtual_receivers:
        for rec_name, virt_x, virt_z, row_idx, col_idx in virtual_receivers:
            output['all_receivers'][rec_name] = {
                'elem_id': source_elem,
                'x': float(virt_x),
                'z': float(virt_z),
                'method': 'virtual',
                'grid_shape': all_coords[source_elem]['grid_shape'],
                'grid_index': [row_idx, col_idx],
            }

    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2)

    n_virt = len(virtual_receivers) if virtual_receivers else 0
    print(f"  Wrote coords JSON: {json_path} (9 element receivers + {n_virt} virtual)")


def main():
    parser = argparse.ArgumentParser(
        description="Compute DFDM element coordinates and write SPECFEM2D STATIONS file")
    parser.add_argument('--dfdm-data-dir', default=DEFAULT_DFDM_DATA_DIR,
                        help='DFDM output directory containing grid_x_N files')
    parser.add_argument('--config', default=None,
                        help='DFDM config.toml (to read source_element, receiver_elements)')
    parser.add_argument('--stations-out', default=None,
                        help='Path to write SPECFEM2D STATIONS file')
    parser.add_argument('--coords-json-out', default=None,
                        help='Path to write coordinates JSON file')
    args = parser.parse_args()

    # Parse config to get source/receiver element IDs
    if args.config:
        source_elem, receiver_elems = parse_config_toml(args.config)
    else:
        source_elem = DEFAULT_SOURCE_ELEM
        receiver_elems = set(DEFAULT_RECEIVER_ELEMS)

    print("=" * 90)
    print("DFDM Element Center Coordinates")
    print(f"  data_dir:         {args.dfdm_data_dir}")
    print(f"  source_element:   {source_elem}")
    print(f"  receiver_elems:   {sorted(receiver_elems)}")
    print("=" * 90)

    all_coords = compute_all_coordinates(args.dfdm_data_dir, source_elem, receiver_elems)

    # Print summary
    src_x, src_z = get_source_bilinear_coordinates(all_coords, source_elem)
    for eid in sorted(all_coords.keys()):
        c = all_coords[eid]
        dist = np.sqrt((c['x'] - src_x)**2 + (c['z'] - src_z)**2)
        marker = " <-- SOURCE" if eid == source_elem else ""
        print(f"  E{eid}  grid={c['grid_shape'][0]}x{c['grid_shape'][1]}"
              f"  ({c['x']:14.2f}, {c['z']:14.2f}) m"
              f"  dist={dist/1e3:8.1f} km  [{c['method']}]{marker}")

    # Find virtual receiver grid points within source element
    print(f"\n  Finding virtual receivers on element {source_elem} grid...")
    virtual_receivers = find_virtual_receivers_on_grid(
        args.dfdm_data_dir, source_elem, src_x, src_z,
        target_distances_km=(50, 100, 200, 300),
    )

    # Write output files if requested
    if args.stations_out:
        write_stations_file(all_coords, source_elem, args.stations_out, virtual_receivers)

    if args.coords_json_out:
        write_coords_json(all_coords, source_elem, args.coords_json_out, virtual_receivers)

    # Always print STATIONS content for reference
    elem_order = [1, 2, 3, 0, 4, 5, 6, 7, 8]
    rec_names = ["REC01", "REC02", "REC03", "REC04", "REC05",
                 "REC06", "REC07", "REC08", "REC09"]
    print("\nSPECFEM2D STATIONS content:")
    for rec_name, eid in zip(rec_names, elem_order):
        if eid in all_coords:
            c = all_coords[eid]
            print(f"  {rec_name:<10s} SY  {c['x']:16.4f}  {c['z']:16.4f}  0.0  0.0")
    for rec_name, vx, vz, ri, ci in virtual_receivers:
        print(f"  {rec_name:<10s} SY  {vx:16.4f}  {vz:16.4f}  0.0  0.0  [grid {ri},{ci}]")


if __name__ == "__main__":
    main()
