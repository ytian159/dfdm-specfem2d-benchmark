#!/usr/bin/env python3
"""
Compare DFDM and SPECFEM2D pressure snapshots and pressure waveforms.

Generates one pressure figure per time step:
1. Pressure snapshots (DFDM vs SPECFEM2D + difference)
2. Pressure waveform comparison at all receivers
"""

import json
import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from scipy.interpolate import griddata

# ============================================================
# Configuration  (overridable via environment variables)
# ============================================================
DFDM_DATA_DIR = os.environ.get(
    "BENCH_DFDM_DATA_DIR",
    "/Users/yuantian/Desktop/Berkeley_work/DFDM_2D_1.0_bench/build/sample_out",
)
SPECFEM_OUTPUT_DIR = os.environ.get(
    "BENCH_SPECFEM_OUTPUT_DIR",
    "/Users/yuantian/Desktop/Berkeley_work/specfem2d_ak135f/OUTPUT_FILES",
)
FIGURE_OUTPUT_DIR = os.environ.get("BENCH_FIGURE_DIR", SPECFEM_OUTPUT_DIR)
DFDM_CONFIG = os.environ.get(
    "BENCH_DFDM_CONFIG",
    "/Users/yuantian/Desktop/Berkeley_work/DFDM_2D_1.0_bench/config/config.toml",
)

DT = float(os.environ.get("BENCH_DT", "0.1"))
TOTAL_NSTEP = int(os.environ.get("BENCH_NSTEP", "16000"))
DFDM_TOTAL_ELEM = 9
DFDM_SNAP_INTERVAL = int(os.environ.get("BENCH_DFDM_SNAP_INTERVAL", "0"))


def _auto_detect_snap_interval(data_dir, total_elem=9):
    """Auto-detect DFDM snapshot interval from available files."""
    import glob
    files = sorted(glob.glob(f"{data_dir}/elem_0_*.out"))
    if len(files) < 3:
        return 50  # fallback
    steps = []
    for f in files:
        base = os.path.basename(f)
        step = int(base.replace("elem_0_", "").replace(".out", ""))
        if step > 0:
            steps.append(step)
    steps.sort()
    if len(steps) >= 2:
        diffs = [steps[i+1] - steps[i] for i in range(min(10, len(steps)-1))]
        from collections import Counter
        most_common = Counter(diffs).most_common(1)[0][0]
        return most_common
    return 50


if DFDM_SNAP_INTERVAL == 0:
    DFDM_SNAP_INTERVAL = _auto_detect_snap_interval(DFDM_DATA_DIR)
    print(f"  Auto-detected DFDM snapshot interval: {DFDM_SNAP_INTERVAL}")

DFDM_SOURCE_PEAK = 100.0
SPECFEM_SOURCE_PEAK = 0.0
SPECFEM_T0 = 180.0
SPECFEM_DUMP_INTERVAL = int(os.environ.get("BENCH_SPECFEM_DUMP_INTERVAL", "1000"))
TRACE_TAIL_TRIM = float(os.environ.get("BENCH_TRACE_TAIL_TRIM", "100.0"))

VP = 2750.0
RHO = 2000.0

COORDS_JSON = os.environ.get("BENCH_COORDS_JSON", "")



def _load_coordinates():
    if not COORDS_JSON:
        raise FileNotFoundError(
            "BENCH_COORDS_JSON environment variable is not set. "
            "Please set it to the path of the coordinates JSON file."
        )
    if not os.path.exists(COORDS_JSON):
        raise FileNotFoundError(
            f"Coordinates JSON file not found: {COORDS_JSON}"
        )
    with open(COORDS_JSON) as f:
        data = json.load(f)
    source = (data["source_x"], data["source_z"])
    receivers = {}
    for name, info in data["all_receivers"].items():
        receivers[name] = (info["elem_id"], info["x"], info["z"])
    print(f"  Loaded coordinates from: {COORDS_JSON}")
    return source, receivers


SOURCE, ALL_RECEIVERS = _load_coordinates()

DFDM_ELEM_TO_SPECFEM = {
    1: "REC01",
    2: "REC02",
    3: "REC03",
    0: "REC04",
    4: "REC05",
    5: "REC06",
    6: "REC07",
    7: "REC08",
    8: "REC09",
}

# Virtual receivers map to source element for SPECFEM lookup
VIRTUAL_RECEIVER_SPECFEM_NAMES = {
    "REC10": "REC10",
    "REC11": "REC11",
    "REC12": "REC12",
    "REC13": "REC13",
}


def _load_dfdm_receiver_elements(config_path):
    receiver_elems = {1, 2, 3}
    if not os.path.exists(config_path):
        print(f"  WARNING: DFDM config not found: {config_path}, using {sorted(receiver_elems)}")
        return receiver_elems

    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("receiver_elements"):
                continue
            value = line.split("=", 1)[1].strip().strip("[]")
            if not value:
                receiver_elems = set()
            else:
                receiver_elems = {int(x.strip()) for x in value.split(",")}
            break

    print(f"  Receiver elements from config: {sorted(receiver_elems)}")
    return receiver_elems


DFDM_RECEIVER_ELEMS = _load_dfdm_receiver_elements(DFDM_CONFIG)

ELEM_COLORS = {
    0: "tab:brown",
    1: "tab:green",
    2: "tab:orange",
    3: "tab:purple",
    4: "tab:cyan",
    5: "tab:pink",
    6: "tab:red",
    7: "tab:blue",
    8: "tab:olive",
}
ELEM_MARKERS = {
    0: "v",
    1: "s",
    2: "^",
    3: "o",
    4: "D",
    5: "p",
    6: "h",
    7: "*",
    8: "X",
}

R_DOMAIN = 1700.0


def build_comparison_pairs():
    pairs = []
    for sem_step in range(SPECFEM_DUMP_INTERVAL, TOTAL_NSTEP + 1, SPECFEM_DUMP_INTERVAL):
        sem_time = (sem_step - 1) * DT - SPECFEM_T0
        time_after_peak = sem_time - SPECFEM_SOURCE_PEAK

        if time_after_peak < 50:
            continue

        dfdm_time = time_after_peak + DFDM_SOURCE_PEAK
        dfdm_step_exact = dfdm_time / DT
        dfdm_step = int(round(dfdm_step_exact / DFDM_SNAP_INTERVAL)) * DFDM_SNAP_INTERVAL

        if dfdm_step < 0 or dfdm_step > TOTAL_NSTEP:
            continue

        pairs.append((time_after_peak, dfdm_step, sem_step))
    return pairs


def load_dfdm_pressure_snapshot(data_dir, step, total_elem, snap_interval, dt, rho):
    step_prev = step - snap_interval
    step_next = step + snap_interval
    snap_dt = snap_interval * dt

    all_x, all_z, all_p = [], [], []
    for elem_id in range(total_elem):
        f_prev = f"{data_dir}/elem_{elem_id}_{step_prev}.out"
        f_curr = f"{data_dir}/elem_{elem_id}_{step}.out"
        f_next = f"{data_dir}/elem_{elem_id}_{step_next}.out"
        if not (os.path.exists(f_prev) and os.path.exists(f_curr) and os.path.exists(f_next)):
            continue
        gx = np.genfromtxt(f"{data_dir}/grid_x_{elem_id}", delimiter=",")
        gz = np.genfromtxt(f"{data_dir}/grid_z_{elem_id}", delimiter=",")
        v_prev = np.genfromtxt(f_prev, delimiter=",")
        v_curr = np.genfromtxt(f_curr, delimiter=",")
        v_next = np.genfromtxt(f_next, delimiter=",")
        d2u = (v_next - 2.0 * v_curr + v_prev) / (snap_dt ** 2)
        pressure = rho * d2u
        all_x.append(gx.flatten())
        all_z.append(gz.flatten())
        all_p.append(pressure.flatten())

    if not all_x:
        return None, None, None
    return np.concatenate(all_x), np.concatenate(all_z), np.concatenate(all_p)


def load_specfem_pressure_snapshot(output_dir, step):
    filename = f"{output_dir}/wavefield{step:07d}_01.txt"
    if not os.path.exists(filename):
        return None
    return np.loadtxt(filename).flatten()


def load_specfem_pressure_seismogram(output_dir, receiver):
    filename = f"{output_dir}/SY.{receiver}.PRE.semp"
    if not os.path.exists(filename):
        return None, None
    data = np.loadtxt(filename)
    return data[:, 0], data[:, 1]


def load_dfdm_receiver_pressure(data_dir, elem_id, rho):
    """Load DFDM receiver data and return (time, pressure) arrays.

    Reads the recorded values file, computes the scalar field as the
    average of the two recorded components, then derives pressure as
    rho * d^2u/dt^2 using second-order finite differences.
    """
    ur_file = f"{data_dir}/receiver_elem_{elem_id}_recorded_values.out"
    if not os.path.exists(ur_file):
        raise FileNotFoundError(
            f"Missing DFDM receiver file: {ur_file}. "
            "Ensure receiver_elements in config.toml includes this element and rerun DFDM."
        )
    ur = np.genfromtxt(ur_file, delimiter=",")
    t = ur[:, 0]
    u = (ur[:, 1] + ur[:, 2]) / 2.0

    dt = t[1] - t[0] if len(t) > 1 else 1.0
    n = len(u)
    d2u_dt2 = np.zeros(n)
    d2u_dt2[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / (dt ** 2)
    if n > 2:
        d2u_dt2[0] = (u[2] - 2 * u[1] + u[0]) / (dt ** 2)
        d2u_dt2[-1] = (u[-1] - 2 * u[-2] + u[-3]) / (dt ** 2)

    return t, dt, rho * d2u_dt2


def load_dfdm_virtual_receiver_pressure(data_dir, elem_id, virt_x, virt_z, rho, dt, snap_interval):
    """Extract DFDM pressure time series from snapshots at nearest grid point.

    For virtual receivers that don't have dedicated receiver files,
    extract data from wavefield snapshots at the nearest grid point.
    The actual snapshot interval is auto-detected from the element's files,
    since the source element may have every-timestep snapshots while other
    elements are saved less frequently.
    """
    # Load grid for the element
    gx = np.genfromtxt(f"{data_dir}/grid_x_{elem_id}", delimiter=",")
    gz = np.genfromtxt(f"{data_dir}/grid_z_{elem_id}", delimiter=",")

    # Find nearest grid point to virtual receiver
    dist = np.sqrt((gx - virt_x)**2 + (gz - virt_z)**2)
    min_idx = np.unravel_index(np.argmin(dist), gx.shape)

    # Read snapshots at this grid point
    import glob
    snapshot_files = sorted(glob.glob(f"{data_dir}/elem_{elem_id}_*.out"))
    if len(snapshot_files) < 3:
        print(f"  WARNING: Not enough snapshots for elem {elem_id}, falling back to receiver file")
        return load_dfdm_receiver_pressure(data_dir, elem_id, rho)

    # Extract step numbers from filenames
    steps = []
    for f in snapshot_files:
        base = os.path.basename(f)
        step = int(base.replace(f"elem_{elem_id}_", "").replace(".out", ""))
        steps.append(step)
    steps = np.array(sorted(steps))

    # Auto-detect actual snapshot interval for this element
    if len(steps) >= 2:
        step_diffs = np.diff(steps)
        from collections import Counter
        actual_interval = Counter(step_diffs.tolist()).most_common(1)[0][0]
    else:
        actual_interval = snap_interval  # fallback to global
    actual_snap_dt = actual_interval * dt

    # Load potential values at the grid point from each snapshot
    u_timeseries = []
    for step in steps:
        snap_file = f"{data_dir}/elem_{elem_id}_{step}.out"
        if os.path.exists(snap_file):
            u_snap = np.genfromtxt(snap_file, delimiter=",")
            u_timeseries.append(u_snap[min_idx])
        else:
            u_timeseries.append(0.0)

    u_timeseries = np.array(u_timeseries)
    t = steps * dt

    # Compute pressure from second time derivative using actual interval
    n = len(u_timeseries)
    d2u_dt2 = np.zeros(n)
    d2u_dt2[1:-1] = (u_timeseries[2:] - 2 * u_timeseries[1:-1] + u_timeseries[:-2]) / (actual_snap_dt ** 2)
    if n > 2:
        d2u_dt2[0] = (u_timeseries[2] - 2 * u_timeseries[1] + u_timeseries[0]) / (actual_snap_dt ** 2)
        d2u_dt2[-1] = (u_timeseries[-1] - 2 * u_timeseries[-2] + u_timeseries[-3]) / (actual_snap_dt ** 2)

    return t, actual_snap_dt, rho * d2u_dt2


def plot_snapshot(ax, x_km, z_km, values, title, vmin, vmax):
    mask_r = np.sqrt(x_km ** 2 + z_km ** 2) <= R_DOMAIN * 1.05
    x_plot, z_plot, v_plot = x_km[mask_r], z_km[mask_r], values[mask_r]

    triang = tri.Triangulation(x_plot, z_plot)
    edge_max = 50
    xi, zi = x_plot[triang.triangles], z_plot[triang.triangles]
    mask = (np.ptp(xi, axis=1) > edge_max) | (np.ptp(zi, axis=1) > edge_max)
    triang.set_mask(mask)

    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    im = ax.tripcolor(triang, v_plot, cmap="seismic", norm=norm, shading="gouraud")

    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(R_DOMAIN * np.cos(theta), R_DOMAIN * np.sin(theta), "k-", linewidth=0.8)

    ax.plot(SOURCE[0] / 1e3, SOURCE[1] / 1e3, "r*", markersize=14,
            markeredgecolor="k", markeredgewidth=0.5, zorder=10)

    for name, (eid, rx, rz) in ALL_RECEIVERS.items():
        color = ELEM_COLORS[eid]
        marker = ELEM_MARKERS[eid]
        ax.plot(rx / 1e3, rz / 1e3, marker, color=color,
                markersize=8, markeredgecolor="k", markeredgewidth=0.5, zorder=10)
        dx_off = 30 if rx >= 0 else -30
        dz_off = 30 if rz >= 0 else -30
        ha = "left" if rx >= 0 else "right"
        va = "bottom" if rz >= 0 else "top"
        ax.text(rx / 1e3 + dx_off, rz / 1e3 + dz_off, name,
                fontsize=7, fontweight="bold", color=color, ha=ha, va=va, zorder=11)

    ax.set_xlim(-R_DOMAIN * 1.05, R_DOMAIN * 1.05)
    ax.set_ylim(-R_DOMAIN * 1.05, R_DOMAIN * 1.05)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("X (km)", fontsize=9)
    ax.set_ylabel("Z (km)", fontsize=9)
    ax.tick_params(labelsize=8)
    return im


def _add_time_markers(ax, time_after_peak, expected):
    if time_after_peak is not None:
        ax.axvline(x=time_after_peak, color="gray", linestyle=":", linewidth=2, alpha=0.5)
    ax.axvline(x=expected, color="red", linestyle="--", linewidth=1.5, alpha=0.7, zorder=5)
    ax.text(expected + 5, 0.92, f"t_arr={expected:.0f}s",
            fontsize=7, color="red", fontweight="bold", alpha=0.8,
            transform=ax.get_xaxis_transform())


def _common_trace_window(t_d, t_s, tail_trim):
    t_start = max(np.min(t_d), np.min(t_s))
    t_end = min(np.max(t_d), np.max(t_s)) - tail_trim
    if t_end <= t_start:
        t_end = min(np.max(t_d), np.max(t_s))
    return t_start, t_end


def compute_pressure_difference(dfdm_data, specfem_data):
    t_d = dfdm_data["t"]
    y_d = dfdm_data["pressure"]
    if specfem_data is None:
        return t_d, y_d, np.zeros_like(y_d), np.zeros_like(y_d), 0.0, np.min(t_d), np.max(t_d)

    t_s = specfem_data["t_pre"]
    y_s = specfem_data["pressure"]
    if t_s is None or y_s is None:
        return t_d, y_d, np.zeros_like(y_d), np.zeros_like(y_d), 0.0, np.min(t_d), np.max(t_d)

    t_start, t_end = _common_trace_window(t_d, t_s, TRACE_TAIL_TRIM)
    mask = (t_d >= t_start) & (t_d <= t_end)
    if not np.any(mask):
        t_start = max(np.min(t_d), np.min(t_s))
        t_end = min(np.max(t_d), np.max(t_s))
        mask = (t_d >= t_start) & (t_d <= t_end)

    t_common = t_d[mask]
    y_d_common = y_d[mask]
    y_s_interp = np.interp(t_common, t_s, y_s)
    diff = y_d_common - y_s_interp

    max_ref = np.max(np.abs(y_s_interp))
    rel_err = (np.max(np.abs(diff)) / max_ref * 100.0) if max_ref > 0 else 0.0
    return t_common, y_d_common, y_s_interp, diff, rel_err, t_start, t_end


def plot_pressure_panel(ax, rec_name, elem_id, time_after_peak, dfdm_data, specfem_data):
    color = ELEM_COLORS[elem_id]
    _, rx, rz = ALL_RECEIVERS[rec_name]
    dist = np.sqrt((rx - SOURCE[0]) ** 2 + (rz - SOURCE[1]) ** 2) / 1e3
    expected = dist / (VP / 1e3)

    t_common, y_d, y_s, _, _, t_start, t_end = compute_pressure_difference(dfdm_data, specfem_data)
    ax.plot(t_common, y_d, "-", color=color, linewidth=1,
            label=f"DFDM PRE (max={np.max(np.abs(y_d)):.2e})")
    ax.plot(t_common, y_s, "--", color="black", linewidth=1, alpha=0.8,
            label=f"SPECFEM PRE (max={np.max(np.abs(y_s)):.2e})")

    _add_time_markers(ax, time_after_peak, expected)
    ax.set_title(f"{rec_name} PRE  d={dist:.0f}km  t_arr={expected:.0f}s", fontsize=9, color=color)
    ax.legend(fontsize=6, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(t_start, t_end)


def plot_pressure_difference_panel(ax, rec_name, elem_id, time_after_peak, dfdm_data, specfem_data):
    color = ELEM_COLORS[elem_id]
    _, rx, rz = ALL_RECEIVERS[rec_name]
    dist = np.sqrt((rx - SOURCE[0]) ** 2 + (rz - SOURCE[1]) ** 2) / 1e3
    expected = dist / (VP / 1e3)

    t_common, _, _, diff, rel_err, t_start, t_end = compute_pressure_difference(dfdm_data, specfem_data)
    ax.plot(t_common, diff, "-", color="red", linewidth=1.0, alpha=0.9,
            label=f"Diff (max={np.max(np.abs(diff)):.2e})")
    ax.axhline(y=0, color="gray", linewidth=0.5, alpha=0.5)

    _add_time_markers(ax, time_after_peak, expected)
    ax.set_title(f"{rec_name} PRE diff  d={dist:.0f}km  rel={rel_err:.1f}%", fontsize=9, color=color)
    ax.legend(fontsize=6, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(t_start, t_end)


def main():
    print("=" * 60)
    print("DFDM vs SPECFEM2D Pressure Comparison")
    print("Output: one pressure figure per time step")
    print(f"  DT={DT}, NSTEP={TOTAL_NSTEP}, SNAP_INTERVAL={DFDM_SNAP_INTERVAL}")
    print(f"  DFDM config: {DFDM_CONFIG}")
    print(f"  Trace tail trim: {TRACE_TAIL_TRIM:.1f} s")
    print(f"  DFDM dir:    {DFDM_DATA_DIR}")
    print(f"  SPECFEM dir: {SPECFEM_OUTPUT_DIR}")
    print(f"  Figure dir:  {FIGURE_OUTPUT_DIR}")
    print("=" * 60)

    pairs = build_comparison_pairs()
    print("\nComparison pairs:")
    for tap, d_step, s_step in pairs:
        print(f"  after_peak={tap:7.1f}s  DFDM={d_step:5d}  SPECFEM={s_step:5d}")

    print("\nLoading grids...")
    specfem_grid = np.loadtxt(f"{SPECFEM_OUTPUT_DIR}/wavefield_grid_for_dumps.txt")
    sx_km, sz_km = specfem_grid[:, 0] / 1e3, specfem_grid[:, 1] / 1e3

    needed_receiver_elems = {elem_id for elem_id, _, _ in ALL_RECEIVERS.values()}
    missing_from_config = needed_receiver_elems - DFDM_RECEIVER_ELEMS
    if missing_from_config:
        raise RuntimeError(
            "DFDM config receiver_elements is missing required receiver elements: "
            f"{sorted(missing_from_config)}. "
            "Update config.toml and rerun DFDM so all receiver files exist."
        )

    print("\nLoading DFDM receiver data...")
    dfdm_rec = {}
    for rec_name, (elem_id, rx, rz) in ALL_RECEIVERS.items():
        # Check if this is a virtual receiver
        is_virtual = rec_name in VIRTUAL_RECEIVER_SPECFEM_NAMES

        if is_virtual:
            # Extract from snapshots at the virtual receiver location
            t, rec_dt, pressure = load_dfdm_virtual_receiver_pressure(
                DFDM_DATA_DIR, elem_id, rx, rz, RHO, DT, DFDM_SNAP_INTERVAL
            )
        else:
            # Use dedicated receiver file
            t, rec_dt, pressure = load_dfdm_receiver_pressure(DFDM_DATA_DIR, elem_id, RHO)

        t_aligned = t - DFDM_SOURCE_PEAK

        dfdm_rec[rec_name] = {"t": t_aligned, "pressure": pressure}
        dist = np.sqrt((rx - SOURCE[0]) ** 2 + (rz - SOURCE[1]) ** 2) / 1e3
        method_str = "[virtual from snapshots]" if is_virtual else ""
        print(
            f"  {rec_name} (elem {elem_id}, dt={rec_dt:.2f}s): "
            f"dist={dist:.0f}km, pressure max={np.max(np.abs(pressure)):.3e} {method_str}"
        )

    print("\nLoading SPECFEM2D pressure seismograms...")
    specfem_rec = {}
    for rec_name, (elem_id, _, _) in ALL_RECEIVERS.items():
        # Check if this is a virtual receiver or element receiver
        if rec_name in VIRTUAL_RECEIVER_SPECFEM_NAMES:
            sem_name = VIRTUAL_RECEIVER_SPECFEM_NAMES[rec_name]
        else:
            sem_name = DFDM_ELEM_TO_SPECFEM.get(elem_id, None)
            if sem_name is None:
                print(f"  {rec_name}: No SPECFEM mapping for element {elem_id}")
                continue

        t_pre, v_pre = load_specfem_pressure_seismogram(SPECFEM_OUTPUT_DIR, sem_name)
        if t_pre is None:
            specfem_rec[rec_name] = None
            print(f"  {rec_name} ({sem_name}): PRE not found")
            continue

        specfem_rec[rec_name] = {"t_pre": t_pre, "pressure": v_pre}
        print(f"  {rec_name} ({sem_name}): pressure max={np.max(np.abs(v_pre)):.3e}")

    rec_sorted = sorted(
        ALL_RECEIVERS.keys(),
        key=lambda r: np.sqrt((ALL_RECEIVERS[r][1] - SOURCE[0]) ** 2 + (ALL_RECEIVERS[r][2] - SOURCE[1]) ** 2),
    )

    for time_after_peak, dfdm_step, specfem_step in pairs:
        print(f"\nTime after source peak = {time_after_peak:.1f}s ...")
        tag = int(time_after_peak)

        dfdm_x, dfdm_z, dfdm_p = load_dfdm_pressure_snapshot(
            DFDM_DATA_DIR, dfdm_step, DFDM_TOTAL_ELEM, DFDM_SNAP_INTERVAL, DT, RHO
        )
        if dfdm_x is None:
            print(
                f"  DFDM pressure snapshot at step {dfdm_step} not found "
                f"(need steps {dfdm_step-DFDM_SNAP_INTERVAL} and {dfdm_step+DFDM_SNAP_INTERVAL}), skipping"
            )
            continue

        specfem_v = load_specfem_pressure_snapshot(SPECFEM_OUTPUT_DIR, specfem_step)
        if specfem_v is None:
            print(f"  SPECFEM pressure snapshot at step {specfem_step} not found, skipping")
            continue

        dx_km, dz_km = dfdm_x / 1e3, dfdm_z / 1e3
        dfdm_t_phys = dfdm_step * DT
        sem_t_phys = (specfem_step - 1) * DT - SPECFEM_T0

        dfdm_absmax = np.percentile(np.abs(dfdm_p[dfdm_p != 0]), 99) if np.any(dfdm_p != 0) else 1
        specfem_absmax = np.percentile(np.abs(specfem_v[specfem_v != 0]), 99) if np.any(specfem_v != 0) else 1

        n_rec = len(rec_sorted)
        NCOLS = 3
        n_wav_rows = (n_rec + NCOLS - 1) // NCOLS
        n_rows = 2 + 2 * n_wav_rows
        height_ratios = [1.3, 0.05] + [0.75, 0.35] * n_wav_rows

        fig = plt.figure(figsize=(26, 3 + 4.0 * n_wav_rows + 5))
        gs = gridspec.GridSpec(n_rows, NCOLS, figure=fig,
                               height_ratios=height_ratios,
                               hspace=0.45, wspace=0.25)

        ax_dfdm = fig.add_subplot(gs[0, 0])
        ax_spec = fig.add_subplot(gs[0, 1])
        ax_diff_snap = fig.add_subplot(gs[0, 2])

        im1 = plot_snapshot(
            ax_dfdm,
            dx_km,
            dz_km,
            dfdm_p,
            f"DFDM pressure  t={time_after_peak:.0f}s \n(raw t_dfdm={dfdm_t_phys:.0f}s, max={np.max(np.abs(dfdm_p)):.2e})",
            -dfdm_absmax,
            dfdm_absmax,
        )
        im2 = plot_snapshot(
            ax_spec,
            sx_km,
            sz_km,
            specfem_v,
            f"SPECFEM2D pressure  t={time_after_peak:.0f}s \n(raw t_seis={sem_t_phys:.0f}s, max={np.max(np.abs(specfem_v)):.2e})",
            -specfem_absmax,
            specfem_absmax,
        )

        dfdm_pts = np.column_stack((dx_km, dz_km))
        specfem_pts = np.column_stack((sx_km, sz_km))
        dfdm_on_specfem = griddata(dfdm_pts, dfdm_p, specfem_pts, method="linear", fill_value=0.0)
        diff_snap = dfdm_on_specfem - specfem_v
        diff_absmax = np.percentile(np.abs(diff_snap[diff_snap != 0]), 99) if np.any(diff_snap != 0) else 1

        im3 = plot_snapshot(
            ax_diff_snap,
            sx_km,
            sz_km,
            diff_snap,
            f"Difference (DFDM - SPECFEM)  t={time_after_peak:.0f}s\nmax|diff|={np.max(np.abs(diff_snap)):.2e}",
            -diff_absmax,
            diff_absmax,
        )

        cax1 = fig.add_subplot(gs[1, 0])
        cax2 = fig.add_subplot(gs[1, 1])
        cax3 = fig.add_subplot(gs[1, 2])
        fig.colorbar(im1, cax=cax1, orientation="horizontal", label="DFDM pressure (+rho*d2U/dt2)")
        fig.colorbar(im2, cax=cax2, orientation="horizontal", label="SPECFEM2D pressure")
        fig.colorbar(im3, cax=cax3, orientation="horizontal", label="Difference (DFDM - SPECFEM)")

        for idx, rec_name in enumerate(rec_sorted):
            row_overlay = 2 + (idx // NCOLS) * 2
            row_diff = row_overlay + 1
            col = idx % NCOLS
            elem_id = ALL_RECEIVERS[rec_name][0]

            ax_ovl = fig.add_subplot(gs[row_overlay, col])
            plot_pressure_panel(ax_ovl, rec_name, elem_id, time_after_peak,
                                dfdm_rec[rec_name], specfem_rec.get(rec_name))
            ax_ovl.set_ylabel("PRE", fontsize=8)

            ax_dif = fig.add_subplot(gs[row_diff, col], sharex=ax_ovl)
            plot_pressure_difference_panel(ax_dif, rec_name, elem_id, time_after_peak,
                                           dfdm_rec[rec_name], specfem_rec.get(rec_name))
            ax_dif.set_ylabel("Diff PRE", fontsize=8)
            if row_diff >= n_rows - 2:
                ax_dif.set_xlabel("Time (s)", fontsize=9)

        fig.suptitle(
            f"DFDM vs SPECFEM2D PRE  |  {time_after_peak:.0f}s   |  Vp={VP:.0f} m/s, rho={RHO:.0f} kg/m^3\n"
            "Top: Pressure snapshots  |  Bottom: PRE waveforms "
            "(Solid=DFDM, Dashed=SPECFEM, Red=Diff)  |  Receivers sorted by distance",
            fontsize=12,
            fontweight="bold",
        )

        outfile = f"{FIGURE_OUTPUT_DIR}/combined_t{tag:04d}.png"
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved pressure figure: {outfile}")

    print("\nDone!")


if __name__ == "__main__":
    main()
