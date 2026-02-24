#!/bin/bash
# ==============================================================================
# DFDM vs SPECFEM2D Benchmarking Script
#
# Runs both solvers on the same 1700 km uniform acoustic domain (Vp=2750 m/s,
# rho=2000 kg/m^3, f0=0.015 Hz) and compares their wavefield snapshots and
# receiver waveforms.
#
# Usage:
#   ./run_benchmark.sh                    # Run full benchmark (baseline profile)
#   ./run_benchmark.sh dfdm               # Run DFDM only
#   ./run_benchmark.sh specfem [NPROC]    # Run SPECFEM2D only (default NPROC=4)
#   ./run_benchmark.sh compare            # Run comparison only
#   ./run_benchmark.sh report             # Print performance summary
#   ./run_benchmark.sh specfem_mesh       # Build SPECFEM2D external mesh
#   ./run_benchmark.sh clean              # Remove output files
#
# Precision profiles (set via PROFILE environment variable):
#   PROFILE=baseline  ./run_benchmark.sh all      # Default settings
#   PROFILE=highres   ./run_benchmark.sh all      # Higher spatial accuracy
#   PROFILE=highres_dt ./run_benchmark.sh all     # Above + halved time step
#
# Individual parameter overrides (environment variables):
#   DFDM_PPW=5 DFDM_ORDER=7 ./run_benchmark.sh dfdm
#   SPECFEM_TIME_SCHEME=2 ./run_benchmark.sh specfem 1
#   BENCH_DT=0.05 ./run_benchmark.sh all
#
# Prerequisites:
#   - DFDM built via setup.sh (or pre-built binaries present)
#   - SPECFEM2D compiled with MPI (xmeshfem2D and xspecfem2D in bin/)
#   - Python 3 with numpy, matplotlib, scipy
#   - OpenMPI (mpirun)
# ==============================================================================

set -e

# ===================== Base Paths (relative to this script) =====================
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# DFDM paths
DFDM_PROJECT_DIR="${SCRIPT_DIR}/DFDM_2D_1.0"
DFDM_BUILD_DIR="${DFDM_PROJECT_DIR}/build"
DFDM_EXEC="${DFDM_BUILD_DIR}/dfdm"
DFDM_CONFIG="${DFDM_PROJECT_DIR}/config/config.toml"
DFDM_MESH_EXEC="${DFDM_PROJECT_DIR}/mesh_gen_ak135/build/mesh_gen"
DFDM_MESH_OUTPUT="${DFDM_PROJECT_DIR}/mesh_gen_ak135/build/output_test"

# SPECFEM2D paths
SPECFEM_DIR="${SCRIPT_DIR}/specfem2d"
SPECFEM_BIN="${SPECFEM_DIR}/bin"
SPECFEM_PAR_FILE="${SPECFEM_DIR}/DATA/Par_file"

# Comparison and coordinate scripts
COMPARE_SCRIPT="${SCRIPT_DIR}/scripts/compare_snapshots.py"
COORD_SCRIPT="${SCRIPT_DIR}/scripts/compute_dfdm_coordinates.py"
SPECFEM_STATIONS="${SPECFEM_DIR}/DATA/STATIONS"

# DFDM MPI tasks
DFDM_NPROC=${DFDM_NPROC:-4}
# Always record all 9 benchmark receivers in DFDM output files
DFDM_RECEIVER_ELEMS=${DFDM_RECEIVER_ELEMS:-"[0, 1, 2, 3, 4, 5, 6, 7, 8]"}

# SPECFEM2D MPI tasks (can be overridden via command line $2)
SPECFEM_NPROC=${2:-4}

# ===================== Precision Parameters =====================
# These control numerical accuracy without recompilation.
# Override via environment variables or use PROFILE presets.
#
# DFDM spatial parameters (changes trigger automatic mesh regeneration):
#   DFDM_PPW           Points per wavelength          (baseline: 3.0)
#   DFDM_ORDER         Polynomial order (b1 and b2)   (baseline: 5)
#   DFDM_GAUSS_ORDER   Gauss quadrature order         (baseline: 6)
#
# SPECFEM2D parameters:
#   SPECFEM_TIME_SCHEME  Time stepping scheme          (baseline: 1)
#       1 = Newmark (2nd order)
#       2 = LDDRK4-6 (4th-order 6-stage Runge-Kutta)
#       3 = classical RK4 (4th-order 4-stage Runge-Kutta)
#
# Shared parameters:
#   BENCH_DT           Time step for both codes (s)   (baseline: 0.1)
#
# Profiles set sensible defaults; env vars override profile values.
# ================================================================

PROFILE=${PROFILE:-baseline}

# Set profile defaults
case "$PROFILE" in
    "baseline")
        _PPW=3.0; _ORDER=5; _GAUSS=6; _SCHEME=1; _DT=0.1
        _SPECFEM_REF_DIR="${SCRIPT_DIR}/specfem2d_mf8"
        ;;
    "highres")
        _PPW=5.0; _ORDER=7; _GAUSS=6; _SCHEME=1; _DT=0.1
        _SPECFEM_REF_DIR="${SCRIPT_DIR}/specfem2d_mf8"
        ;;
    "highres_dt")
        _PPW=5.0; _ORDER=7; _GAUSS=6; _SCHEME=1; _DT=0.05
        ;;
    "reference")
        # Use DFDM baseline params but SPECFEM2D with finer mesh (multiplication_factor=8)
        _PPW=3.0; _ORDER=5; _GAUSS=6; _SCHEME=1; _DT=0.1
        _SPECFEM_REF_DIR="${SCRIPT_DIR}/specfem2d_mf8"
        ;;
    *)
        echo "ERROR: Unknown profile '$PROFILE'. Use: baseline, highres, highres_dt, reference"
        exit 1
        ;;
esac

# Environment variable overrides (take precedence over profile)
DFDM_PPW=${DFDM_PPW:-$_PPW}
DFDM_ORDER=${DFDM_ORDER:-$_ORDER}
DFDM_GAUSS_ORDER=${DFDM_GAUSS_ORDER:-$_GAUSS}
SPECFEM_TIME_SCHEME=${SPECFEM_TIME_SCHEME:-$_SCHEME}
BENCH_DT=${BENCH_DT:-$_DT}

# Allow profile to override SPECFEM2D directory (e.g., reference profile uses mf8 mesh)
if [ -n "${_SPECFEM_REF_DIR:-}" ]; then
    SPECFEM_DIR="$_SPECFEM_REF_DIR"
    SPECFEM_BIN="${SPECFEM_DIR}/bin"
    SPECFEM_PAR_FILE="${SPECFEM_DIR}/DATA/Par_file"
    SPECFEM_STATIONS="${SPECFEM_DIR}/DATA/STATIONS"
fi

# Compute derived parameters
TOTAL_SIM_TIME=1800.0
DT=$BENCH_DT
NSTEP=$(python3 -c "print(int(${TOTAL_SIM_TIME} / ${DT}))")
TOTAL_TIME=$TOTAL_SIM_TIME

# SPECFEM2D snapshot/seismogram intervals (maintain ~100s physical interval)
SPECFEM_SNAP_INTERVAL=$(python3 -c "print(int(100.0 / ${DT}))")

# ===================== Output Directories (profile-tagged) =====================
if [ "$PROFILE" = "baseline" ] || [ "$PROFILE" = "reference" ]; then
    # baseline: default dirs; reference: same DFDM params as baseline, and
    # mf8 SPECFEM_DIR is already a separate folder so use its default OUTPUT_FILES
    DFDM_OUTPUT_DIR="${DFDM_BUILD_DIR}/sample_out"
    SPECFEM_OUTPUT_DIR="${SPECFEM_DIR}/OUTPUT_FILES"
    FIGURE_OUTPUT_DIR="${SPECFEM_OUTPUT_DIR}"
else
    DFDM_OUTPUT_DIR="${DFDM_BUILD_DIR}/sample_out_${PROFILE}"
    SPECFEM_OUTPUT_DIR="${SPECFEM_DIR}/OUTPUT_FILES_${PROFILE}"
    FIGURE_OUTPUT_DIR="${SPECFEM_OUTPUT_DIR}"
fi
DFDM_LOG_DIR="${DFDM_OUTPUT_DIR}/domain/logs"

# Coordinates JSON (profile-tagged, written by compute_dfdm_coordinates.py)
COORDS_JSON="${DFDM_OUTPUT_DIR}/dfdm_coords.json"

# Performance log (profile-tagged)
PERF_LOG="${SCRIPT_DIR}/benchmark_performance_${PROFILE}.log"

# ===================== Utilities =====================
function log_step() {
    echo ""
    echo "=================================================================="
    echo "  $1"
    echo "=================================================================="
}

function log_time() {
    local label="$1"
    local t_start="$2"
    local t_end="$3"
    local elapsed=$(echo "$t_end - $t_start" | bc)
    echo "${label}: ${elapsed} s" >> "$PERF_LOG"
    echo "  Elapsed: ${elapsed} s"
}

function get_time() {
    python3 -c "import time; print(f'{time.time():.3f}')"
}

function check_executable() {
    if [ ! -f "$1" ]; then
        echo "ERROR: Executable not found: $1"
        echo "       Run ./setup.sh to build DFDM, or compile SPECFEM2D separately."
        exit 1
    fi
}

function print_profile_info() {
    echo ""
    echo "--- Precision Profile: ${PROFILE} ---"
    echo "  DFDM:   ppw=${DFDM_PPW}, order=${DFDM_ORDER}, gauss_order=${DFDM_GAUSS_ORDER}"
    echo "          receiver_elements=${DFDM_RECEIVER_ELEMS}"
    echo "  SPECFEM: time_scheme=${SPECFEM_TIME_SCHEME} (1=Newmark, 2=LDDRK4-6, 3=RK4)"
    echo "  DT=${DT}, NSTEP=${NSTEP}, total=${TOTAL_SIM_TIME}s"
    echo "  DFDM output:    ${DFDM_OUTPUT_DIR}"
    echo "  SPECFEM output: ${SPECFEM_OUTPUT_DIR}"
    echo ""
}

# macOS / Linux compatible sed -i
function sed_inplace() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "$@"
    else
        sed -i "$@"
    fi
}

# ===================== DFDM Config + Mesh =====================
function apply_dfdm_config() {
    # Update config.toml with current precision parameters
    log_step "Updating DFDM config (ppw=${DFDM_PPW}, order=${DFDM_ORDER}, gauss=${DFDM_GAUSS_ORDER}, dt=${DT})"

    sed_inplace "s/^time_steps = .*/time_steps = $NSTEP/" "$DFDM_CONFIG"
    sed_inplace "s/^NCPUs = .*/NCPUs = $DFDM_NPROC/" "$DFDM_CONFIG"
    sed_inplace "s/^delta_t = .*/delta_t = $DT/" "$DFDM_CONFIG"
    sed_inplace "s/^ppw = .*/ppw = $DFDM_PPW/" "$DFDM_CONFIG"
    sed_inplace "s/^order_b1 = .*/order_b1 = $DFDM_ORDER/" "$DFDM_CONFIG"
    sed_inplace "s/^order_b2 = .*/order_b2 = $DFDM_ORDER/" "$DFDM_CONFIG"
    sed_inplace "s/^gauss_order = .*/gauss_order = $DFDM_GAUSS_ORDER/" "$DFDM_CONFIG"
    sed_inplace "s/^receiver_elements = .*/receiver_elements = $DFDM_RECEIVER_ELEMS/" "$DFDM_CONFIG"

    echo "  Config updated: $DFDM_CONFIG"
    grep -E "^(ppw|order_b|gauss_order|delta_t|time_steps|NCPUs|receiver_elements)" "$DFDM_CONFIG"
}

function regenerate_dfdm_mesh() {
    # Regenerate mesh with current config (needed when ppw/order change)
    log_step "Regenerating DFDM mesh (ppw=${DFDM_PPW}, order=${DFDM_ORDER})"

    check_executable "$DFDM_MESH_EXEC"

    # Clean old mesh output and regenerate
    rm -rf "$DFDM_MESH_OUTPUT"
    mkdir -p "$DFDM_MESH_OUTPUT"

    "$DFDM_MESH_EXEC" "$DFDM_CONFIG" "$DFDM_MESH_OUTPUT" || {
        echo "ERROR: Mesh generation failed"
        exit 1
    }
    echo "  Mesh regenerated in $DFDM_MESH_OUTPUT"
}

# ===================== SPECFEM2D Config =====================
function apply_specfem_config() {
    log_step "Updating SPECFEM2D config (DT=${DT}, time_scheme=${SPECFEM_TIME_SCHEME}, NSTEP=${NSTEP})"

    sed_inplace "s/^NPROC[[:space:]]*=.*/NPROC                           = ${1:-$SPECFEM_NPROC}/" "$SPECFEM_PAR_FILE"
    sed_inplace "s/^NSTEP[[:space:]]*=.*/NSTEP                           = ${NSTEP}/" "$SPECFEM_PAR_FILE"
    sed_inplace "s/^DT[[:space:]]*=.*/DT                              = ${DT}/" "$SPECFEM_PAR_FILE"
    sed_inplace "s/^time_stepping_scheme[[:space:]]*=.*/time_stepping_scheme            = ${SPECFEM_TIME_SCHEME}/" "$SPECFEM_PAR_FILE"
    sed_inplace "s/^NTSTEP_BETWEEN_OUTPUT_IMAGES[[:space:]]*=.*/NTSTEP_BETWEEN_OUTPUT_IMAGES    = ${SPECFEM_SNAP_INTERVAL}/" "$SPECFEM_PAR_FILE"
    sed_inplace "s/^NTSTEP_BETWEEN_OUTPUT_SEISMOS[[:space:]]*=.*/NTSTEP_BETWEEN_OUTPUT_SEISMOS   = ${SPECFEM_SNAP_INTERVAL}/" "$SPECFEM_PAR_FILE"

    echo "  Par_file updated: $SPECFEM_PAR_FILE"
    grep -E "^(NPROC|NSTEP|DT |time_stepping_scheme|NTSTEP_BETWEEN_OUTPUT)" "$SPECFEM_PAR_FILE"
}

# ===================== Coordinate Sync =====================
# After DFDM runs, extract grid coordinates and update SPECFEM2D STATIONS
# so both codes measure at the same physical locations.
function sync_coordinates() {
    log_step "Syncing coordinates: DFDM grid -> SPECFEM2D STATIONS"

    if [ ! -f "$COORD_SCRIPT" ]; then
        echo "ERROR: Coordinate script not found: $COORD_SCRIPT"
        exit 1
    fi

    # Check that DFDM grid files exist
    if [ ! -f "${DFDM_OUTPUT_DIR}/grid_x_0" ]; then
        echo "ERROR: No DFDM grid files in ${DFDM_OUTPUT_DIR}. Run DFDM first."
        exit 1
    fi

    python3 "$COORD_SCRIPT" \
        --dfdm-data-dir "$DFDM_OUTPUT_DIR" \
        --config "$DFDM_CONFIG" \
        --stations-out "$SPECFEM_STATIONS" \
        --coords-json-out "$COORDS_JSON"

    echo "  STATIONS updated: $SPECFEM_STATIONS"
    echo "  Coords JSON:      $COORDS_JSON"
    echo "  STATIONS content:"
    cat "$SPECFEM_STATIONS" | while read line; do echo "    $line"; done
}

# ===================== DFDM =====================
function run_dfdm() {
    print_profile_info
    log_step "Running DFDM (${DFDM_NPROC} MPI tasks, ${NSTEP} steps, profile=${PROFILE})"

    check_executable "$DFDM_EXEC"

    # Apply precision config and regenerate mesh
    apply_dfdm_config
    regenerate_dfdm_mesh

    # Clean old output for this profile
    rm -rf "$DFDM_OUTPUT_DIR"
    mkdir -p "$DFDM_OUTPUT_DIR"
    mkdir -p "$DFDM_LOG_DIR"

    # Save config snapshot with output
    cp "$DFDM_CONFIG" "$DFDM_OUTPUT_DIR/config_used.toml"

    # Set up OpenBLAS if needed
    OPENBLAS_INSTALL="${DFDM_PROJECT_DIR}/openblas_install"
    if [ -d "$OPENBLAS_INSTALL" ]; then
        export CMAKE_PREFIX_PATH="${OPENBLAS_INSTALL}"
    fi

    echo "  Config: ${DFDM_CONFIG}"
    echo "  Output: ${DFDM_OUTPUT_DIR}"
    echo "  Steps:  ${NSTEP} (${TOTAL_TIME} s)"

    local t_start=$(get_time)

    # Compute relative output path from build dir
    local rel_output
    rel_output=$(python3 -c "import os; print(os.path.relpath('${DFDM_OUTPUT_DIR}', '${DFDM_BUILD_DIR}'))")

    cd "$DFDM_BUILD_DIR"
    mpirun -np $DFDM_NPROC ./dfdm \
        --config-file ../config/config.toml \
        --output-directory "./${rel_output}/" \
        --mesh-input-directory ../mesh_gen_ak135/build/output_test/ \
        --domain-output "./${rel_output}/domain/" \
        --log-files "./${rel_output}/domain/logs" \
        2>&1 | tee "${DFDM_OUTPUT_DIR}/dfdm_run.log"
    cd - > /dev/null

    local t_end=$(get_time)
    log_time "DFDM (${DFDM_NPROC} procs, ${NSTEP} steps, ${PROFILE})" "$t_start" "$t_end"

    # Verify output
    local n_snaps=$(ls -1 "$DFDM_OUTPUT_DIR"/elem_0_*.out 2>/dev/null | wc -l | tr -d ' ')
    local n_recv=$(ls -1 "$DFDM_OUTPUT_DIR"/receiver_elem_*_recorded_values.out 2>/dev/null | wc -l | tr -d ' ')
    echo "  Snapshots: ${n_snaps}, Receiver files: ${n_recv}"

    # Extract coordinates and update SPECFEM2D STATIONS to match DFDM grid
    sync_coordinates
}

# ===================== SPECFEM2D =====================
function run_specfem() {
    local nproc=${1:-$SPECFEM_NPROC}
    print_profile_info
    log_step "Running SPECFEM2D (${nproc} MPI tasks, ${NSTEP} steps, profile=${PROFILE})"

    check_executable "${SPECFEM_BIN}/xmeshfem2D"
    check_executable "${SPECFEM_BIN}/xspecfem2D"

    # Apply precision config
    apply_specfem_config "$nproc"

    # SPECFEM2D hardcodes output to OUTPUT_FILES/ relative to its working dir.
    # We always run in the default location, then copy to the profile directory.
    local SPECFEM_DEFAULT_OUTPUT="${SPECFEM_DIR}/OUTPUT_FILES"

    # Clean default output directory for fresh run
    rm -rf "${SPECFEM_DEFAULT_OUTPUT}"/*
    mkdir -p "$SPECFEM_DEFAULT_OUTPUT"

    # Copy config files to default output for reference
    cp "$SPECFEM_PAR_FILE" "$SPECFEM_DEFAULT_OUTPUT/"
    cp "${SPECFEM_DIR}/DATA/SOURCE" "$SPECFEM_DEFAULT_OUTPUT/"
    cp "${SPECFEM_DIR}/DATA/STATIONS" "$SPECFEM_DEFAULT_OUTPUT/" 2>/dev/null

    cd "$SPECFEM_DIR"

    echo "  Par_file: ${SPECFEM_PAR_FILE}"
    echo "  Output:   ${SPECFEM_OUTPUT_DIR}"
    echo "  Steps:    ${NSTEP} (${TOTAL_TIME} s)"

    # --- Mesher ---
    log_step "SPECFEM2D Mesher (NPROC=${nproc})"
    local t_mesh_start=$(get_time)
    if [ "$nproc" -eq 1 ]; then
        ./bin/xmeshfem2D 2>&1 | tee "${SPECFEM_DEFAULT_OUTPUT}/output_mesher.log"
    else
        mpirun -np $nproc ./bin/xmeshfem2D 2>&1 | tee "${SPECFEM_DEFAULT_OUTPUT}/output_mesher.log"
    fi
    local t_mesh_end=$(get_time)
    log_time "SPECFEM2D mesher (${nproc} procs)" "$t_mesh_start" "$t_mesh_end"

    # --- Solver ---
    log_step "SPECFEM2D Solver (NPROC=${nproc})"
    local t_solver_start=$(get_time)
    if [ "$nproc" -eq 1 ]; then
        ./bin/xspecfem2D 2>&1 | tee "${SPECFEM_DEFAULT_OUTPUT}/output_solver.log"
    else
        mpirun -np $nproc ./bin/xspecfem2D 2>&1 | tee "${SPECFEM_DEFAULT_OUTPUT}/output_solver.log"
    fi
    local t_solver_end=$(get_time)
    log_time "SPECFEM2D solver (${nproc} procs, ${NSTEP} steps, ${PROFILE})" "$t_solver_start" "$t_solver_end"

    cd - > /dev/null

    # For non-baseline profiles, copy output to the profile directory
    if [ "$SPECFEM_OUTPUT_DIR" != "$SPECFEM_DEFAULT_OUTPUT" ]; then
        log_step "Copying SPECFEM2D output to profile directory"
        rm -rf "$SPECFEM_OUTPUT_DIR"
        cp -r "$SPECFEM_DEFAULT_OUTPUT" "$SPECFEM_OUTPUT_DIR"
        echo "  Copied: $SPECFEM_DEFAULT_OUTPUT -> $SPECFEM_OUTPUT_DIR"
    fi

    # Verify output
    local n_seis=$(ls -1 "$SPECFEM_OUTPUT_DIR"/SY.*.se* 2>/dev/null | wc -l | tr -d ' ')
    local n_snaps=$(ls -1 "$SPECFEM_OUTPUT_DIR"/wavefield*.txt 2>/dev/null | wc -l | tr -d ' ')
    echo "  Seismograms: ${n_seis}, Wavefield dumps: ${n_snaps}"
}

# ===================== Comparison =====================
function run_comparison() {
    print_profile_info
    log_step "Running DFDM vs SPECFEM2D Comparison (profile=${PROFILE})"

    if [ ! -f "$COMPARE_SCRIPT" ]; then
        echo "ERROR: Comparison script not found: $COMPARE_SCRIPT"
        exit 1
    fi

    # Create figure output directory if needed
    mkdir -p "$FIGURE_OUTPUT_DIR"

    # Pass configuration to comparison script via environment variables
    export BENCH_DFDM_DATA_DIR="$DFDM_OUTPUT_DIR"
    export BENCH_SPECFEM_OUTPUT_DIR="$SPECFEM_OUTPUT_DIR"
    export BENCH_FIGURE_DIR="$FIGURE_OUTPUT_DIR"
    export BENCH_DT="$DT"
    export BENCH_NSTEP="$NSTEP"
    export BENCH_SPECFEM_DUMP_INTERVAL="$SPECFEM_SNAP_INTERVAL"
    export BENCH_COORDS_JSON="$COORDS_JSON"
    export BENCH_DFDM_CONFIG="$DFDM_CONFIG"
    export BENCH_TRACE_TAIL_TRIM="100.0"

    echo "  Compare script env:"
    echo "    BENCH_DFDM_DATA_DIR=$DFDM_OUTPUT_DIR"
    echo "    BENCH_SPECFEM_OUTPUT_DIR=$SPECFEM_OUTPUT_DIR"
    echo "    BENCH_FIGURE_DIR=$FIGURE_OUTPUT_DIR"
    echo "    BENCH_DT=$DT"
    echo "    BENCH_NSTEP=$NSTEP"
    echo "    BENCH_SPECFEM_DUMP_INTERVAL=$SPECFEM_SNAP_INTERVAL"
    echo "    BENCH_COORDS_JSON=$COORDS_JSON"
    echo "    BENCH_DFDM_CONFIG=$DFDM_CONFIG"
    echo "    BENCH_TRACE_TAIL_TRIM=100.0"

    local t_start=$(get_time)
    python3 "$COMPARE_SCRIPT" 2>&1 | tee "${FIGURE_OUTPUT_DIR}/comparison.log"
    local t_end=$(get_time)
    log_time "Comparison plots (${PROFILE})" "$t_start" "$t_end"

    local n_figs=$(ls -1 "$FIGURE_OUTPUT_DIR"/combined_t*.png 2>/dev/null | wc -l | tr -d ' ')
    echo "  Generated ${n_figs} comparison figures in ${FIGURE_OUTPUT_DIR}/"
}

# ===================== Performance Report =====================
function print_report() {
    log_step "Performance Summary (profile=${PROFILE})"

    if [ ! -f "$PERF_LOG" ]; then
        echo "No performance data found for profile '${PROFILE}'. Run the benchmark first."
        return
    fi

    echo ""
    echo "--- Timing Results ---"
    cat "$PERF_LOG"
    echo ""

    # Extract SPECFEM solver timing from its log
    local specfem_log="${SPECFEM_OUTPUT_DIR}/output_solver.log"
    if [ -f "$specfem_log" ]; then
        echo "--- SPECFEM2D Solver Details ---"
        grep -E "Total duration|Average duration|Total number of time steps" "$specfem_log" 2>/dev/null
        echo ""
    fi

    # Summary
    echo "--- Configuration ---"
    echo "  Profile:    ${PROFILE}"
    echo "  Domain:     1700 km uniform acoustic sphere"
    echo "  Vp:         2750 m/s"
    echo "  Density:    2000 kg/m^3"
    echo "  Source:     Ricker wavelet, f0=0.015 Hz"
    echo "  Time steps: ${NSTEP} (dt=${DT} s, total=${TOTAL_TIME} s)"
    echo "  DFDM:       ${DFDM_NPROC} MPI tasks, ppw=${DFDM_PPW}, order=${DFDM_ORDER}, gauss=${DFDM_GAUSS_ORDER}"
    echo "  SPECFEM2D:  ${SPECFEM_NPROC} MPI tasks, time_scheme=${SPECFEM_TIME_SCHEME}"
    echo ""
}

# ===================== Clean =====================
function clean_outputs() {
    log_step "Cleaning Output Files (profile=${PROFILE})"

    echo "Removing DFDM outputs: ${DFDM_OUTPUT_DIR}"
    rm -rf "$DFDM_OUTPUT_DIR"
    mkdir -p "$DFDM_OUTPUT_DIR"

    echo "Removing SPECFEM2D outputs: ${SPECFEM_OUTPUT_DIR}"
    rm -rf "${SPECFEM_OUTPUT_DIR}"
    mkdir -p "$SPECFEM_OUTPUT_DIR"

    echo "Removing performance log: ${PERF_LOG}"
    rm -f "$PERF_LOG"

    echo "Done."
}

# ===================== Full Benchmark =====================
function run_full_benchmark() {
    print_profile_info
    log_step "Full DFDM vs SPECFEM2D Benchmark (profile=${PROFILE})"

    # Initialize performance log
    echo "=== Benchmark: $(date) ===" > "$PERF_LOG"
    echo "Profile=${PROFILE}" >> "$PERF_LOG"
    echo "DFDM: ppw=${DFDM_PPW}, order=${DFDM_ORDER}, gauss=${DFDM_GAUSS_ORDER}" >> "$PERF_LOG"
    echo "SPECFEM: time_scheme=${SPECFEM_TIME_SCHEME}" >> "$PERF_LOG"
    echo "NSTEP=${NSTEP}, DT=${DT}, Total=${TOTAL_TIME} s" >> "$PERF_LOG"
    echo "" >> "$PERF_LOG"

    local t_total_start=$(get_time)

    run_dfdm
    run_specfem "$SPECFEM_NPROC"
    run_comparison

    local t_total_end=$(get_time)
    log_time "Total benchmark (${PROFILE})" "$t_total_start" "$t_total_end"

    print_report
}

# ===================== Main =====================
case "${1:-all}" in
    "all")
        run_full_benchmark
        ;;
    "dfdm")
        echo "=== DFDM Run (${PROFILE}): $(date) ===" >> "$PERF_LOG"
        run_dfdm
        ;;
    "specfem")
        echo "=== SPECFEM2D Run (${PROFILE}): $(date) ===" >> "$PERF_LOG"
        # Sync STATIONS from DFDM grid if grid files exist
        if [ -f "${DFDM_OUTPUT_DIR}/grid_x_0" ]; then
            sync_coordinates
        else
            echo "  WARNING: No DFDM grid files found. STATIONS not updated."
        fi
        run_specfem "$SPECFEM_NPROC"
        ;;
    "compare")
        echo "=== Comparison (${PROFILE}): $(date) ===" >> "$PERF_LOG"
        run_comparison
        ;;
    "report")
        print_report
        ;;
    "specfem_mesh")
        # Build SPECFEM2D external mesh: compile F90 generator, run it, cut to 1700km
        log_step "Generating SPECFEM2D mesh (${SPECFEM_DIR})"
        MESH_SRC="${SPECFEM_DIR}/mesh_gen/create_mesh_AK135F_2D_with_central_cube_no_PML_whole_earth.F90"
        MESH_BIN="${SPECFEM_BIN}/xcreate_mesh_files"
        CUT_SCRIPT="${SPECFEM_DIR}/mesh_gen/cut_mesh_1700km.py"

        if [ ! -f "$MESH_SRC" ]; then
            echo "ERROR: Mesh generator source not found: $MESH_SRC"
            exit 1
        fi

        echo "  Compiling mesh generator..."
        gfortran -O3 -o "$MESH_BIN" "$MESH_SRC"
        echo "  Compiled: $MESH_BIN"

        echo "  Running mesh generator (full-Earth mesh)..."
        cd "$SPECFEM_DIR"
        ./bin/xcreate_mesh_files
        cd - > /dev/null
        echo "  Full-Earth mesh generated in ${SPECFEM_DIR}/MESH/"

        if [ -f "$CUT_SCRIPT" ]; then
            echo "  Cutting mesh to 1700 km domain..."
            python3 "$CUT_SCRIPT"
            echo "  Cut mesh complete."
        else
            echo "  WARNING: cut_mesh_1700km.py not found, skipping mesh cut."
        fi
        ;;
    "clean")
        clean_outputs
        ;;
    *)
        echo "DFDM vs SPECFEM2D Benchmarking Script"
        echo ""
        echo "Usage: $0 [all|dfdm|specfem|compare|report|clean] [SPECFEM_NPROC]"
        echo ""
        echo "Actions:"
        echo "  all      - Run full benchmark: DFDM + SPECFEM2D + comparison (default)"
        echo "  dfdm     - Run DFDM solver only"
        echo "  specfem  - Run SPECFEM2D solver only (pass NPROC as 2nd arg, default=4)"
        echo "  compare  - Run comparison/plotting only (requires prior solver runs)"
        echo "  report   - Print performance summary from last benchmark"
        echo "  specfem_mesh - Build SPECFEM2D external mesh (compile generator + cut to 1700km)"
        echo "  clean    - Remove all output files for current profile"
        echo ""
        echo "Precision Profiles (set via PROFILE env var):"
        echo "  baseline    - Default settings (ppw=3, order=5, Newmark, dt=0.1)"
        echo "  highres     - Higher spatial accuracy (ppw=5, order=7, dt=0.1)"
        echo "  highres_dt  - Above + halved time step (ppw=5, order=7, dt=0.05)"
        echo "  reference   - DFDM baseline + SPECFEM2D finer mesh (mf=8)"
        echo ""
        echo "Environment Variable Overrides:"
        echo "  DFDM_PPW=5.0         DFDM points per wavelength"
        echo "  DFDM_ORDER=7         DFDM polynomial order (b1 and b2)"
        echo "  DFDM_GAUSS_ORDER=8   DFDM Gauss quadrature order"
        echo "  SPECFEM_TIME_SCHEME=2 SPECFEM time stepping (1=Newmark, 2=LDDRK, 3=RK4)"
        echo "  BENCH_DT=0.05        Time step for both codes"
        echo ""
        echo "Examples:"
        echo "  $0                                      # Full highres benchmark"
        echo "  PROFILE=baseline $0 all                 # Full baseline benchmark"
        echo "  PROFILE=highres $0 dfdm                 # DFDM only with highres"
        echo "  PROFILE=highres $0 specfem 1            # SPECFEM2D with highres, serial"
        echo "  PROFILE=highres $0 compare              # Compare highres results"
        echo "  DFDM_PPW=5 DFDM_ORDER=7 $0 dfdm        # Custom DFDM params"
        echo "  SPECFEM_TIME_SCHEME=2 $0 specfem 1      # SPECFEM with LDDRK4-6"
        echo ""
        echo "Current profile: ${PROFILE}"
        print_profile_info
        echo "Output directories:"
        echo "  DFDM:    ${DFDM_OUTPUT_DIR}"
        echo "  SPECFEM: ${SPECFEM_OUTPUT_DIR}"
        echo "  Figures: ${FIGURE_OUTPUT_DIR}"
        echo ""
        echo "Setup:"
        echo "  Run ./setup.sh to check/build prerequisites"
        exit 1
        ;;
esac
