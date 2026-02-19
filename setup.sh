#!/bin/bash
# ==============================================================================
# DFDM Benchmark Setup Script
#
# Checks for pre-built binaries and builds them if missing.
# Prerequisites: cmake (>=3.22), MPI, GCC (with OpenMP), OpenBLAS
#
# Usage:
#   ./setup.sh              # Build everything that's missing
#   ./setup.sh --force      # Force rebuild even if binaries exist
#   ./setup.sh --check      # Only check, don't build
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DFDM_DIR="${SCRIPT_DIR}/DFDM_2D_1.0"
DFDM_BUILD_DIR="${DFDM_DIR}/build"
DFDM_EXEC="${DFDM_BUILD_DIR}/dfdm"
MESH_GEN_DIR="${DFDM_DIR}/mesh_gen_ak135"
MESH_GEN_BUILD_DIR="${MESH_GEN_DIR}/build"
MESH_GEN_EXEC="${MESH_GEN_BUILD_DIR}/mesh_gen"
OPENBLAS_INSTALL="${DFDM_DIR}/openblas_install"

FORCE=false
CHECK_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --force) FORCE=true ;;
        --check) CHECK_ONLY=true ;;
    esac
done

function log_step() {
    echo ""
    echo "=================================================================="
    echo "  $1"
    echo "=================================================================="
}

function check_prerequisites() {
    log_step "Checking prerequisites"
    local missing=()

    command -v cmake >/dev/null 2>&1 || missing+=("cmake (>=3.22)")
    command -v mpirun >/dev/null 2>&1 || missing+=("MPI (mpirun)")
    command -v python3 >/dev/null 2>&1 || missing+=("python3")

    if [ ${#missing[@]} -gt 0 ]; then
        echo "ERROR: Missing prerequisites:"
        for m in "${missing[@]}"; do
            echo "  - $m"
        done
        echo ""
        echo "Install on macOS:  brew install cmake open-mpi gcc openblas python3"
        echo "Install on Ubuntu: apt install cmake libopenmpi-dev gfortran libopenblas-dev python3"
        return 1
    fi

    echo "  All prerequisites found."
    return 0
}

function check_python_packages() {
    log_step "Checking Python packages"
    local missing=()

    python3 -c "import numpy" 2>/dev/null || missing+=("numpy")
    python3 -c "import matplotlib" 2>/dev/null || missing+=("matplotlib")
    python3 -c "import scipy" 2>/dev/null || missing+=("scipy")

    if [ ${#missing[@]} -gt 0 ]; then
        echo "  Missing Python packages: ${missing[*]}"
        echo "  Install with: pip3 install -r ${SCRIPT_DIR}/requirements.txt"
        return 1
    fi

    echo "  All Python packages found."
    return 0
}

function build_dfdm() {
    log_step "Building DFDM solver"

    if [ ! -d "$DFDM_DIR/src" ]; then
        echo "ERROR: DFDM source not found at ${DFDM_DIR}/src"
        echo "       If this is a fresh clone, ensure DFDM_2D_1.0/ is present."
        exit 1
    fi

    mkdir -p "$DFDM_BUILD_DIR"

    local cmake_args="-DCMAKE_BUILD_TYPE=Release"

    # Use bundled OpenBLAS install if available
    if [ -d "$OPENBLAS_INSTALL" ]; then
        cmake_args="${cmake_args} -DCMAKE_PREFIX_PATH=${OPENBLAS_INSTALL}"
        echo "  Using bundled OpenBLAS: ${OPENBLAS_INSTALL}"
    fi

    cd "$DFDM_BUILD_DIR"
    cmake .. ${cmake_args}
    make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
    cd - > /dev/null

    if [ -f "$DFDM_EXEC" ]; then
        echo "  DFDM built successfully: ${DFDM_EXEC}"
    else
        echo "ERROR: DFDM build failed"
        exit 1
    fi
}

function build_mesh_gen() {
    log_step "Building DFDM mesh generator"

    if [ ! -d "$MESH_GEN_DIR/src" ]; then
        echo "ERROR: Mesh generator source not found at ${MESH_GEN_DIR}/src"
        exit 1
    fi

    mkdir -p "$MESH_GEN_BUILD_DIR"

    local cmake_args="-DCMAKE_BUILD_TYPE=Release"

    if [ -d "$OPENBLAS_INSTALL" ]; then
        cmake_args="${cmake_args} -DCMAKE_PREFIX_PATH=${OPENBLAS_INSTALL}"
    fi

    cd "$MESH_GEN_BUILD_DIR"
    cmake .. ${cmake_args}
    make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
    cd - > /dev/null

    if [ -f "$MESH_GEN_EXEC" ]; then
        echo "  Mesh generator built successfully: ${MESH_GEN_EXEC}"
    else
        echo "ERROR: Mesh generator build failed"
        exit 1
    fi
}

# ===================== Main =====================

check_prerequisites
check_python_packages || true  # warn but don't fail

echo ""
echo "--- Binary Status ---"

DFDM_OK=false
MESH_OK=false
SPECFEM_MESH_OK=false
SPECFEM_SOLVE_OK=false

[ -f "$DFDM_EXEC" ] && DFDM_OK=true
[ -f "$MESH_GEN_EXEC" ] && MESH_OK=true
[ -f "${SCRIPT_DIR}/specfem2d/bin/xmeshfem2D" ] && SPECFEM_MESH_OK=true
[ -f "${SCRIPT_DIR}/specfem2d/bin/xspecfem2D" ] && SPECFEM_SOLVE_OK=true

echo "  DFDM solver:     $([ "$DFDM_OK" = true ] && echo 'OK' || echo 'MISSING')"
echo "  DFDM mesh_gen:   $([ "$MESH_OK" = true ] && echo 'OK' || echo 'MISSING')"
echo "  SPECFEM mesher:  $([ "$SPECFEM_MESH_OK" = true ] && echo 'OK' || echo 'MISSING (must compile SPECFEM2D separately)')"
echo "  SPECFEM solver:  $([ "$SPECFEM_SOLVE_OK" = true ] && echo 'OK' || echo 'MISSING (must compile SPECFEM2D separately)')"

if [ "$CHECK_ONLY" = true ]; then
    echo ""
    echo "Check-only mode. Exiting."
    exit 0
fi

# Build DFDM if missing or forced
if [ "$FORCE" = true ] || [ "$DFDM_OK" = false ]; then
    build_dfdm
else
    echo "  DFDM solver already built. Use --force to rebuild."
fi

if [ "$FORCE" = true ] || [ "$MESH_OK" = false ]; then
    build_mesh_gen
else
    echo "  Mesh generator already built. Use --force to rebuild."
fi

if [ "$SPECFEM_MESH_OK" = false ] || [ "$SPECFEM_SOLVE_OK" = false ]; then
    echo ""
    echo "NOTE: SPECFEM2D binaries are not present."
    echo "  To use SPECFEM2D, compile it separately and copy the binaries:"
    echo "    cp /path/to/specfem2d/bin/xmeshfem2D ${SCRIPT_DIR}/specfem2d/bin/"
    echo "    cp /path/to/specfem2d/bin/xspecfem2D ${SCRIPT_DIR}/specfem2d/bin/"
    echo "  (Also copy to specfem2d_mf8/bin/ if using highres profile)"
fi

echo ""
echo "Setup complete."
