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

# ===================== Platform Detection =====================
if [ -n "$NERSC_HOST" ] || hostname | grep -qE '(login|nid).*\.perlmutter'; then
    ON_PERLMUTTER=true
else
    ON_PERLMUTTER=false
fi

# Load modules on Perlmutter
if [ "$ON_PERLMUTTER" = true ]; then
    module load cmake PrgEnv-gnu cray-mpich python 2>/dev/null || true
    export OPEN_BLAS_CMAKE_PATH=/pscratch/sd/m/mgawan/openblas_install/OpenBLAS/install/lib/cmake
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-16}
    CXX_COMPILER="CC"  # Cray wrapper (includes MPI)
else
    CXX_COMPILER=""  # let cmake auto-detect
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DFDM_DIR="${SCRIPT_DIR}/DFDM_2D_1.0"
DFDM_BUILD_DIR="${DFDM_DIR}/build"
DFDM_EXEC="${DFDM_BUILD_DIR}/dfdm"
MESH_GEN_DIR="${DFDM_DIR}/mesh_gen_ak135"
MESH_GEN_BUILD_DIR="${MESH_GEN_DIR}/build"
MESH_GEN_EXEC="${MESH_GEN_BUILD_DIR}/mesh_gen"
OPENBLAS_INSTALL="${DFDM_DIR}/openblas_install"
SPECFEM2D_SRC_DIR="${SCRIPT_DIR}/specfem2d_src"

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
    if [ "$ON_PERLMUTTER" = true ]; then
        command -v srun >/dev/null 2>&1 || missing+=("srun (SLURM)")
        command -v CC >/dev/null 2>&1 || missing+=("Cray C++ wrapper (CC)")
    else
        command -v mpirun >/dev/null 2>&1 || missing+=("MPI (mpirun)")
    fi
    command -v python3 >/dev/null 2>&1 || missing+=("python3")

    if [ ${#missing[@]} -gt 0 ]; then
        echo "ERROR: Missing prerequisites:"
        for m in "${missing[@]}"; do
            echo "  - $m"
        done
        echo ""
        echo "Install on macOS:       brew install cmake open-mpi gcc openblas python3"
        echo "Install on Ubuntu:      apt install cmake libopenmpi-dev gfortran libopenblas-dev python3"
        echo "On NERSC Perlmutter:    module load cmake PrgEnv-gnu cray-mpich"
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

    # On Perlmutter, use Cray C++ wrapper and system OpenBLAS
    if [ -n "$CXX_COMPILER" ]; then
        cmake_args="${cmake_args} -DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
    fi
    if [ -n "$OPEN_BLAS_CMAKE_PATH" ]; then
        cmake_args="${cmake_args} -DCMAKE_PREFIX_PATH=${OPEN_BLAS_CMAKE_PATH}"
        echo "  Using OpenBLAS: ${OPEN_BLAS_CMAKE_PATH}"
    elif [ -d "$OPENBLAS_INSTALL" ]; then
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

    if [ -n "$CXX_COMPILER" ]; then
        cmake_args="${cmake_args} -DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
    fi
    if [ -n "$OPEN_BLAS_CMAKE_PATH" ]; then
        cmake_args="${cmake_args} -DCMAKE_PREFIX_PATH=${OPEN_BLAS_CMAKE_PATH}"
    elif [ -d "$OPENBLAS_INSTALL" ]; then
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

function clone_mesh_gen() {
    log_step "Cloning DFDM mesh generator (mesh_gen_ak135)"

    local MESH_GEN_REPO="https://github.com/mgawan/mesh_gen.git"
    git clone -b ak135 "$MESH_GEN_REPO" "$MESH_GEN_DIR" || {
        echo "ERROR: Failed to clone mesh_gen. Trying SSH..."
        git clone -b ak135 "git@github.com:mgawan/mesh_gen.git" "$MESH_GEN_DIR" || {
            echo "ERROR: Failed to clone mesh_gen via SSH"
            exit 1
        }
    }

    cd "$MESH_GEN_DIR"
    git submodule init && git submodule update
    cd - > /dev/null

    echo "  Mesh generator cloned: ${MESH_GEN_DIR}"
}

function check_specfem_platform() {
    # Check if SPECFEM2D binaries match the current platform
    local binary="${SCRIPT_DIR}/specfem2d/bin/xspecfem2D"
    if [ ! -f "$binary" ]; then
        return 1  # missing
    fi
    # Check if the binary is for this platform
    local fmt
    fmt=$(file "$binary" 2>/dev/null)
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64)
            echo "$fmt" | grep -q "ELF 64-bit.*x86-64" && return 0 ;;
        Darwin-arm64)
            echo "$fmt" | grep -q "Mach-O.*arm64" && return 0 ;;
        Darwin-x86_64)
            echo "$fmt" | grep -q "Mach-O.*x86_64" && return 0 ;;
    esac
    return 1  # wrong platform
}

function build_specfem() {
    log_step "Building SPECFEM2D from source"

    if [ ! -d "$SPECFEM2D_SRC_DIR" ]; then
        echo "  Cloning SPECFEM2D source..."
        git clone "https://github.com/SPECFEM/specfem2d.git" "$SPECFEM2D_SRC_DIR" 2>/dev/null || \
        git clone "git@github.com:SPECFEM/specfem2d.git" "$SPECFEM2D_SRC_DIR" || {
            echo "ERROR: Failed to clone SPECFEM2D"
            exit 1
        }
    fi

    cd "$SPECFEM2D_SRC_DIR"

    if [ "$ON_PERLMUTTER" = true ]; then
        MPI_INC="$CRAY_MPICH_DIR/include" \
        MPI_LIBS="$CRAY_MPICH_DIR/lib" \
        ./configure FC=ftn CC=cc CXX=CC MPIFC=ftn --with-mpi
    else
        ./configure FC=gfortran MPIFC=mpif90
    fi

    make all

    # Copy binaries to benchmark directories
    for dest in "${SCRIPT_DIR}/specfem2d/bin" "${SCRIPT_DIR}/specfem2d_mf8/bin"; do
        mkdir -p "$dest"
        cp bin/xmeshfem2D "$dest/"
        cp bin/xspecfem2D "$dest/"
    done

    cd - > /dev/null
    echo "  SPECFEM2D built and installed."
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

# Clone mesh_gen_ak135 if missing
if [ ! -d "$MESH_GEN_DIR/src" ]; then
    if [ "$CHECK_ONLY" = true ]; then
        echo "  Mesh generator source: MISSING (run setup without --check to clone)"
    else
        clone_mesh_gen
    fi
fi

if [ "$FORCE" = true ] || [ "$MESH_OK" = false ]; then
    if [ "$CHECK_ONLY" != true ]; then
        build_mesh_gen
    fi
else
    echo "  Mesh generator already built. Use --force to rebuild."
fi

# Build SPECFEM2D if binaries are missing or wrong platform
if ! check_specfem_platform || [ "$FORCE" = true ]; then
    if [ "$CHECK_ONLY" = true ]; then
        if [ "$SPECFEM_MESH_OK" = false ] || [ "$SPECFEM_SOLVE_OK" = false ]; then
            echo "  SPECFEM2D binaries: MISSING"
        else
            echo "  SPECFEM2D binaries: WRONG PLATFORM (run setup without --check to rebuild)"
        fi
    else
        build_specfem
    fi
else
    echo "  SPECFEM2D binaries already present for this platform."
fi

echo ""
echo "Setup complete."
