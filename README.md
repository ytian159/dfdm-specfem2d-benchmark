# DFDM vs SPECFEM2D Benchmark

Standalone benchmark comparing [DFDM](https://github.com/ytian159/DFDM_2D_1.0) and [SPECFEM2D](https://github.com/SPECFEM/specfem2d) solvers on a 1700 km uniform acoustic domain (Vp=2750 m/s, ρ=2000 kg/m³, f₀=0.015 Hz).

## Quick Start

```bash
# Clone with submodules
git clone --recursive https://github.com/ytian159/dfdm-specfem2d-benchmark.git
cd dfdm-specfem2d-benchmark

# Install Python dependencies
pip3 install -r requirements.txt

# Check prerequisites and build DFDM if needed
./setup.sh

# Run full benchmark (DFDM + SPECFEM2D + comparison)
./run_benchmark.sh
```

## Repository Structure

```
dfdm-specfem2d-benchmark/
├── DFDM_2D_1.0/              # DFDM solver (git submodule → ytian159/DFDM_2D_1.0)
│   ├── src/, include/        # C++ source code
│   ├── config/               # Configuration files
│   ├── mesh_gen_ak135/       # Mesh generator
│   ├── build/                # Compiled binaries (if built)
│   └── openblas_install/     # Bundled OpenBLAS (optional)
├── specfem2d/                # SPECFEM2D setup (mf=5, baseline mesh)
│   ├── bin/                  # Pre-compiled binaries
│   ├── DATA/                 # Par_file, SOURCE, STATIONS
│   ├── MESH/                 # 1700 km mesh files
│   └── mesh_gen/             # Mesh generation tools
├── specfem2d_mf8/            # High-resolution SPECFEM2D (mf=8)
├── scripts/                  # Python analysis scripts
│   ├── compare_snapshots.py
│   └── compute_dfdm_coordinates.py
├── run_benchmark.sh          # Main benchmark script
├── setup.sh                  # Build/check prerequisites
└── requirements.txt          # Python dependencies
```

## Prerequisites

- **Compilers**: cmake (≥3.22), MPI (OpenMPI), GCC (with OpenMP), gfortran
- **Libraries**: OpenBLAS (or system BLAS)
- **Python 3**: numpy, matplotlib, scipy
- **Tools**: git

### Installation (macOS)
```bash
brew install cmake open-mpi gcc openblas python3
pip3 install -r requirements.txt
```

### Installation (Ubuntu/Debian)
```bash
apt install cmake libopenmpi-dev gfortran libopenblas-dev python3 python3-pip
pip3 install -r requirements.txt
```

## Usage

### 1. Setup (first time only)
```bash
./setup.sh
```
Checks prerequisites and builds DFDM if binaries are missing.

### 2. Run Benchmark

**Full benchmark (default: highres profile)**
```bash
./run_benchmark.sh
```

**Run individual solvers**
```bash
./run_benchmark.sh dfdm              # DFDM only
./run_benchmark.sh specfem 4         # SPECFEM2D only (4 MPI tasks)
./run_benchmark.sh compare           # Comparison plots only
```

**Different precision profiles**
```bash
PROFILE=baseline ./run_benchmark.sh all      # ppw=3, order=5, dt=0.1
PROFILE=highres ./run_benchmark.sh all       # ppw=5, order=7, dt=0.1 (default)
PROFILE=highres_dt ./run_benchmark.sh all    # ppw=5, order=7, dt=0.05
```

**Custom parameters**
```bash
DFDM_PPW=5 DFDM_ORDER=7 ./run_benchmark.sh dfdm
BENCH_DT=0.05 ./run_benchmark.sh all
```

### 3. View Results

**Performance summary**
```bash
./run_benchmark.sh report
```

**Output locations**
- DFDM: `DFDM_2D_1.0/build/sample_out_<profile>/`
- SPECFEM2D: `specfem2d/OUTPUT_FILES_<profile>/`
- Comparison plots: Same as SPECFEM2D output directory

## Precision Profiles

| Profile      | DFDM (ppw, order) | SPECFEM2D (mesh, scheme) | dt (s) |
|--------------|-------------------|--------------------------|--------|
| `baseline`   | 3.0, 5            | mf=5, Newmark            | 0.1    |
| `highres`    | 5.0, 7            | mf=8, Newmark            | 0.1    |
| `highres_dt` | 5.0, 7            | mf=5, Newmark            | 0.05   |
| `reference`  | 3.0, 5            | mf=8, Newmark            | 0.1    |

Override via environment variables: `DFDM_PPW`, `DFDM_ORDER`, `DFDM_GAUSS_ORDER`, `SPECFEM_TIME_SCHEME`, `BENCH_DT`

## Working with DFDM Source Code

The DFDM source is included as a git submodule pointing to your fork ([ytian159/DFDM_2D_1.0](https://github.com/ytian159/DFDM_2D_1.0)), which tracks the upstream [mgawan/DFDM_2D_1.0](https://github.com/mgawan/DFDM_2D_1.0).

### Making changes to DFDM

```bash
cd DFDM_2D_1.0

# Make changes to source code
vim src/simulation.cpp

# Commit and push to your fork
git add src/simulation.cpp
git commit -m "Update simulation logic"
git push origin main

# Update benchmark repo to track new commit
cd ..
git add DFDM_2D_1.0
git commit -m "Update DFDM submodule to latest"
git push
```

### Submitting changes upstream

1. Push changes to your fork (`ytian159/DFDM_2D_1.0`)
2. Create a pull request on GitHub from your fork to `mgawan/DFDM_2D_1.0`

### Updating DFDM from upstream

```bash
cd DFDM_2D_1.0
git remote add upstream https://github.com/mgawan/DFDM_2D_1.0.git
git fetch upstream
git merge upstream/main
git push origin main

cd ..
git add DFDM_2D_1.0
git commit -m "Update DFDM to upstream main"
git push
```

## SPECFEM2D Binaries

Pre-compiled SPECFEM2D binaries are included for convenience. To recompile (if needed on a different platform):

1. Clone and build SPECFEM2D separately
2. Copy binaries to `specfem2d/bin/` and `specfem2d_mf8/bin/`

```bash
# Example
git clone https://github.com/SPECFEM/specfem2d.git /tmp/specfem2d
cd /tmp/specfem2d
./configure FC=gfortran MPIFC=mpif90
make all

# Copy to benchmark repo
cp bin/xmeshfem2D /path/to/dfdm_benchmark/specfem2d/bin/
cp bin/xspecfem2D /path/to/dfdm_benchmark/specfem2d/bin/
cp bin/xmeshfem2D /path/to/dfdm_benchmark/specfem2d_mf8/bin/
cp bin/xspecfem2D /path/to/dfdm_benchmark/specfem2d_mf8/bin/
```

## Output Files

**DFDM**
- Receiver waveforms: `receiver_elem_N_recorded_values.out`
- Wavefield snapshots: `elem_N_STEP.out`
- Grid coordinates: `grid_x_N`, `grid_z_N`

**SPECFEM2D**
- Seismograms: `SY.RECNN.PRE.semp`, `SY.RECNN.BXX.semd`, etc.
- Wavefield dumps: `wavefield*.txt`
- Grid definition: `wavefield_grid_for_dumps.txt`

**Comparison**
- Combined snapshot plots: `combined_tNNNN.png`
- Waveform comparisons: Overlaid in snapshot plots
- Logs: `comparison.log`

## Troubleshooting

**"Executable not found" errors**
```bash
./setup.sh              # Build missing binaries
./setup.sh --force      # Force rebuild
```

**"OpenBLAS not found" during DFDM build**
```bash
# Use system OpenBLAS
export CMAKE_PREFIX_PATH=/path/to/openblas

# Or use Homebrew OpenBLAS (macOS)
export CMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
```

**SPECFEM2D segfault**
- Check MPI is properly installed (`mpirun --version`)
- Ensure binaries are compiled for your platform (see "SPECFEM2D Binaries" above)

**Coordinate mismatch warnings**
- Normal: SPECFEM2D uses bilinear interpolation from DFDM's GLL grid
- Documented in `compute_dfdm_coordinates.py`

## References

- **DFDM**: Yuan, Y.O., Simons, F.J., & Tromp, J. (2020). "Double-difference adjoint seismic tomography". *Geophysical Journal International*, 206(3), 1599-1618.
- **SPECFEM2D**: Tromp, J., Komatitsch, D., & Liu, Q. (2008). "Spectral-element and adjoint methods in seismology". *Communications in Computational Physics*, 3(1), 1-32.

## License

- DFDM: See [ytian159/DFDM_2D_1.0](https://github.com/ytian159/DFDM_2D_1.0)
- SPECFEM2D: GPL v3 (see [SPECFEM/specfem2d](https://github.com/SPECFEM/specfem2d))
- This benchmark: MIT (see LICENSE)

## Contact

Yuan Tian - [@ytian159](https://github.com/ytian159)

Repository: [https://github.com/ytian159/dfdm-specfem2d-benchmark](https://github.com/ytian159/dfdm-specfem2d-benchmark)
