### 2D usntructurd mesh solver 

You are currently on `unstructured` branch, this branch forked from the structured 2D mesh solver code. The solver on this branch accepts unstructured mesh as input and simulates accoustic wavepropagation using DFD Method. Below are simple instructions to do a simple run.

### building the solver
DFDM depends on MPI and OpenBLAS, MPI is available by default on Perlmutter and a copy of OpenBLAS is available in Muaaz' account. To build DFDM on Perlmutter follow the below instructions, if you are building on your personal computer then jump ahead to section `building in environments other perlmutter`.

```bash
git clone --single-branch --branch unstructured https://github.com/mgawan/DFDM_2D_1.0.git
cd DFDM_2D_1.0
export DFDM_ROOT=${PWD}
git submodule init
git submodule update

source config/env_perlmutter.sh

mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ../
make -j4
```
The above will build the solver, then before doing a test run you will need to generate a mesh. Follow below steps to build and run the DFDM mesh generator.

### CLI Options
The solver uses the following command-line flags:

- `--config-file` (required): Path to the TOML simulation config.
- `--output-directory` (required): Directory for all simulation outputs.
- `--mesh-input-directory` (required): Directory containing the generated mesh.
- `--domain-output` (required): Directory to store a copy of domain decomposition.
- `--log-files` (optional): Directory for per-rank log files.
- `--print-debug` (optional): `true|false`; when `true`, detailed logs are written to per-rank files.
- `--enable-tests-run` (optional): When `true`, the solver exits before the time step loop (useful for CI or quick checks).

### Building Mesh Generator
Clone and build the mesh_gen code:

```bash
cd ${DFDM_ROOT}
mkdir mesh_gen
cd mesh_gen
git clone https://github.com/mgawan/mesh_gen.git .
git submodule init
git submodule update

mkdir build 
cd build
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ../
make -j4
```
### Test Simulation run
Once the mesh generator has been built successfully, do the below run for generating simple mesh that will be input to the DFD Solver:

```bash
mkdir ${DFDM_ROOT}/mesh_input/
cd ${DFDM_ROOT}/mesh_gen/build
./mesh_gen ${DFDM_ROOT}/config/config.toml ${DFDM_ROOT}/mesh_input/
```

The above will generate a mesh in `${DFDM_ROOT}/mesh_input/` directory which can then be used for simulation:

for running the solver first load the cpu resources:


```bash
salloc -N 1 -C cpu -A <account> -t 30 -q interactive
```

```bash
cd ${DFDM_ROOT}/build
mkdir sample_out
mkdir -p sample_out/domain/logs
srun -n <CPUS> ./dfdm \
    --config-file ../config/config.toml \
    --output-directory ./sample_out/ \
    --mesh-input-directory ${DFDM_ROOT}/mesh_input/ \
    --domain-output ./sample_out/domain/ \
    --log-files ./sample_out/domain/logs \
    --print-debug false
```
In the above line, CPUS must be replaced by the number of CPUs listed in the config.toml file, CPU numbers must match in the config.toml file because distribution of elements depends on the number of CPUs.

The above will produce a large number of output files for the simulation in the directory `./sample_out/` this directory can be changed by the users.

Note: To perform an early exit (no time-step loop) for testing or CI, append `--enable-tests-run true` to the `srun` command.

### Building and Running on environments other than Perlmutter
Before installing DFDM, OpenBLAS and MPI needs to be available. MPI can be installed using package manager of your choice, for OpenBLAS follow below instructions:

```bash
git clone --single-branch --branch unstructured https://github.com/mgawan/DFDM_2D_1.0.git
cd DFDM_2D_1.0
export DFDM_ROOT=${PWD}
git submodule init
git submodule update

mkdir ${DFDM_ROOT}/openblas_install/
cd ${DFDM_ROOT}
git clone https://github.com/OpenMathLib/OpenBLAS.git
export BLAS_DIR=${DFDM_ROOT}/OpenBLAS/
export BLAS_INSTALL=${DFDM_ROOT}/openblas_install/
cd ${BLAS_DIR}
make USE_OPENMP=1 NO_SHARED=0 MAKE_NB_JOBS=12  CC="gcc" FC="gfortran" DYNAMIC_ARCH=1 COMMON_OPT="-O3 -fPIC -pthread" FCOMMON_OPT="-O3 -fPIC -pthread" LDFLAGS="-lpthread -fopenmp"
make PREFIX ${BLAS_INSTALL} install
export OPEN_BLAS_CMAKE_PATH=${BLAS_INSTALL}/lib/cmake

cd ${DFDM_ROOT}

mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ../
make -j4
```

#### Building Mesh Generator
Clone and build the mesh_gen code:

```bash
cd ${DFDM_ROOT}
mkdir mesh_gen
cd mesh_gen
git clone https://github.com/mgawan/mesh_gen.git .
git submodule init
git submodule update

mkdir build 
cd build
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release ../
make -j4
```
#### Test Simulation run
Once the mesh generator has been built successfully, do the below run for generating simple mesh that will be input to the DFD Solver:

```bash
mkdir ${DFDM_ROOT}/mesh_input/
cd ${DFDM_ROOT}/mesh_gen/build
./mesh_gen ${DFDM_ROOT}/config/config.toml ${DFDM_ROOT}/mesh_input/
```

The above will generate a mesh in `${DFDM_ROOT}/mesh_input/` directory which can then be used for simulation:

for running the solver first load the cpu resources:

```bash
salloc -N 1 -C cpu -A <account> -t 30 -q interactive
```

```bash
cd ${DFDM_ROOT}/build
mkdir sample_out
srun -n <CPUS> ./dfdm \
    --config-file ../config/config.toml \
    --output-directory ./sample_out/ \
    --mesh-input-directory ${DFDM_ROOT}/mesh_input/ \
    --domain-output ./sample_out/domain/ \
    --log-files ./sample_out/domain/logs \
    --print-debug false
```
In the above line, CPUS must be replaced by the number of CPUs listed in the config.toml file, CPU numbers must match in the config.toml file because distribution of elements depends on the number of CPUs.

The above will produce a large number of output files for the simulation in the directory `./sample_out/` this directory can be changed by the users.

Note: To perform an early exit (no time-step loop) for testing or CI, append `--enable-tests-run true` to the `srun` command.

# Visualizing Output using Visualizer.py

The `visualizer.py` script is designed to plot 2D unstructured data from simulation outputs. It can visualize grids, element IDs, and time-step data.

---

## Usage

Run the script from the command line with the required arguments:

```bash
python visualizer.py --data_dir <path_to_data_directory> --total_elem <number_of_elements> --time_steps <number_of_time_steps> [options]
```
### Visualizer.py Documentation

The `visualizer.py` script is a Python-based tool designed to visualize the output data generated by the DFDM solver. It supports plotting 2D unstructured mesh data, including grids, element IDs, and time-step data.

---

#### Features

- Visualize 2D unstructured mesh grids.
- Display element IDs for debugging and analysis.
- Plot simulation data across multiple time steps.
- Configurable options for customizing the visualization.

---

#### Command-Line Arguments

The script accepts the following arguments:

| Argument            | Description                                                                 | Required | Default Value |
|---------------------|-----------------------------------------------------------------------------|----------|---------------|
| `--data_dir`        | Path to the directory containing simulation output data.                   | Yes      | N/A           |
| `--total_elem`      | Total number of elements in the simulation.                                | Yes      | N/A           |
| `--time_steps`      | Number of time steps to visualize.                                         | Yes      | N/A           |
| `--plot_grid`       | Flag to enable grid visualization.                                         | No       | False         |
| `--plot_elem_ids`   | Flag to enable plotting of element IDs.                                    | No       | False         |

---

#### Example Usage

1. **Basic Visualization**:
    ```bash
    python visualizer.py --data_dir ./sample_out --total_elem 9 --time_steps 10000
    ```

2. **Visualize Grid and Element IDs**:
    ```bash
    python visualizer.py --data_dir ./sample_out --total_elem 9 --time_steps 10000 --plot_grid --plot_elem_ids
    ```

---

#### Dependencies

Ensure the following Python libraries are installed before running the script:

- `matplotlib`
- `numpy`
- `argparse`

Install them using pip if not already available:

```bash
pip install matplotlib numpy
```

---

#### Notes

- The `--data_dir` must contain the simulation output files in the expected format.
- Ensure the `--total_elem` and `--time_steps` values match the simulation configuration.
- For large datasets, consider increasing system memory or reducing the number of time steps to visualize.








## Testing
Currently a series of tests are being added to ensure continuous development and testing of this solver. To run the existing tests `BUILD_TESTS` option in CMake can be turned on, this will build the tests along with enabling the test code.

```bash
cmake -DBUILD_TESTS=ON ../
```

After building as usual, `ctest` can be called to run all the tests. It must be noted that the parallel simulation test will require a node allocation when running tests on a slurm equipped machine.

You can also run the solver with an early-exit for quick checks or CI pipelines by adding the flag:

```bash
srun -n <CPUS> ./dfdm \
    --config-file ../config/config.toml \
    --output-directory ./sample_out/ \
    --mesh-input-directory ${DFDM_ROOT}/mesh_input/ \
    --domain-output ./sample_out/domain/ \
    --enable-tests-run true
```

### Next Steps:

1. Ensure multiple ranks are not competing for the same file. Directory creation, IO handling,
2. Right now, all the elements for each rank are defined in main and then each domain has its element numbers, the element list is passed around.
3. Add the rotation field in elements.
4. add the conforming boundary condition and T matrix.
5. Element neighbors from Matlab code.
6. CPU neighbors, figure out how to do this. 
7. checkpoint: domain.cpp line 83, element.cpp line 3
