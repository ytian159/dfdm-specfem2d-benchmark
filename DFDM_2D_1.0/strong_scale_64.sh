#!/bin/bash

#SBATCH -C cpu
#SBATCH -N 1
#SBATCH -A m4661
#SBATCH -t 90
#SBATCH -q regular
#SBATCH -o output_%j.txt

echo "Running on host $(hostname) with $SLURM_JOB_NUM_NODES nodes"
echo "Running on node $(scontrol show hostnames $SLURM_JOB_NODELIST)"

NB=8
CPUS=64

echo "Running strong scaling test for NB=${NB} and CPUS=${CPUS}"

DFDM_BUILD_DIR=/pscratch/sd/m/mgawan/dfdm_scaling_studies/solver/build
INPUT_DATA_DIR=/pscratch/sd/m/mgawan/dfdm_scaling_studies/input_data
STRONG_SCALE_DIR=${INPUT_DATA_DIR}/strong_scaling
INPUT_MESH_DIR=${STRONG_SCALE_DIR}/NB_${NB}/CPUS_${CPUS}
CONFIG_FILE=${INPUT_MESH_DIR}/config_file_${CPUS}.toml
OUTPUT_DIR=/pscratch/sd/m/mgawan/dfdm_scaling_studies/jobs/NB_${NB}/CPUS_${CPUS}/output
DOMAIN_OUTPUT_DIR=${OUTPUT_DIR}/domain
LOG_FILES_DIR=${OUTPUT_DIR}/logs
PRINT_DEBUG=false

if [ "$PRINT_DEBUG" = true ]; then
  echo "Debugging trace is enabled"
else
  echo "Debugging trace is disabled"
fi

echo "Entering build directory: ${DFDM_BUILD_DIR}"

cd ${DFDM_BUILD_DIR}
source ../config/env_perlmutter.sh

TOTAL_TASKS=${CPUS}
TASKS_PER_NODE=$((TOTAL_TASKS / SLURM_JOB_NUM_NODES))
CPUS_PER_TASK=$((256 / TASKS_PER_NODE))



set -x
srun --cpus-per-task=${CPUS_PER_TASK} \
     --tasks-per-node=${TASKS_PER_NODE} \
     --cpu_bind=cores \
     ./dfdm \
     --config-file ${CONFIG_FILE} \
     --output-directory ${OUTPUT_DIR}/ \
     --mesh-input-directory ${INPUT_MESH_DIR}/ \
     --domain-output ${DOMAIN_OUTPUT_DIR}/ \
     --log-files ${LOG_FILES_DIR}/ \
     --print-debug ${PRINT_DEBUG}
