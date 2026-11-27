#!/bin/bash
# run_mpi_cuda.sh: Slurm job for the MPI+CUDA and CUDA codes (strong scaling)
#   Settings are for the AiMOS DCS cluster (IBM Power9, 4 V100 GPUs per node); change the
#   partition and modules for another cluster.
#
# Usage (from src/, after "make POWER9=1 cuda mpi-cuda"):
#   mkdir -p results
#   sbatch run_mpi_cuda.sh
# Output: results/mpi_cuda_<jobid>.out, results/mpi_cuda_err_<jobid>.out

#SBATCH --job-name=poisson2d-mpi-cuda
#SBATCH -N 4                    # 4 compute nodes
#SBATCH -n 16                   # 16 ranks: one per GPU
#SBATCH --partition=dcs-2024
#SBATCH --time=00:05:00
#SBATCH --gres=gpu:4
#SBATCH --output=results/mpi_cuda_%j.out
#SBATCH --error=results/mpi_cuda_err_%j.out

cd "$SLURM_SUBMIT_DIR"

module purge
module load xl_r spectrum-mpi cuda/11.2

ARGS="2.0 1e-10 20000"          # omega tol max_iter
RUN="mpirun --bind-to core"

echo "Starting at $(date)"
echo "Running on hosts: $SLURM_NODELIST ($SLURM_NNODES nodes), job ID: $SLURM_JOB_ID"

echo "----------------------- STRONG SCALING -----------------------------"
for N in 256 512 1024 2048 4096
do
    echo "=============== STRONG SCALING FOR N = $N ==============="

    # -N 4: 4 ranks per node (one per GPU)
    echo "=== NP = 16 (4 nodes) ==="
    $RUN -N 4 -np 16 ./poisson2d_mpi_cuda $N $ARGS

    echo "=== NP = 8 (2 nodes) ==="
    $RUN -N 4 -np 8 ./poisson2d_mpi_cuda $N $ARGS

    echo "=== NP = 4 and 2 (1 node) ==="
    $RUN -np 4 ./poisson2d_mpi_cuda $N $ARGS
    $RUN -np 2 ./poisson2d_mpi_cuda $N $ARGS

    echo "=== NP = 1 (CUDA only) ==="
    ./poisson2d_cuda $N $ARGS
done

echo "Finished at $(date)"
