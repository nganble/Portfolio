#!/bin/bash
# run_mpi.sh: Slurm job for the MPI code (strong and weak scaling)
#   Settings are for the AiMOS DCS cluster (IBM Power9, 32 cores per node); change the
#   partition and modules for another cluster.
#
# Usage (from src/, after "make POWER9=1 mpi"):
#   mkdir -p results
#   sbatch run_mpi.sh
# Output: results/mpi_<jobid>.out, results/mpi_err_<jobid>.out

#SBATCH --job-name=poisson2d_mpi
#SBATCH -N 4                    # 4 compute nodes
#SBATCH -n 128                  # 128 cores
#SBATCH --partition=dcs-2024
#SBATCH --time=00:05:00
#SBATCH --gres=gpu:4
#SBATCH --output=results/mpi_%j.out
#SBATCH --error=results/mpi_err_%j.out

cd "$SLURM_SUBMIT_DIR"

module purge
module load xl_r spectrum-mpi cuda/11.2

ARGS="2.0 1e-10 20000"          # omega tol max_iter
RUN="mpirun --bind-to core --report-bindings"

echo "Starting at $(date)"
echo "Running on hosts: $SLURM_NODELIST ($SLURM_NNODES nodes), job ID: $SLURM_JOB_ID"

# --------------------------- STRONG SCALING -----------------------------
echo "----------------------- STRONG SCALING -----------------------------"
for N in 128 256 512 1024
do
    echo "=============== STRONG SCALING FOR N = $N ==============="

    echo "=== NP = 128 (4 nodes) ==="
    $RUN -N 32 -np 128 ./poisson2d_mpi $N $ARGS

    echo "=== NP = 64 (2 nodes) ==="
    $RUN -N 32 -np 64 ./poisson2d_mpi $N $ARGS

    for nprocs in 32 16 8 4 2
    do
        echo "=== NP = $nprocs (1 node) ==="
        $RUN -np $nprocs ./poisson2d_mpi $N $ARGS
    done
done

# --------------------------- WEAK SCALING -----------------------------
# N^2 / np is constant: (N, np) = (256, 2), (512, 8), (1024, 32), (2048, 128)
echo "----------------------- WEAK SCALING -------------------------"
$RUN -N 2  -np 2   ./poisson2d_mpi 256  $ARGS
$RUN -N 2  -np 8   ./poisson2d_mpi 512  $ARGS
$RUN -N 8  -np 32  ./poisson2d_mpi 1024 $ARGS
$RUN -N 32 -np 128 ./poisson2d_mpi 2048 $ARGS

echo "Finished at $(date)"
