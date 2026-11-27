#!/bin/bash

#Slurm Directives
#SBATCH --job-name=poisson2d_mpi
#SBATCH -N 4				# Request 4 compute nodes
#SBATCH -n 128				# Reserve 128 cores
#SBATCH --partition=dcs-2024
#SBATCH --time=00:05:00          
#SBATCH --gres=gpu:4          
#SBATCH --output=/gpfs/u/home/PCPG/PCPGlngn/scratch/project/results/mpi_%j.out
#SBATCH --error=/gpfs/u/home/PCPG/PCPGlngn/scratch/project/results/mpi_err_%j.out


# Change directory
cd /gpfs/u/home/PCPG/PCPGlngn/scratch/project

# Clear the environment from any previously loaded modules
module purge

# Load the module environment suitable for the job
module load xl_r spectrum-mpi cuda/11.2

# Debugging information
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
# echo "Running $SLURM_NTASKS tasks."
# echo "Current working directory is `pwd`"
echo "Job ID: $SLURM_JOB_ID"

# mpirun --bind-to core --report-bindings -np $SLURM_NTASKS ./poisson2d_mpi 1024 2.0 1e-10 20000


# --------------------------- STRONG SCALING -----------------------------

echo "----------------------- STRONG SCALING  -----------------------------"
echo "---------------------------------------------------------------------"

for N in 128 256 512 1024
do
    echo "=============================================================="
    echo "RUNNING STRONG SCALING FOR N = $N"
    echo "=============================================================="

    # 4 nodes
    echo "=== NP = 128 (4 Nodes) ==="
    mpirun --bind-to core --report-bindings -N 32 -np 128 \
        ./poisson2d_mpi $N 2.0 1e-10 20000

    # 2 nodes
    echo "=== NP = 64 (2 Nodes) ==="
    mpirun --bind-to core --report-bindings -N 32 -np 64 \
        ./poisson2d_mpi $N 2.0 1e-10 20000

    # 1 node scaling
    for nprocs in 32 16 8 4 2
    do
        echo "------------------------------------------"
        echo "RUNNING WITH NP = $nprocs, N = $N"
        echo "------------------------------------------"

        mpirun --bind-to core --report-bindings -np $nprocs \
            ./poisson2d_mpi $N 2.0 1e-10 20000

        echo "Finished NP = $nprocs, N = $N at $(date)"
    done

done

# UNCOMMENT HERE


# echo "----------------------- STRONG SCALING  -----------------------------"
# echo "---------------------------------------------------------------------"
# N=1024
# # RUN ON FOUR NODES
# echo "=== Running Multi-Node: NP = 128 (4 Nodes) ==="
# # We specify -N 32 to ensure it spreads across 4 nodes (32 per node)
# mpirun --bind-to core --report-bindings -N 32 -np 128 ./poisson2d_mpi 1024 2.0 1e-10 20000


# # RUN ON TWO NODES
# echo "=== Running Multi-Node: NP = 64 (2 Nodes) ==="
# # We specify -N 32 to ensure it spreads across 2 nodes (32 per node)
# mpirun --bind-to core --report-bindings -N 32 -np 64 ./poisson2d_mpi 1024 2.0 1e-10 20000


# # RUN ONE A SINGLE NODES
# # Loop through process counts
# for nprocs in 32 16 8 4 2
# do
# 	echo "------------------------------------------"
#     echo "RUNNING WITH NP = $nprocs"
#     echo "------------------------------------------"

#     # Run the job
#     mpirun --bind-to core --report-bindings -np $nprocs ./poisson2d_mpi 1024 2.0 1e-10 20000
#     echo "Finished NP = $nprocs at $(date)"

# done




# mpirun --bind-to core --report-bindings -np 2 ./poisson2d_mpi 1024 2.0 1e-10 20000
# mpirun --bind-to core --report-bindings -np 4 ./poisson2d_mpi 1024 2.0 1e-10 20000
# mpirun --bind-to core --report-bindings -np 8 ./poisson2d_mpi 1024 2.0 1e-10 20000
# mpirun --bind-to core --report-bindings -np 16 ./poisson2d_mpi 1024 2.0 1e-10 20000
# mpirun --bind-to core --report-bindings -np 32 ./poisson2d_mpi 1024 2.0 1e-10 20000


# --------------------------- WEAK SCALING -----------------------------
echo "----------------------- RUNNING WEAK SCALING  -------------------------"

# mpirun --bind-to core --report-bindings -N 32 -np 128 ./poisson2d_mpi 2048 2.0 1e-10 20000
# mpirun --bind-to core --report-bindings -N 32 -np 128 ./poisson2d_mpi 2048 2.0 1e-10 20000
# mpirun --bind-to core --report-bindings -N 8 -np 32 ./poisson2d_mpi 1024 2.0 1e-10 20000
# mpirun --bind-to core --report-bindings -N 2 -np 8 ./poisson2d_mpi 512 2.0 1e-10 20000
mpirun --bind-to core --report-bindings -N 2 -np 2 ./poisson2d_mpi 256 2.0 1e-10 20000




