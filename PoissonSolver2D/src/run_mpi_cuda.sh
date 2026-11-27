#!/bin/bash

#Slurm Directives
#SBATCH --job-name=poisson2d-mpi-cuda
#SBATCH -N 4				# Request 4 compute nodes
#SBATCH -n 16				# Reserve 16 cores
#SBATCH --partition=dcs-2024
#SBATCH --time=00:5:00         
#SBATCH --gres=gpu:4          
#SBATCH --output=/gpfs/u/home/PCPG/PCPGlngn/scratch/project/results/mpi_cuda_%j.out
#SBATCH --error=/gpfs/u/home/PCPG/PCPGlngn/scratch/project/results/mpi_cuda_err_%j.out


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


echo "----------------------- STRONG SCALING  -----------------------------"
echo "---------------------------------------------------------------------"

# ------------------------------------------------------------------------------- #
for N in 256 512 1024 2048 4096
do
	echo "=============================================================="
    echo "RUNNING STRONG SCALING FOR N = $N"
    echo "=============================================================="

	# RUN ON FOUR NODES
	echo "========== Running Multi-Node: NP = 16 (4 Nodes) =========="
	# We specify -N 4 to ensure it spreads across 4 nodes (4 processes per node)
	mpirun --bind-to core -N 4 -np 16 ./poisson2d_mpi_cuda $N 2.0 1e-10 20000


	# RUN ON TWO NODES
	echo "========== Running Multi-Node: NP = 8 (2 Nodes) =========="
	# We specify -N 4 to ensure it spreads across 2 nodes (4 processes per node)
	mpirun --bind-to core -N 4 -np 8 ./poisson2d_mpi_cuda $N 2.0 1e-10 20000


	# RUN ONE A SINGLE NODES
	echo "========== Running Sigle-Node: NP = 4 and 2 (1 Node) =========="
	mpirun --bind-to core -np 4 ./poisson2d_mpi_cuda $N 2.0 1e-10 20000
	mpirun --bind-to core -np 2 ./poisson2d_mpi_cuda $N 2.0 1e-10 20000

	# Just run CUDA version only
	echo "========== Running Sigle-Node: NP = 1 (CUDA ONLY) =========="
	./poisson2d_cuda $N 2.0 1e-10 20000

done

