#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob_hybrid.out.%j
#SBATCH -e ./tjob_hybrid.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J test_slurm
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# for OpenMP:
#SBATCH --cpus-per-task=40
#
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#
# Wall clock limit:
#SBATCH --time=00:30:00

# Load compiler and MPI modules with explicit version specifications,
# consistently with the versions used to build the executable.
module purge
module load intel/21.2.0 impi/2021.2

# Run the program:
srun ./build_dir/panda < inp_sc >log.log

