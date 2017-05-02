#! /bin/bash

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=cnstahl@princeton.edu

#SBATCH --job-name=test
#SBATCH --output=res.txt
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1 
#SBATCH  --mem=5000
#SBATCH --time=24:00:00

srun python entropy.py