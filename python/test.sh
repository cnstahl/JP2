#! /bin/bash

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=cnstahl@princeton.edu

#SBATCH --job-name=test
#SBATCH --output=res.txt
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1 
#SBATCH  --mem=40
#SBATCH --time=0:01:00

srun python test.py