#! /bin/bash

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=cnstahl@princeton.edu

#SBATCH --job-name=test
#SBATCH --output=res.txt
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1 
#SBATCH  --mem=4000
#SBATCH --time=01:00:00

for N in 11 12 13 14 15
do 
    printf "\nComputing N = " $N
    python count_ground_states.py $N > ../data/ground\ state\ counts$N.txt
done