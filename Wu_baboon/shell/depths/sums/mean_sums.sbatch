#!/bin/bash

#SBATCH -p general

#SBATCH --ntasks=20
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --chdir=/N/project/Wuetal_baboon/depths/sums/

for i in {1..20}
do
echo Running: chr$i
./mean_sums.sh ../samtools_depth.chr$i.txt > dpcounts.chr$i.txt &
done
wait
echo Done
