#!/bin/bash

#SBATCH -p general

#SBATCH --ntasks=20
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --chdir=/N/project/ayeaye/depths/sums/

for i in {0..19}
do
echo Running: batch$i
./mean_sums.sh ../samtools_depth.batch$i.txt > dpcounts.batch$i.txt &
done
wait
echo Done
