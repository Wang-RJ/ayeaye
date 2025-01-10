#!/bin/bash

#SBATCH -p general

#SBATCH --ntasks=20
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --chdir=/N/project/Wuetal_baboon/depths/

for i in {1..20}
do
echo Running: chr$i
./count_tbodepth_wu.sh samtools_depth.chr$i.txt > counts/tbocounts.chr$i.txt &
done
wait
echo Done

