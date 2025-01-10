#!/bin/bash

#SBATCH -p general

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:00:00
#SBATCH --mem=1G
#SBATCH --chdir=/N/project/Wuetal_baboon/depths/

module unload gcc
module unload python
module load gcc/9.3.0
module load python/3.9.8
module load samtools

samtools depth -r $chr ../bam_files/*.bam > samtools_depth.$chr.txt

