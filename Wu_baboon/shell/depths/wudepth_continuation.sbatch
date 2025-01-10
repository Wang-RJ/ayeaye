#!/bin/bash

#SBATCH -p general

#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --chdir=/N/project/Wuetal_baboon/depths/

module unload gcc
module unload python
module load gcc/9.3.0
module load python/3.9.8
module load samtools

samtools depth -r chr2:188193158-188283584 ../bam_files/*.bam > samtools_depth.chr2_continuation.txt &
samtools depth -r chr3:162523916-181338541 ../bam_files/*.bam > samtools_depth.chr3_continuation.txt &
samtools depth -r chr4:144873507-166074386 ../bam_files/*.bam > samtools_depth.chr4_continuation.txt

