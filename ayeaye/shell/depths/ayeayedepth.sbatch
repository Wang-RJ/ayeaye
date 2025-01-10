#!/bin/bash

#SBATCH -p general

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:00:00
#SBATCH --mem=3G
#SBATCH --chdir=/N/project/ayeaye/depths/

module unload gcc
module unload python
module load gcc/9.3.0
module load python/3.9.8
module load samtools

contigs=$(cat batch_lists/depth_batch$batch)
echo Contigs:$contigs

for contig in $contigs
do
samtools depth -r $contig ../bam_files/*.bam >> samtools_depth.batch$batch.txt
done
