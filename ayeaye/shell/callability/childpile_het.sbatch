#!/bin/bash

#SBATCH -p general

#SBATCH --ntasks=1
#SBATCH --time=01:30:00
#SBATCH --mem=2G
#SBATCH --chdir=/N/project/ayeaye/callability/

module unload gcc
module unload python
module load gcc/9.3.0
module load python/3.9.8
module load bcftools

childbam=$(cut -d',' -f3 ../trio_tableMFC.csv | head -n $i | tail -n 1).dups.BGI.bam

bcftools mpileup -Oz -A -f /N/project//ayeaye/ref/GCA_023783475.1_ASM2378347v1_genomic.fna -Q 10 \
-R <(bcftools query hets/hetxmitDP.trio$i.vcf.gz -f'%CHROM\t%POS\n') \
-a ADF,ADR,AD,DP /N/project/ayeaye/bam_files/$childbam > hets/pileup/hetpileup.trio$i.vcf.gz

