#!/bin/bash

# Run bcftools mpileup to compare with GATK calls at Mendelian violation candidates
bcftools mpileup -Ou -A -f baboon/reference/GCA_000264685.2_Panu_3.0_genomic.fna -Q 10 -R $file -a ADF,ADR,AD,DP ../bam_files/*.bam | \
bcftools call -m -Ov -o ../bcfpileup/bcfpileup_wubaboon_$file.vcf

# Run bcftools to extract GATK genotypes for each trio
bcftools annotate wu1/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.vcf.gz -R ../cpos_DP_GQ.txt > gatk_candidates_DPGQ.trio1.vcf