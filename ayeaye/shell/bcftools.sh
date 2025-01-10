#!/bin/bash

# Run bcftools mpileup to compare with GATK calls at Mendelian violation candidates

bcftools mpileup -Ou -A -f ../ref/GCA_023783475.1_ASM2378347v1_genomic.fna -Q 10 -R cpos_DP_GQ.txt -a ADF,ADR,AD,DP bam_files/*.bam | \
bcftools call -m -Ov -o bcfpileup/bcfpileup_ayeaye.vcf

grep -v '##' bcfpileup_ayeaye.vcf > bcfpileup_ayeaye.decap.vcf

# Run bcftools to extract GATK genotypes for each trio
bcftools 1/output.filtered.snps.CalculateGenotypePosteriors.deNovo.vcf.gz -R cpos_DP_GQ.txt > gatk_candidates_DPGQ.trio1.vcf