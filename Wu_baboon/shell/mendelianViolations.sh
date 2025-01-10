# Use bcftools plugin to create initial candidate list 
# !/bin/bash

bcftools annotate -Ou wu1/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.vcf.gz | bcftools +mendelian -T trio_tableMFC.csv -m x > ../mvfs/MV_bcfmendelian.vcf

for i in {1..19}
do
sample=$(head -n $i ../trio_tableMFC.csv | tail -n 1)
bcftools view -Ou -s $sample MV_bcfmendelian.vcf | bcftools +mendelian -t $sample -m x > MVF_trio$i.vcf
done

for i in {1..19}
do
bcftools view -M2 -m2 -i'GT[0]="RR" && GT[1]="RR" && GT[2]="het"' MVF_trio$i.vcf | grep -v "##" > MVF_trio$i.decap.vcf
done