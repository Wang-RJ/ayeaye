# Use bcftools plugin to create initial candidate list 
# !/bin/bash

bcftools annotate -Ou 1/output.filtered.snps.CalculateGenotypePosteriors.deNovo.vcf.gz | bcftools +mendelian -T ../trio_tableMFC.csv -m x > ../mvfs/MV_bcfmendelian.trio1.vcf

# trios 3 and 6 have mixed up parentage, will need to drop in downstream analysis

for i in {1..6}
do
sample=$(head -n $i ../trio_tableMFC.csv | tail -n 1)
bcftools view -Ou -s $sample MV_bcfmendelian.trio1.vcf | bcftools +mendelian -t $sample -m x > MVF_trio$i.vcf
done

for i in {1,2,4,5}
do
bcftools view -M2 -m2 -i'GT[0]="RR" && GT[1]="RR" && GT[2]="het"' MVF_trio$i.vcf | grep -v "##" > MVF_trio$i.decap.vcf
done