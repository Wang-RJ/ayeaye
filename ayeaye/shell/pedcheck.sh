# Generate sample for pedigree checking

bcftools view -r $(tr '\n' ',' < ../candidate_contigs.txt) -i 'MIN(FMT/DP)>=20 & MAX(FMT/DP)<=60 & MIN(FMT/GQ)>=20 & GT[*]!="mis"' 1/output.filtered.snps.CalculateGenotypePosteriors.deNovo.vcf.gz -Oz -o ayeaye.pedcheck.preproc.vcf.gz

shuf -n 100000 <(zcat ayeaye.pedcheck.preproc.vcf.gz | grep -v '#')> ayeaye100k_shufSNPs.vcf

# zcat ayeaye.pedcheck.preproc.vcf.gz | grep -n '#CHROM'
## line 2751

cat <(zcat ayeaye.pedcheck.preproc.vcf.gz | head -n 2751 | tail -n 1) ayeaye100k_shufSNPs.vcf > ayeaye100k_shufSNPs.header.vcf
mv ayeaye100k_shufSNPs.header.vcf ayeaye100k_shufSNPs.vcf
