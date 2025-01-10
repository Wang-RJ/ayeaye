# Create PSC for candidate sex chromosome analysis

## wc -l contigs_100kb.txt
# 398 contigs_100kb.txt

for i in {1..398};
do contig=$(head -$i BGI_contigs_100kb.txt | tail -1 | cut -f1)
echo $contig >> bcftools_PSC_bycontig.txt
bcftools stats -s - -r $contig 1/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.vcf.gz | grep 'PSC' >> bcftools_PSC_bycontig.txt
done

grep -v '#' bcftools_PSC_bycontig.txt | grep PSC > bcftools_PSC_bycontig.decap.txt