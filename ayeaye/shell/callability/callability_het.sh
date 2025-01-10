# Grab all "true" heterozygotes in child:
# High quality genotypes in all trio individuals and where parents are homozygous for different alleles

for i in {1..12};
do
sample=$(head -n $i ../trio_tableMFC.csv | tail -n 1)
bcftools view -Ou -M2 -m2 -r $(tr '\n' ',' < candidate_contigs.txt) -s $sample ../vcf_filesBGI/8/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.vcf.gz | \
bcftools annotate -Oz -i'(GT[0]="RR" && GT[1]="AA") || (GT[0]="AA" && GT[1]="RR")' -o hets/hetcandidates.trio$i.vcf.gz &
done

# Create denominator by filtering parent genotypes to high quality sites and sampling roughly 250k
for i in {1..12}
do
bcftools view -Oz -i'FORMAT/DP[0] > 15 && FORMAT/DP[1] > 15 && FORMAT/DP[0] < 60 && FORMAT/DP[1] < 60 && FORMAT/GQ[0] > 60 && FORMAT/GQ[1] > 60' \
hets/hetcandidates.trio$i.vcf.gz -o hets/hetcandidates.filtered.trio$i.vcf.gz &
done

for i in {1..12}
do
lines=$(zcat hets/hetcandidates.filtered.trio$i.vcf.gz | grep -v '#' | wc -l)
zcat hets/hetcandidates.filtered.trio$i.vcf.gz | awk -v lines=$lines 'BEGIN {srand()} (rand() * lines < 250000 || /^#/ )' | bcftools view -Oz -o hets/hetdenom.trio$i.vcf.gz &
done

# Filtering children for transmission
for i in {1..12}
do
bcftools view -Oz -i'GT[2]="het"' hets/hetdenom.trio$i.vcf.gz -o hets/hettransmitted.trio$i.vcf.gz &
done

# Filtering children for depth
for i in {1..12}
do
child=$(cut -d',' -f3 ../trio_tableMFC.csv | head -n $i | tail -n 1)
max_idx=$(grep -n ^$child.dups.BGI.bam ../depthsBGI/bam_order | cut -d: -f1)
max_dp=$(head -n $max_idx ../depthsBGI/sums/max_depth_by_id.txt | tail -n 1)
bcftools view -Oz -i'FORMAT/DP[2] > 15 && FORMAT/DP[2] < '$max_dp hets/hettransmitted.trio$i.vcf.gz -o hets/hetxmitDP.trio$i.vcf.gz
done

# Create bcf mpileup of filtered sites
for trio in {1..12}
do
echo $trio
sbatch --job-name=cpile_$trio --export=i=$trio childpile_het.sbatch
done

for i in {1..12}
do
bcftools index hets/pileup/hetpileup.trio$i.vcf.gz &
bcftools index hets/hetxmitDP.trio$i.vcf.gz
done

for GQ in {20,30,40,50,60,70,80}
do
  for i in {1..12}
  do
    bcftools query hets/pileup/hetpileup.trio$i.vcf.gz -R <(bcftools query -i'FORMAT/GQ[2] > '$GQ hets/hetxmitDP.trio$i.vcf.gz -f '%CHROM\t%POS\n') -i'TYPE="snp" && FORMAT/ADF[:1] > 0 && FORMAT/ADR[:1] > 0' -f'%CHROM\t%POS\n' > hets/hetpiled_positionsGQ$GQ.trio$i.txt &
  done
  wait
done

# Build callability table for each trio
header="trio\tdenominator\ttransmit_filter\tdp_filter\tdpgq_filter\tbam_filter\tab_filtered"
for i in {1..12}; do bcftools index hets/hettransmitted.trio$i.vcf.gz; done
for i in {1..12}; do zgrep -v '#' hets/hetdenom.trio$i.vcf.gz | wc -l >> hets/columns/het_denom.column; done
for i in {1..12}; do zgrep -v '#' hets/hettransmitted.trio$i.vcf.gz | wc -l >> hets/columns/het_transmit.column; done
for i in {1..12}; do zgrep -v '#' hets/hetxmitDP.trio$i.vcf.gz | wc -l >> hets/columns/het_DP.column; done

for GQ in {20,30,40,50,60,70,80}
do
sbatch -J callGQ$GQ --export=GQ=$GQ het_dpgq_column.script
done

# table columns:
# denom, transmit, dp, dpgq, bam, ab30, ab35, ab40, ab45, ab50, ab55, ab60
for GQ in {20,30,40,50,60,70,80}
do
paste hets/columns/het_denom.column hets/columns/het_transmit.column hets/columns/het_DP.column hets/columns/het_dpgq$GQ.column \
hets/columns/het_bamgq$GQ.column hets/columns/het_ab30gq$GQ.column hets/columns/het_ab35gq$GQ.column \
hets/columns/het_ab40gq$GQ.column hets/columns/het_ab45gq$GQ.column hets/columns/het_ab50gq$GQ.column \
hets/columns/het_ab55gq$GQ.column hets/columns/het_ab60gq$GQ.column > het_callabilityGQ$GQ.table
done
