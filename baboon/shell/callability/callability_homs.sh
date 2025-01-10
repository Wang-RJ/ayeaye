# Parent callability assessed separately, callability from parents is a square of individual probabilities of passing homozygous filter

for i in {1..15}
do
sample=$(head -n $i ../trio_tableMFC.csv | tail -n 1)
bcftools view -Ou -M2 -m2 -r $(tr '\n' ',' < ../autosomes.list) -i'INFO/AC < 2' -s $sample ../vcf_files/wu19/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.vcf.gz | \
bcftools annotate -x INFO/InbreedingCoeff,INFO/BaseQRankSum,INFO/ExcessHet,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQRankSum,INFO/ReadPosRankSum,INFO/SOR \
-i'(GT[0]="RR" && GT[1]="RR")' | awk 'BEGIN {srand()} (rand() <= 0.5 || /^#/ )' | bcftools view -Oz -o homs/homcandidates.trio$i.vcf.gz &
done

for i in {16..19}
do
sample=$(head -n $i ../trio_tableMFC.csv | tail -n 1)
bcftools view -Ou -M2 -m2 -r $(tr '\n' ',' < ../autosomes.list) -i'INFO/AC < 2' -s $sample ../vcf_files/wu1/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.vcf.gz | \
bcftools annotate -x INFO/InbreedingCoeff,INFO/BaseQRankSum,INFO/ExcessHet,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQRankSum,INFO/ReadPosRankSum,INFO/SOR \
-i'(GT[0]="RR" && GT[1]="RR")' | awk 'BEGIN {srand()} (rand() <= 0.5 || /^#/ )' | bcftools view -Oz -o homs/homcandidates.trio$i.vcf.gz &
done

# Filter to high quality sites
for i in {1..19}
do
bcftools view -Oz -i'FORMAT/DP[0] > 15 && FORMAT/DP[1] > 15 && FORMAT/DP[0] < 60 && FORMAT/DP[1] < 60 && FORMAT/GQ[0] > 60 && FORMAT/GQ[1] > 60' \
homs/homcandidates.trio$i.vcf.gz -o homs/homcandidates.filtered.trio$i.vcf.gz &
done

for i in {1..19}
do
lines=$(zcat homs/homcandidates.filtered.trio$i.vcf.gz | grep -v '#' | wc -l)
zcat homs/homcandidates.filtered.trio$i.vcf.gz | awk -v lines=$lines 'BEGIN {srand()} (rand() * lines < 250000 || /^#/ )' | \
bcftools view -Oz -o homs/homdenom.trio$i.vcf.gz &
done

# Filtering child for transmission
for i in {1..19}
do
bcftools view -Oz -i'GT[2]="RR"' homs/homdenom.trio$i.vcf.gz -o homs/homtransmitted.trio$i.vcf.gz &
done

# Filtering child for depth
for i in {1..19}
do
child=$(cut -d',' -f3 ../trio_tableMFC.csv | head -n $i | tail -n 1)
max_idx=$(grep -n ^$child.dups.bam ../depths/sums/bam_order | cut -d: -f1)
max_dp=$(head -n $max_idx ../depths/sums/max_depth_by_id.txt | tail -n 1)
bcftools view -Oz -i'FORMAT/DP[2] > 15 && FORMAT/DP[2] < '$max_dp homs/homtransmitted.trio$i.vcf.gz -o homs/homxmitDP.trio$i.vcf.gz &
done

# Create pileup of filtered sites

### Deprecated
# for i in {1..19}
# do
# childbam=$(cut -d',' -f3 ../trio_tableMFC.csv | head -n $i | tail -n 1).bam
# bcftools mpileup -Oz -A -f /N/project/baboon/reference/GCA_000264685.2_Panu_3.0_genomic.fna -Q 10 \
# -R <(bcftools query homs/homxmitDP.trio$i.vcf.gz -f'%CHROM\t%POS\n') \
# -a ADF,ADR,AD,DP /N/project/baboon/bam_files/$childbam > homs/pileup/hompileup.trio$i.vcf.gz &
# done

for trio in {1..19}
do
echo $trio
sbatch --job-name=cpile_$trio --export=i=$trio childpile_homs.sbatch
done

for i in {1..19}
do
bcftools index homs/pileup/hompileup.trio$i.vcf.gz
bcftools index homs/homxmitDP.trio$i.vcf.gz &
done

for GQ in {20,30,40,50,60,70,80}
do
  for i in {1..19}
  do
    bcftools query homs/pileup/hompileup.trio$i.vcf.gz -R <(bcftools query -i'FORMAT/GQ[2] > '$GQ homs/homxmitDP.trio$i.vcf.gz -f '%CHROM\t%POS\n') \
-e'TYPE="snp" && FORMAT/AD[:1] > 0' -f'%CHROM\t%POS\n' > homs/hompiledAD0_positionsGQ$GQ.trio$i.txt &
  done
  wait
done

# Build callability table for each trio
header="trio\tdenominator\ttransmit_filter\tdpgq_filter\tbam_filter\tab_filtered"
for i in {1..19}; do bcftools index homs/homtransmitted.trio$i.vcf.gz; done
for i in {1..19}; do zgrep -v '#' homs/homdenom.trio$i.vcf.gz | wc -l >> homs/columns/hom_denom.column; done
for i in {1..19}; do zgrep -v '#' homs/homtransmitted.trio$i.vcf.gz | wc -l >> homs/columns/hom_transmit.column; done
for i in {1..19}; do zgrep -v '#' homs/homxmitDP.trio$i.vcf.gz | wc -l >> homs/columns/hom_DP.column; done

for GQ in {20,30,40,50,60,70,80}
do
  for i in {1..19}
  do
  bcftools view -i'FORMAT/GQ[2] > '$GQ homs/homxmitDP.trio$i.vcf.gz | grep -v '#' | wc -l >> homs/columns/hom_dpgq$GQ.column
  wc -l < homs/hompiledAD0_positionsGQ$GQ.trio$i.txt >> homs/columns/hom_AD0_bamgq$GQ.column
  done
done

# denom, transmit, dp, dpgq, bam
for GQ in {20,30,40,50,60,70,80}
do
paste homs/columns/hom_denom.column homs/columns/hom_transmit.column homs/columns/hom_DP.column homs/columns/hom_dpgq$GQ.column homs/columns/hom_AD0_bamgq$GQ.column > hom_callabilityGQ$GQ.table
done