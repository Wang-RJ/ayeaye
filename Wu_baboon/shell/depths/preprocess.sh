# Index bam files
sbatch --job-name=index_10173.dups.bam --export=bamfile=10173.dups.bam bamindex.sbatch

for bamfile in $(ls *.bam | tail -n +12)
do
echo $bamfile
sbatch --job-name=idx_$bamfile --export=bamfile=$bamfile bamindex.sbatch
done

# Run samtools depth
sbatch --job-name=dp_chr1 --export=chr=chr1 wudepth.sbatch
for i in {2..20}
do
echo chr$i
sbatch --job-name=dp_chr$i --export=chr=chr$i wudepth.sbatch
done

# Didn't finish depth for 2,3,4 in 8 hours. Paste on continuation after finding breakpoint
cat <(head -n -1 samtools_depth.chr2.txt) <(tail -n +2 samtools_depth.chr2_continuation.txt) > samtools_depth.chr2.tmp
cat <(head -n -1 samtools_depth.chr3.txt) <(tail -n +2 samtools_depth.chr3_continuation.txt) > samtools_depth.chr3.tmp
cat <(head -n -1 samtools_depth.chr4.txt) <(tail -n +2 samtools_depth.chr4_continuation.txt) > samtools_depth.chr4.tmp

# Make tables for mean coverage
for file in $(ls dpcounts.chr* | sort -V)
do
head -n 1 $file >> sums_by_chr.txt
tail -n +2 $file | head -n 1 >> sqsums_by_chr.txt
tail -n 1 $file >> totpos_by_chr.txt
done

awk '{for(i=1;i<=NF;i++) {sum[i] += $i}} END {for(i=1;i<=NF;i++) print sum[i]}' sums_by_chr.txt > sums_by_id.txt
paste -d/ sums_by_id.txt <(yes $(paste -sd+ totpos_by_chr.txt | bc) | head -n 11) | bc -l > mean_depth_by_id.txt
awk '{print $1 + 4 *sqrt($1)}' mean_depth_by_id.txt > max_depth_by_id.txt

# Paste tbo counts together after running tbo_wu.sbatch which calls count_tbo
paste tbocounts* | tr ' ' '\t' | awk '{for(i=2;i<=NF;i+=2) { sum[NR] += $i } print sum[NR]}'

# bcfpileup slurm call

sbatch --job-name=cposaa --export=file=cposaa.txt bcfpileup/bcfpileup.sbatch

for file in $(ls cpos*.txt | tail -n +2)
do
echo $file
sbatch --job-name=pile_$file --export=file=$file ../bcfpileup/bcfpileup.sbatch
done

for file in $(tr '_' '\t' < ../bcfpileup/incomplete_1hr | cut -f2)
do
echo $file
sbatch --job-name=pile_$file --export=file=$file ../bcfpileup/bcfpileup.sbatch
done

# merge bcfpile

for file in $(ls bcfpileup_wubaboon_cpos*)
do
echo $file
grep -v '#' $file >> bcfpileup_wubaboon_merged.vcf
done

#

bcftools annotate -x INFO/CSQ,INFO/InbreedingCoeff,INFO/BaseQRankSum,INFO/ExcessHet,INFO/ClippingRankSum,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQRankSum,INFO/ReadPosRankSum,INFO/SOR \
wu19/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.vcf.gz -R ../cpos_DP_GQ.txt > gatk_candidates_DPGQ.trio19.vcf

bcftools annotate -x INFO/CSQ,INFO/InbreedingCoeff,INFO/BaseQRankSum,INFO/ExcessHet,INFO/ClippingRankSum,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQRankSum,INFO/ReadPosRankSum,INFO/SOR \
wu1/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.vcf.gz -R ../cpos_DP_GQ.txt > gatk_candidates_DPGQ.trio1.vcf

excl_samples=$(tail -n +16 ../trio_tableMFC.csv | tr ',' '\n' | sort | uniq -u | tr '\n' ',')
bcftools view -Oz -s ^${excl_samples::-1} gatk_candidates_DPGQ.trio19.vcf -o gatk_candidates_DPGQ_unphased.trio19.vcf.gz
bcftools view -Oz -s ${excl_samples::-1} gatk_candidates_DPGQ.trio1.vcf -o gatk_candidates_DPGQ_unphased.trio1.vcf.gz

bcftools merge -m snps gatk_candidates_DPGQ_unphased.trio19.vcf.gz gatk_candidates_DPGQ_unphased.trio1.vcf.gz -o gatk_candidates_DPGQ_unphased.vcf
grep -v '##' gatk_candidates_DPGQ_unphased.vcf > gatk_candidates_DPGQ_unphased.decap.vcf