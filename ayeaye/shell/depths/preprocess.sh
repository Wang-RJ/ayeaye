# Index bam files
for bamfile in $(ls *.bam)
do
echo $bamfile
sbatch --job-name=idx_$bamfile --export=bamfile=$bamfile bamindex.sbatch
done

# Generate raw depths from samtools
for i in {0..19}
do
sbatch --job-name=dpbatch$i --export=batch=$i ayeayedepth.sbatch
done

# Create trio ordered table for bams, ie. trio_bam_order
#
#!/bin/bash
#ls -l ../bam_files/*.bam | awk '{print $9}' | sed 's/\.\.\/bam_files\///g' > bam_order
#awk '{print $0, NR}' <(sed 's/\.bam//g' bam_order) > bam_order.idx
#awk 'NR==FNR {a[$1]=$2;next}{print a[$1], a[$2], a[$3]}' bam_order.idx ../trio_tableMFC.csv > trio_bam_order

# Calculate mean/sums from batches 0-19
sbatch --job-name=dp_batchmeans mean_sums.sbatch

# Create tables, including for max depths by trio
for file in $(ls dpcounts.batch* | sort -V)
do
head -n 1 $file >> sums_by_batch.txt
tail -n +2 $file | head -n 1 >> sqsums_by_batch.txt
tail -n 1 $file >> totpos_by_batch.txt
done

awk '{for(i=1;i<=NF;i++) {sum[i] += $i}} END {for(i=1;i<=NF;i++) print sum[i]}' sums_by_batch.txt > sums_by_id.txt
paste -d/ sums_by_id.txt <(yes $(paste -sd+ totpos_by_batch.txt | bc) | head -n 18) | bc -l > mean_depth_by_id.txt
awk '{print $1 + 4 *sqrt($1)}' mean_depth_by_id.txt > max_depth_by_id.txt

paste tbocounts* | tr ' ' '\t' | awk '{for(i=2;i<=NF;i+=2) { sum[NR] += $i } print sum[NR]}' > callable_sites.txt

# Send table creation to slurm
sbatch --job-name=tbo_counts tbo_ayeaye.sbatch