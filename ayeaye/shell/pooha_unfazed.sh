#!/bin/bash

# Define paths to files and directories
TRIO_FILE="/N/project//ayeaye/trio_tableMFC.csv"
VCF_FILE="/N/project//ayeaye/vcf_files/1/output.filtered.snps.removed.CalculateGenotypePosteriors.deNovo.only.gz"
BAM_DIR="/N/project//ayeaye/bam_files/"
POOHA_OUT_DIR="./pooha/"
UNFAZED_OUT_DIR="./unfazed/"

# Generate a BED file with regions spanning ±1000 bp around each mutation
awk '{start=$2-1000; end=$2+1000; print $1"\t"start"\t"end}' all_mutations.txt > mutations_1kb.bed

# Generate a VCF with regions spanning ±1000 bp around each mutation
bcftools view -R mutations_1kb.bed "$VCF_FILE" -Oz -o 1kb.vcf.gz

# index
bcftools index 1kb.vcf.gz

# create output directories for each program
mkdir -p "$POOHA_OUT_DIR" "$UNFAZED_OUT_DIR"

# run POOHA
conda activate pooha
for i in {1..12}; do
    trio=$(awk -F ',' "NR==$i {print \$2, \$1, \$3}" "$TRIO_FILE")
    bam=$(awk -F ',' "NR==$i {print \$3}" "$TRIO_FILE")
    POOHA 1kb.vcf $trio "$BAM_DIR/$bam.bam" >> "$POOHA_OUT_DIR/phased_mut.txt"
done

# Count phased mutations in POOHA
grep -f <(awk '{print $2}' all_mutations.txt) "$POOHA_OUT_DIR/phased_mut.txt" | grep 'phased' | wc -l

# Prepare mutation BED files for each trio for Unfazed
for i in {1..12}; do
    grep -w "trio$i" all_mutations.txt | awk '{print $1, $2}' OFS='\t' > "mut_trio$i.txt"
    awk '{print $1, $2-1, $2}' OFS="\t" "mut_trio$i.txt" > "mut_trio$i.bed"
done

for i in {1..12}; do
    awk -F ',' -v line=$i -v mut_file="mut_trio$i.bed" 'NR==line {val=$3}
        END {while ((getline line < mut_file) > 0) print line "\t" val "\tPOINT"}' "$TRIO_FILE" > "dnm_trio$i.bed"
done

# Run Unfazed
for i in {1..12}; do
    unfazed -d dnm_trio"$i".bed \
        -s 1kb.vcf.gz \
        -p /N/project//ayeaye/vcf_files/"$i"/trio.ped \
        -b "$BAM_DIR" \
        -g na \
        -o bed > "$UNFAZED_OUT_DIR/phased_trio$i.txt"
done

# Combine Unfazed results
for i in {1..12}; do
    grep -v "#" "$UNFAZED_OUT_DIR/phased_trio$i.txt" >> "$UNFAZED_OUT_DIR/phased_mut.txt"
done