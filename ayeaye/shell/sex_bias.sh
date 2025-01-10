#!/bin/bash

# Phylogenetic Sex Bias in Strepsirrhines

set -e
set -o pipefail
set -u


ALIGNMENT_URL="https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.hal"
ALIGNMENT_FILE="241-mammalian-2020v2.hal"
SPECIES_FILE="species.txt"
SPECIES_LIST=$(cat "$SPECIES_FILE")
REF_GENOME="Microcebus_murinus"
AUTOSOMES=$(seq 1 32)
WORK_DIR=$(pwd)
MAF_DIR="$WORK_DIR/maf_chr"
FILTERED_DIR="$WORK_DIR/filtered"
LIKELIHOODS_DIR="$WORK_DIR/likelihoods"
TREES_FILE="$WORK_DIR/trees.txt"


# Download 241-way mammalian alignment
# Date: September 13, 2023
download_alignment() {
  if [[ ! -f $ALIGNMENT_FILE ]]; then
    echo "Downloading alignment..."
    wget -q "$ALIGNMENT_URL" -O "$ALIGNMENT_FILE"
  else
    echo "Alignment file already exists."
  fi
}


# Convert alignment to MAF format for each chromosome
convert_to_maf() {
	echo "Converting alignments to MAF format..."
	mkdir -p maf_chr
	for i in {61..92}; do
	    hal2maf $ALIGNMENT_FILE maf_chr/chr"$i".maf \
	      --targetGenomes $SPECIES_LIST \
	      --refGenome $REF_GENOME \
	      --refSequence CM0076"$i".1 \
	      --noDupes --noAncestors --onlyOrthologs &
	  done
  wait
}


# Determine the quality of the assembly of each strepsirrhines species
## stats.bpp
SPECIES=(hg38,Mirza_coquereli,Microcebus_murinus,Cheirogaleus_medius,Indri_indri,Propithecus_coquereli,Eulemur_flavifrons,Eulemur_fulvus,Lemur_catta,Daubentonia_madagascariensis,Otolemur_garnetti,Nycticebus_coucang)
input.file=$(DATA).maf.gz
input.file.compression=gzip
input.format=Maf
input.dots=as_gaps
maf.filter=                                 \
    SequenceStatistics(                     \
        statistics=(\
            BlockLength,                    \
            SequenceLength(                 \
                species=$(SP))),                \
                ref_species=hg38,                   \
        file=$(DATA).$(SP).statistics.csv.gz,       \
                compression=gzip)

## Stasts function
run_maffilter_stats() {
  echo "Running maffilter statistics..."
  for CHR in "${AUTOSOMES[@]}"; do
    for SPECIES in $(cat "$SPECIES_FILE" | tr ',' ' '); do
      maffilter param=stats.bpp DATA=chr"$CHR" SP="$SPECIES" &
    done
    wait
  done
}
# Select species based on quality, if needed


# Filter alignments with species of interest
## filter_sp.bpp,
SPECIES=(Mirza_coquereli,Microcebus_murinus,Cheirogaleus_medius,Indri_indri,Propithecus_coquereli,Eulemur_flavifrons,Lemur_catta,Daubentonia_madagascariensis,Otolemur_garnettii,Nycticebus_coucang,Homo_sapiens)

input.file=$(DATA).maf
input.file.compression=none
input.format=Maf
input.dots=as_gaps
output.log=species/$(DATA).maffilter.log
maf.filter=                                 \
    Subset(species=$(SPECIES), strict=yes, keep=no, remove_duplicates=yes, verbose=no),\
    Output(file=species/$(DATA).maf.gz, compression=gzip)

## function to extract species of interest
filter_species() {
  echo "Filtering alignments for species of interest..."
  mkdir -p "$FILTERED_DIR"
  for CHR in "${AUTOSOMES[@]}"; do
    maffilter param=filter_sp.bpp DATA=chr"$CHR"
    rm -f species/chr"$CHR".maffilter.log
  done
  maffilter param=filter_sp.bpp DATA=chrX
  rm -f species/chrX.maffilter.log
}


# Filter maf file further processing
for CHR in "${AUTOSOMES[@]}"; do
  gzip -d species/chr$CHR.maf.gz
done

gzip -d species/chrX.maf.gz

# function for exon and PAR removal:
remove_exons_and_PAR() {
  echo "Filtering out exons from autosomes..."
  for CHR in "${AUTOSOMES[@]}"; do
    maf_parse --features features/chr"$CHR".gtf -M "$(cat $SPECIES_FILE)" \
      "$FILTERED_DIR/chr${CHR}.maf" > "$FILTERED_DIR/chr${CHR}_noexon.maf"
  done
  echo "Filtering out exons and PAR from chrX..."
  maf_parse --features features/chrX.gtf -M "$(cat $SPECIES_FILE)" \
      "$FILTERED_DIR/chrX.maf" > "$FILTERED_DIR/chrX_noexon.maf"
  maf_parse --features features/par_coord.bed -M "$(cat $SPECIES_FILE)" \
    "$FILTERED_DIR/chrX_noexon.maf" > "$FILTERED_DIR/chrX_nopar.maf"
}


# Phylogenetic analysis with phyloFit
## the Species tree, sp_tree.txt, looks like this:

# (((Nycticebus_coucang,Otolemur_garnettii),(Daubentonia_madagascariensis,(((Propithecus_coquereli,Indri_indri),(Cheirogaleus_medius,(Microcebus_murinus,Mirza_coquereli))),(Lemur_catta,Eulemur_flavifrons)))),Homo_sapiens);

# Note: phyloFit was run 6 times, and the output with the highest likelihood was chosen

run_phylofit() {
  echo "Running phyloFit..."
  mkdir -p "$LIKELIHOODS_DIR"

  # Loop for autosomes
  for CHR in "${AUTOSOMES[@]}"; do
    for i in {1..6}; do
      phyloFit -r --EM --precision MED --subst-mod UNREST -Z \
        --msa-format MAF "$FILTERED_DIR/chr${CHR}_noexon.maf" \
        --tree sp_tree.txt \
        -e "$LIKELIHOODS_DIR/errors_chr${CHR}.${i}" \
        -o "$LIKELIHOODS_DIR/chr${CHR}_output.${i}" &
    done
  done
  wait

  # Loop for chrX
  for j in {1..6}; do
    phyloFit -r --EM --precision MED --subst-mod UNREST -Z \
      --msa-format MAF "$FILTERED_DIR/chrX_nopar.maf" \
      --tree sp_tree.txt \
      -e "$LIKELIHOODS_DIR/errors_chrX.${j}" \
      -o "$LIKELIHOODS_DIR/chrX_output.${j}"
  done
}


# Main functions
download_alignment
convert_to_maf
run_maffilter_stats
filter_species
remove_exons_and_PAR
run_phylofit


# Likelihood extraction and processing
echo "Extracting likelihood data..."
mkdir likelihoods
for CHR in "${AUTOSOMES[@]}"; do
  for j in {1..6}; do
    head -n 4 chr"$CHR"_output."$j".mod | tail -n 1 | cut -c 15- >> likelihoods/chr$CHR.lnl.txt
  done
done

for j in {1..6}; do
  head -n 4 chrX_output."$j".mod | tail -n 1 | cut -c 15- >> likelihoods/chrX.lnl.txt
done


# Select best likelihood reps
echo "Selecting the best likelihood rep for each chromosome..."
for CHR in "${AUTOSOMES[@]}"; do
  sort -gk2,2 <(paste rank.txt chr$CHR.lnl.txt) | tail -n 1 >> highest_lklhd.txt
done

sort -gk2,2 <(paste rank.txt chrX.lnl.txt) | tail -n 1 >> highest_lklhd.txt

# Save trees with the highest likelihoods
echo "Saving trees with best likelihoods..."
for i in {1..33}; do
  chr=$(head -n "$i" highest_lklhd_chr.txt | tail -n 1 | awk '{print $1}')
  rep=$(head -n "$i" highest_lklhd_chr.txt | tail -n 1 | awk '{print $2}')
  tail -n 1 ../"$chr"_output."$rep".mod >> trees.txt
done

paste chr.txt trees.txt > trees_chr.txt

# Calculate branch lengths across autosomes
echo "Calculating average branch lengths..."
head -n -1 trees.txt | awk '{print $2}' | awk '{ gsub(/:|,|)|:/, " ") } 1' | sed 's/[^0-9 .]*//g' > branch_lengths_chr.txt