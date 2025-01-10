#!/bin/bash
# ls -l ../bam_files/*.bam | awk '{print $9}' | sed 's/\.\.\/bam_files\///g' > bam_order
# awk '{print $0, NR}' <(sed 's/\.bam//g' bam_order) > bam_order.idx
# awk 'NR==FNR {a[$1]=$2;next}{print a[$1], a[$2], a[$3]}' bam_order.idx ../trio_tableMFC.tsv > trio_bam_order

dp_file=$1
awk -v tbo_file="trio_bam_order" '
  BEGIN {
    n = 1;
    while((getline < tbo_file) > 0) {
      split($0, line);
      tbo[n,1] = line[1];
      tbo[n,2] = line[2];
      tbo[n,3] = line[3];
      n++;
    }
    close(listfile)
  }
  {
    for(trio = 1; trio <= 6; trio++) {
      mom_idx = tbo[trio,1] + 2
      dad_idx = tbo[trio,2] + 2
      kid_idx = tbo[trio,3] + 2
      if($mom_idx > 20 && $mom_idx < 60 &&
         $dad_idx > 20 && $dad_idx < 60 &&
         $kid_idx > 20 && $kid_idx < 60) count[trio]++
    }
  }
  END {
    for(trio = 1; trio <= 6; trio++) print trio,count[trio]
  }
' $dp_file
