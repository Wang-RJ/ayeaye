#!/bin/bash
cut --complement -f1-2 $1 | awk '{for(i=1;i<=NF;i++) { sum[i] += $i; sumsq[i] += ($i)^2 }}
END {for (i=1;i<=NF;i++) { printf("%d ", sum[i]) } \
print "";
for (i=1;i<=NF;i++) { printf("%d ", sumsq[i]) }
print "";
print NR}'
