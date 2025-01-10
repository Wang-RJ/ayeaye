library(dplyr)
options(stringsAsFactors = FALSE)

source("denovo_utilities.R")
## Imports set of functions for dealing with lists of dataframes, each a collection of de novo candidates from a different trio
## Querying functions
#   extract_genotype_entry <- function(genotypes_list, column)
#   extract_genotypes <- function(vcf_tables)
#   extract_depths <- function(vcf_tables)
#   extract_gq <- function(vcf_tables)
#   extract_allelic_balance <- function(vcf_tables)
#   extract_neighbor_distance <- function(vcf_tables) 
#   extract_neighbors <- function(vcf_tables, distance) 
#
## Filtering functions
#   filter_depth <- function(vcf_tables, lower_depth, upper_depth)
#   filter_gq <- function(vcf_tables, thresh)

trio_table <- read.csv("../trio_tableMFC.csv", header = FALSE)
names(trio_table) <- c("Mother", "Father", "Child")
trio_table$trio <- paste("trio", 1:12, sep = "")

vcf_files <- list.files("mvf_files/")[grepl("MVF_trio([1-9]|1[0-9]).decap.vcf", list.files("mvf_files/"))]

candidate_table <- list()
for (i in 1:length(vcf_files)) {
  candidate_table[[i]] <- read.table(paste("mvf_files/", vcf_files[i], sep = ""), header = TRUE, comment.char = "", check.names = FALSE)
  candidate_table[[i]]$INFO <- NULL
  candidate_table[[i]]$ID <- NULL
  candidate_table[[i]]$FILTER <- NULL
  candidate_table[[i]]$FORMAT <- NULL
  colnames(candidate_table[[i]])[1] <- "CHROM"
  names(candidate_table)[[i]] <- gsub(".decap.vcf", "", vcf_files[[i]]) %>% gsub("MVF_", "", .)
}

# raw_genocalls <- extract_genotypes(candidate_table)
# momcalls <- raw_genocalls %>% lapply(., select, 1) %>% lapply(., table)
# momcalls       # All homozygous ref
# dadcalls <- raw_genocalls %>% lapply(., select, 2) %>% lapply(., table)
# dadcalls       # All homozygous ref
# kidcalls <- raw_genocalls %>% lapply(., select, 3) %>% lapply(., table)
# kidcalls       # All hets

# raw_allelicbalance <- extract_allelic_balance(candidate_table)
# raw_depths <- extract_depths(candidate_table)
# raw_gq <- extract_gq(candidate_table)

# par(mfrow = c(3,4)); for(i in 1:12) { hist(as.numeric(unlist(raw_gq[[i]])), main = paste("trio", i), xlab = "GQ score") }
# par(mfrow = c(3,4)); for(i in 1:12) { hist(as.numeric(unlist(raw_depths[[i]])), main = paste("trio", i), xlab = "depth") }
# par(mfrow = c(3,4)); for(i in 1:12) { hist(as.numeric(unlist(raw_depths[[i]]) %>% .[. < 80]), main = paste("trio", i), xlab = "depth < 80") }
# par(mfrow = c(3,4)); for(i in 1:12) { hist(as.numeric(raw_allelicbalance[[i]]), breaks = 30, main = paste("trio", i), xlab = "allelic balance (child)") }

# unlist(lapply(candidate_table, nrow))

maxdepths <- read.csv("maxdepths.csv", header = FALSE, check.names = FALSE)
upper_dp_bounds <- setNames(maxdepths$V2, maxdepths$V1)
lower_dp_bounds <- setNames(rep(15, nrow(maxdepths)), maxdepths$V1)

candidates_DP <- filter_depth_var(candidate_table, lower_dp_bounds, upper_dp_bounds)
candidates_DP_GQ <- filter_gq(candidates_DP, 20)

# dpgq_allelicbalance <- extract_allelic_balance(candidates_DP_GQ)
# dpgq_depths <- extract_depths(candidates_DP_GQ)
# dpgq_gq <- extract_gq(candidates_DP_GQ)

# par(mfrow = c(3,4)); for(i in 1:12) { hist(as.numeric(unlist(dpgq_gq[[i]])), main = names(candidates_DP_GQ)[[i]], xlab = "GQ score") }
# par(mfrow = c(3,4)); for(i in 1:12) { hist(as.numeric(unlist(dpgq_depths[[i]]) %>% .[. < 80]), breaks = 20, main = names(candidates_DP_GQ)[[i]], xlab = "depth < 80") }
# par(mfrow = c(3,4)); for(i in 1:12) { hist(as.numeric(dpgq_allelicbalance[[i]]), breaks = 15, xlim = c(0, 1), main = names(candidates_DP_GQ)[[i]], xlab = "allelic balance (child)") }

# par(mfrow = c(1,1)); hist(unlist(dpgq_allelicbalance), breaks = 40, xlab = "allelic balance", main = "allelic balance (child)")
# par(mfrow = c(1,1)); hist(unlist(dpgq_depths), xlab = "depth", main = "filtered [20,60]")
# par(mfrow = c(1,1)); hist(unlist(dpgq_gq), xlab = "GQ", main = "filtered > 20")

# unlist(lapply(candidates_DP_GQ, nrow))

# candidate_positions <- lapply(candidates_DP_GQ, select, 1:2) %>% do.call(rbind, .)
# candidate_positions <- candidate_positions[order(candidate_positions[,1], candidate_positions[,2]),]
# candidate_positions <- unique(candidate_positions)

# write.table(candidate_positions, file = "cpos_DP_GQ.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

source("processBAM.R")
# Imports functions:
# get_bamparents_homobool <- function(mv_frame, trio)   # Takes trio name, e.g. "trio1"
# get_bamchild_hetbool <- function(mv_frame, trio)
# filter_BAM <- function(vcf_tables)

source("processRelated.R")
# Imports functions:
# filter_rel <- function(vcf_tables)
# extract_transflag <- function(vcf_tables)

# candidates_DP_GQ_BAM <- filter_bam(candidates_DP_GQ)
# candidates_DP_GQ_BAM_REL <- filter_rel(candidates_DP_GQ_BAM)

# Quick 2x2 plot of candidates with different allelic balance cutoffs
plot_GQ <- function(AB_cutoff) {
  par(mfrow = c(2,2))
  for(GQ in seq(20,80,20)) {
    tmp <- filter_gq(candidates_DP, GQ) %>% filter_bam %>% filter_rel
    tmp_ab <- candidates_beforeAB$allelic_balance
    hist(tmp_ab, xlab = "allelic balance", ylim = c(0,200), xlim = c(0,.9), breaks = 20,
         main = '',
         xaxt = 'n')
    hist(tmp_ab %>% .[.>AB_cutoff], breaks = 10, col = 'pink', add = TRUE)
    hist(tmp_hifemale$allelic_balance, breaks = 10, col = 'red', add = TRUE)
    axis(1, at = seq(0, 0.9, by = 0.1))
      }
  par(mfrow = c(1,1))
}
# plot_GQ(0.35)

candidates_DP_GQ60 <- filter_gq(candidates_DP, 60)
candidates_DP_GQ60_BAM <- filter_bam(candidates_DP_GQ60) %>% filter_bam_dpvar(., lower_dp_bounds, upper_dp_bounds)
candidates_DP_GQ60_BAM_REL <- candidates_DP_GQ60_BAM %>% filter_rel
candidates_beforeAB <- bind_candidates(candidates_DP_GQ60_BAM_REL)

candidates_DP_GQ60_BAM_REL_AB <- candidates_DP_GQ60_BAM_REL %>% filter_allelicbalance(., 0.35, 1)
candidates_afterAB <- filter_allelicbalance(candidates_DP_GQ60_BAM_REL, 0.35, 1) %>% bind_candidates
denovo_candidates <- candidates_afterAB
