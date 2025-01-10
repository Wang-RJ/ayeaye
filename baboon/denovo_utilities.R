library(dplyr)
options(stringsAsFactors = FALSE)

### Extract functions

extract_genotype_entry <- function(genotypes_list, column) {
  return(lapply(genotypes_list, function(df) {
    return(as.data.frame(lapply(df, function(geno_byindividual) {
      colsplit <- strsplit(geno_byindividual, ":")
      return(sapply(colsplit, "[", column))
    }), check.names = FALSE))
  }))
}

extract_genotypes <- function(vcf_tables) {
  genotypes_list <- lapply(vcf_tables, select, 6:8)
  return(extract_genotype_entry(genotypes_list, 1))
}

extract_depths <- function(vcf_tables) {
  genotypes_list <- lapply(vcf_tables, select, 6:8)
  return(extract_genotype_entry(genotypes_list, 3) %>% lapply(., apply, c(1,2), as.numeric) %>% lapply(., as.data.frame, check.names = FALSE))
}

extract_gq <- function(vcf_tables) {
  genotypes_list <- lapply(vcf_tables, select, 6:8)
  return(extract_genotype_entry(genotypes_list, 4) %>% lapply(., apply, c(1,2), as.numeric) %>% lapply(., as.data.frame, check.names = FALSE))
}

extract_allelic_balance <- function(vcf_tables) {
  genotypes_list <- lapply(vcf_tables, select, 6:8)
  AB_list <- extract_genotype_entry(genotypes_list, 2)
  
  n_allele_ref <- AB_list %>% lapply(., select, 3) %>% lapply(., unlist) %>% lapply(., strsplit, ",") %>% lapply(., sapply, "[", 1) %>% lapply(., as.numeric) # Select 3 for kid
  n_allele_alt <- AB_list %>% lapply(., select, 3) %>% lapply(., unlist) %>% lapply(., strsplit, ",") %>% lapply(., sapply, "[", 2) %>% lapply(., as.numeric)
  n_allele_tot <- mapply("+", n_allele_ref, n_allele_alt)
  
  return(mapply("/", n_allele_alt, n_allele_tot))
}

extract_neighbor_distance <- function(vcf_tables) {
  splitbychr <- mapply(split, vcf_tables, lapply(vcf_tables, "[[", "CHROM"), SIMPLIFY = FALSE)
  diffsbychr <- lapply(splitbychr, function(x) { return(lapply(x, "[[", "POS") %>% lapply(., diff)) })
  
  return(diffsbychr)
}

extract_neighbors <- function(vcf_tables, distance) {
  splitbychr <- mapply(split, vcf_tables, lapply(vcf_tables, "[[", "CHROM"), SIMPLIFY = FALSE)
  diffsbychr <- lapply(splitbychr, function(x) { return(lapply(x, "[[", "POS") %>% lapply(., diff)) })
  
  neighbors <- mapply(function(mvc_split, diff_split) {
    lapply(names(diff_split), function(chr) {
      return(mvc_split[[chr]][sort(unique(c(which(diff_split[[chr]] <= distance), which(diff_split[[chr]] <= distance) + 1))),])
    })
  }, splitbychr, diffsbychr, SIMPLIFY = FALSE)
  
  return(lapply(neighbors, function(x) { do.call(rbind, x) }))
}

##

filter_depth <- function(vcf_tables, lower_depth, upper_depth) {
  depths <- extract_depths(vcf_tables)
  
  depthbool_lower <- lapply(depths, function(depth_frame) { apply(depth_frame, 2, ">", lower_depth) %>% apply(., 1, all) })
  depthbool_upper <- lapply(depths, function(depth_frame) { apply(depth_frame, 2, "<", upper_depth) %>% apply(., 1, all) })
  depthbool <- mapply("&", depthbool_lower, depthbool_upper)
  
  return(mapply(function(candidate_frame, depth_boolean) { return(candidate_frame[depth_boolean,]) },
                vcf_tables, depthbool, SIMPLIFY = FALSE))
}

filter_depth_var <- function(vcf_tables, lower_dp_map, upper_dp_map) {
  depths <- extract_depths(vcf_tables)
  
  depthbool_lower <- lapply(depths, function(depth_frame) {
    mapply(">", depth_frame, lower_dp_map[names(depth_frame)]) %>% apply(., 1, all) })
  depthbool_upper <- lapply(depths, function(depth_frame) {
    mapply("<", depth_frame, upper_dp_map[names(depth_frame)]) %>% apply(., 1, all) })
  depthbool <- mapply("&", depthbool_lower, depthbool_upper)
  
  return(mapply(function(candidate_frame, depth_boolean) { return(candidate_frame[depth_boolean,]) },
                vcf_tables, depthbool, SIMPLIFY = FALSE))
}

filter_gq <- function(vcf_tables, thresh) {
  gqs <- extract_gq(vcf_tables)
  
  gq_bool <- lapply(gqs, function(gq_frame) { apply(gq_frame, 2, ">", thresh) %>% apply(., 1, all) })
  return(mapply(function(candidate_frame, gq_boolean) { return(candidate_frame[gq_boolean,]) },
                vcf_tables, gq_bool, SIMPLIFY = FALSE))
}

filter_allelicbalance <- function(vcf_tables, lower_ab, upper_ab) {
  vcf_tables <- vcf_tables[lapply(vcf_tables, nrow) > 0]
  
  ab_list <- extract_allelic_balance(vcf_tables)
  
  ab_bool <- lapply(ab_list, function(ab_vec) { return(which(ab_vec >= lower_ab & ab_vec <= upper_ab)) } )
  return(mapply(function(vcf_frame, include_vector) { return(vcf_frame[include_vector,]) },
                vcf_tables, ab_bool, SIMPLIFY = FALSE))
}

###

bind_candidates <- function(vcf_tables) {
  vcf_tables <- vcf_tables[lapply(vcf_tables, nrow) > 0]
  
  gq_list <- extract_gq(vcf_tables)
  gq_list <- lapply(gq_list, setNames, c("Mother_GQ", "Father_GQ", "Child_GQ"))
  dp_list <- extract_depths(vcf_tables)
  dp_list <- lapply(dp_list, setNames, c("Mother_DP", "Father_DP", "Child_DP"))
  ab_list <- extract_allelic_balance(vcf_tables)
  transflag_list <- extract_transflag(vcf_tables)
  
  return(data.frame(do.call(rbind, lapply(vcf_tables, select, 1:4)),
                    do.call(rbind, gq_list),
                    do.call(rbind, dp_list),
                    allelic_balance = unlist(ab_list),
                    trans_flag = unlist(transflag_list),
                    trio = unlist(mapply(rep, names(vcf_tables), lapply(vcf_tables, nrow))),
                    row.names = NULL))
}
