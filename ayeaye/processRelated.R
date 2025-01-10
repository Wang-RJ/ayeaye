# Compares genotypes of related individuals

# fix some of the sample names
gatkgt_vcf <- read.table("vcf_files/gatk_candidates_DPGQ_unphased.decap.vcf", header = TRUE, comment.char = "", check.names = FALSE)
names(gatkgt_vcf)[1] <- "CHROM"

# reduce to just the genotype
gatkgt_brief <- gatkgt_vcf[,c(1,2,4,5,10:length(gatkgt_vcf))]
gatkgt_brief[,5:ncol(gatkgt_brief)] <- sapply(gatkgt_brief[,5:ncol(gatkgt_brief)], function(x) { return(substr(x, 1, 3)) })
gatkgt_brief[,5:ncol(gatkgt_brief)] <- apply(gatkgt_brief[,5:ncol(gatkgt_brief)], MARGIN = c(1,2), function(x) { gsub("\\|", "/", x) })

# downstream_descendants pulls all related descendants of that trio (includes descendants of half-sibs to children of the trio, e.g. potential germline mosaicism)
library(memoise)

downstream_descendants <- function(trio, trio_table) {
  trio_idx <- match(trio, trio_table$trio)
  trio_entry <- trio_table[trio_idx,]
  
  if(is.na(trio_entry)[4])
    return(NULL)
  
  else {
    descendant_trios <- unique(c(which(trio_entry$Father == trio_table$Father),
                                 which(trio_entry$Mother == trio_table$Mother),
                                 which(trio_entry$Child == trio_table$Mother),
                                 which(trio_entry$Child == trio_table$Father)))
    return(unique(
      c(as.character(trio_entry[3]),
        unlist(sapply(trio_table$trio[descendant_trios], memo_dd, trio_table[-trio_idx,])))))
  }
}
memo_dd <- memoise(downstream_descendants)

# make a list of downstream descendants and parents, which should both be excluded for each trio
dd_table <- sapply(trio_table$trio, memo_dd, trio_table)
parents_table <- split(trio_table[,1:2], trio_table$trio)
exclude_table <- sapply(trio_table$trio, function(trio) { return(c(dd_table[[trio]], as.character(parents_table[[trio]]))) })

# Cross-reference entries in mv_autosomes to make sure candidates do not appear in unrelated individuals
# exclude parents and downstream descendants from unrelated individuals

# create a unique position ID to match up mvcf and mv_autosome dataframes
gatkgt_brief_upos <- paste(gatkgt_brief$CHROM, gatkgt_brief$POS)

# is there a nonref allele in an unrelated individual?
crossref_mv_exclude <- function(upos, trio) {
  entry <- gatkgt_brief[match(upos, gatkgt_brief_upos),][,-(1:4)]
  entry <- entry[!(names(entry) %in% exclude_table[[trio]])]
  
  return("0/1" %in% entry | "1/1" %in% entry)
}

filter_rel <- function(vcf_tables) {
  rel_exclude_bools <- mapply(function(vcf_frame, trio) { return(sapply(vcf_frame, crossref_mv_exclude, trio)) },
                              lapply(vcf_tables, function(candidate_frame) { paste(candidate_frame$CHROM, candidate_frame$POS) }),
                              names(vcf_tables), SIMPLIFY = FALSE)
  
  return(mapply(function(candidate_frame, excl_bool) { return(candidate_frame[!excl_bool,]) },
                vcf_tables, rel_exclude_bools, SIMPLIFY = FALSE))
}

# Was an allele at a given position transmitted to at least one grandchild
transmission_check <- function(upos, trio) {
  trio_child <- trio_table$Child
  names(trio_child) <- trio_table$trio
  
  grandchildren <- trio_child[which(trio_child[trio] == trio_table$Father | trio_child[trio] == trio_table$Mother)]
  
  if(length(grandchildren) < 1)
    return(NA)
  
  entry <- gatkgt_brief[match(upos, gatkgt_brief_upos),][,as.character(grandchildren)]

  return("0/1" %in% entry | "1/1" %in% entry)  
}

extract_transflag <- function(vcf_tables) {
  mapply(function(candidate_frame, trio) {
    return(mapply(transmission_check, paste(candidate_frame$CHROM, candidate_frame$POS), trio))
  },
  vcf_tables, names(vcf_tables))
}
