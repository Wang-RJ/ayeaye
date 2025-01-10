# Must be in ayeaye/
# Process bam files to remove likely GATK remapping errors

rawbam_bcfcalls <- read.table("vcf_files/bcfpileup_ayeaye.decap.vcf", comment.char = "", header = FALSE, check.names = FALSE)
rawbam_header <- read.table("vcf_files/bcfpileup_ayeaye.decap.header", comment.char = "", check.names = FALSE)

names(rawbam_bcfcalls) <- rawbam_header
names(rawbam_bcfcalls)[1] <- "CHROM"

rb_mvcf <- rawbam_bcfcalls[,c(1,2,10:ncol(rawbam_bcfcalls))]
rb_mvcf_gt <- rb_mvcf[,3:ncol(rb_mvcf)]

rb_mvcf <- rawbam_bcfcalls[,c(1,2,10:ncol(rawbam_bcfcalls))]
rb_mvcf_gt <- as.matrix(rb_mvcf[,3:ncol(rb_mvcf)])
rb_mvcf_AD <- as.data.frame(apply(rb_mvcf_gt, 2, function(col) { return(strsplit(col, ":") %>% sapply(., function(entry) { return(entry[length(entry)]) })) }))
rb_mvcf_DP <- as.data.frame(apply(rb_mvcf_gt, 2, function(col) { return(strsplit(col, ":") %>% sapply(., function(entry) { return(as.numeric(entry[length(entry)-3])) })) }))
rb_mvcf_ADF <- as.data.frame(apply(rb_mvcf_gt, 2, function(col) { return(strsplit(col, ":") %>% sapply(., function(entry) { return(entry[length(entry)-2]) })) }))
rb_mvcf_ADR <- as.data.frame(apply(rb_mvcf_gt, 2, function(col) { return(strsplit(col, ":") %>% sapply(., function(entry) { return(entry[length(entry)-1]) })) }))

rb_mvcf_homoref <- as.data.frame(apply(rb_mvcf_AD, 2, function(col) {
  splitcol <- sapply(col, strsplit, ",")
  single_allele <- sapply(splitcol, length) < 2                      # No ALT entry, e.g. AD = 23
  homo_allele <- as.numeric(sapply(splitcol, tail, 1)) == 0          # No read is alternate, e.g. AD = 23,0
  
  return(single_allele | homo_allele)
}))

rb_mvcf_homorefAD1 <- as.data.frame(apply(rb_mvcf_AD, 2, function(col) {
  splitcol <- sapply(col, strsplit, ",")
  single_allele <- sapply(splitcol, length) < 2                      # No ALT entry, e.g. AD = 23
  homo_allele <- as.numeric(sapply(splitcol, tail, 1)) < 2           # No more than 1 read is alternate, e.g. AD = 22,1
  
  return(single_allele | homo_allele)
}))

rb_mvcf_hetF <- as.data.frame(apply(rb_mvcf_ADF, 2, function(col) {
  splitcol <- sapply(col, strsplit, ",")
  single_allele <- sapply(splitcol, length) < 2                     # Must have an ALT entry
  het_allele <- as.numeric(sapply(splitcol, tail, 1)) > 0           # ALT entry must be > 0
  
  return(!single_allele & het_allele)
}))

rb_mvcf_hetR <- as.data.frame(apply(rb_mvcf_ADR, 2, function(col) {
  splitcol <- sapply(col, strsplit, ",")
  single_allele <- sapply(splitcol, length) < 2
  het_allele <- as.numeric(sapply(splitcol, tail, 1)) > 0
  
  return(!single_allele & het_allele)
}))

rb_mvcf_homoref_df <- rb_mvcf
rb_mvcf_homoref_df[,3:ncol(rb_mvcf_homoref_df)] <- rb_mvcf_homoref
rb_mvcf_homoref_df$upos <- paste(rb_mvcf_homoref_df$CHROM, rb_mvcf_homoref_df$POS)

rb_mvcf_homorefAD1_df <- rb_mvcf
rb_mvcf_homorefAD1_df[,3:ncol(rb_mvcf_homorefAD1_df)] <- rb_mvcf_homorefAD1
rb_mvcf_homorefAD1_df$upos <- paste(rb_mvcf_homorefAD1_df$CHROM, rb_mvcf_homorefAD1_df$POS)

rb_mvcf_hetF_df <- rb_mvcf
rb_mvcf_hetF_df[,3:ncol(rb_mvcf_hetF_df)] <- rb_mvcf_hetF
rb_mvcf_hetF_df$upos <- paste(rb_mvcf_hetF_df$CHROM, rb_mvcf_hetF_df$POS)

rb_mvcf_hetR_df <- rb_mvcf
rb_mvcf_hetR_df[,3:ncol(rb_mvcf_hetR_df)] <- rb_mvcf_hetR
rb_mvcf_hetR_df$upos <- paste(rb_mvcf_hetR_df$CHROM, rb_mvcf_hetR_df$POS)

rb_mvcf_DP_df <- rb_mvcf
rb_mvcf_DP_df[,3:ncol(rb_mvcf_DP_df)] <- rb_mvcf_DP
rb_mvcf_DP_df$upos <- paste(rb_mvcf_DP_df$CHROM, rb_mvcf_DP_df$POS)

momlookup <- as.character(trio_table$Mother)
names(momlookup) <- trio_table$trio
dadlookup <- as.character(trio_table$Father)
names(dadlookup) <- trio_table$trio
kidlookup <- as.character(trio_table$Child)
names(kidlookup) <- trio_table$trio

get_bamparents_homobool <- function(mv_frame, trio) {
  rb_bool_table <- rb_mvcf_homoref_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_homoref_df$upos),]

  return(rb_bool_table[,momlookup[trio]] & rb_bool_table[,dadlookup[trio]])
}

get_bamparents_homoboolAD1 <- function(mv_frame, trio) {
  rb_bool_table <- rb_mvcf_homorefAD1_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_homorefAD1_df$upos),]

  return(rb_bool_table[,momlookup[trio]] & rb_bool_table[,dadlookup[trio]])
}

get_bamchild_hetbool <- function(mv_frame, trio) {
  rb_bool_tableF <- rb_mvcf_hetF_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_hetF_df$upos),]
  rb_bool_tableR <- rb_mvcf_hetR_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_hetR_df$upos),]

  return(rb_bool_tableF[,kidlookup[trio]] & rb_bool_tableR[,kidlookup[trio]])
}

filter_bam <- function(vcf_tables) {
  parentbools <- mapply(get_bamparents_homobool, vcf_tables, names(vcf_tables))
  childbools <- mapply(get_bamchild_hetbool, vcf_tables, names(vcf_tables))
  
  return(mapply(function(candidate_frame, pbool, cbool) { return(candidate_frame[(pbool & cbool),]) },
                vcf_tables, parentbools, childbools, SIMPLIFY = FALSE))
}

filter_bam_dpvar <- function(vcf_tables, lower_depth, upper_depth) {
  get_DPbool <- function(mv_frame, trio) {
    dp_table <- rb_mvcf_DP_df[match(paste(mv_frame$CHROM, mv_frame$POS), rb_mvcf_DP_df$upos),]
    
    mom_id <- momlookup[trio]
    dad_id <- dadlookup[trio]
    kid_id <- kidlookup[trio]
    
    mom_bool <- (dp_table[,mom_id] > lower_depth[mom_id]) & dp_table[,mom_id] < (upper_depth[mom_id])
    dad_bool <- (dp_table[,dad_id] > lower_depth[dad_id]) & dp_table[,dad_id] < (upper_depth[dad_id])
    kid_bool <- (dp_table[,kid_id] > lower_depth[kid_id]) & dp_table[,kid_id] < (upper_depth[kid_id])
    
    return(mom_bool & dad_bool & kid_bool)
  }
  
  dpbools <- mapply(get_DPbool, vcf_tables, names(vcf_tables))
  return(mapply(function(candidate_frame, dpbool) { return(candidate_frame[dpbool,]) }, vcf_tables,
         dpbools, SIMPLIFY = FALSE))
}

filter_bamAD1 <- function(vcf_tables) {
  parentbools <- mapply(get_bamparents_homoboolAD1, vcf_tables, names(vcf_tables))
  childbools <- mapply(get_bamchild_hetbool, vcf_tables, names(vcf_tables))
  
  return(mapply(function(candidate_frame, pbool, cbool) { return(candidate_frame[(pbool & cbool),]) },
                vcf_tables, parentbools, childbools, SIMPLIFY = FALSE))
}
