options(stringsAsFactors = FALSE)

library(kinship2)
ped_input_raw <- read.table("../ayeaye.ped")
ped_input <- ped_input_raw[order(ped_input_raw$V3, decreasing = TRUE),]
ped_input <- ped_input[!duplicated(ped_input$V2),]

# ped_input <- ped_input[-c(1,2,3),]
ayeaye_ped <- pedigree(id = ped_input$V2, dadid = ped_input$V3, momid = ped_input$V4, sex = ped_input$V5)
plot(ayeaye_ped, branch = 1, width = 10, pconnect = 1)

randomSNPs <- read.table("vcf_files/ayeaye100k_shufSNPs.vcf", header = TRUE, comment.char = "", check.names = FALSE)

# Strip to genotypes and remove phase
randomSNPs[,10:ncol(randomSNPs)] <- as.data.frame(lapply(randomSNPs[,10:ncol(randomSNPs)], function(gbi) { colsplit <- strsplit(gbi, ":"); return(sapply(colsplit, "[", 1)) }))
randomSNPs[,10:ncol(randomSNPs)] <- apply(randomSNPs[,10:ncol(randomSNPs)], MARGIN = c(1,2), function(x) { gsub("\\|", "/", x) })

calc_CoR <- function(genotypes, ind_i, ind_j) {
  ind_i <- as.character(ind_i)
  ind_j <- as.character(ind_j)
  het_ij <- sum(genotypes[[ind_i]] == "0/1" & genotypes[[ind_j]] == "0/1")
  het_i <- sum(genotypes[[ind_i]] == "0/1")
  het_j <- sum(genotypes[[ind_j]] == "0/1")
  
  ibs0 <- sum((genotypes[[ind_i]] == "0/0" & genotypes[[ind_j]] == "1/1") |
                (genotypes[[ind_i]] == "1/1" & genotypes[[ind_j]] == "0/0"))
  
  ibs2 <- sum(genotypes[[ind_i]] == genotypes[[ind_j]])
  
  return(c((het_ij - 2 * ibs0) / (sqrt(het_i * het_j)), ibs0, ibs2))
}

all_ayeayes <- as.numeric(names(randomSNPs[,10:ncol(randomSNPs)]))

relation_frame <- data.frame(id1 = c(), id2 = c(), coefficient = c(), ibs0 = c())
for(i in 1:(length(all_ayeayes) - 1)) {
  for(j in (i+1):length(all_ayeayes)) {
    ret <- calc_CoR(randomSNPs, all_ayeayes[i], all_ayeayes[j])
    relation_frame <- rbind(relation_frame,
                            data.frame(id1 = as.character(all_ayeayes[i]), id2 = as.character(all_ayeayes[j]), coefficient = ret[1], ibs0 = ret[2], ibs2 = ret[3])
    )
  }
}

relation_frame$ped_relation <- NA

ped_trios <- read.csv("../trio_tableMFC.csv", header = FALSE)
parent_child_pairs <- rbind(data.frame(id1 = ped_trios$V1, id2 = ped_trios$V3), data.frame(id1 = ped_trios$V2, id2 = ped_trios$V3))

sibling_pairs <- data.frame(id1 = c(100933, 100933, 100933, 100944, 100944, 100945, 100950, 100941, 100942),
                            id2 = c(100945, 100944, 100935, 100945, 100935, 100935, 100940, 100946, 100939))
halfsib_pairs <- data.frame(id1 = c(100950, 100950, 100950, 100940, 100940, 100940, 100938, 100938),
                            id2 = c(100941, 100946, 100938, 100941, 100946, 100938, 100941, 100946))
grandp_pairs <- data.frame(id1 = c(100934, 100934, 100934, 100934, 100934, 100949, 100949, 100949, 100949, 100949, 100936, 100936, 100943, 100943, 100937, 100937, 100948, 100948),
                           id2 = c(100938, 100950, 100940, 100941, 100946, 100938, 100950, 100940, 100941, 100946, 100941, 100946, 100941, 100946, 100942, 100939, 100942, 100939))

zip <- paste(relation_frame[,1], relation_frame[,2])

relation_frame$ped_relation[match(paste(parent_child_pairs[,2], parent_child_pairs[,1]), zip)] <- "parent-child"
relation_frame$ped_relation[match(paste(parent_child_pairs[,1], parent_child_pairs[,2]), zip)] <- "parent-child"
relation_frame$ped_relation[match(paste(sibling_pairs[,1], sibling_pairs[,2]), zip)] <- "full-sib"
relation_frame$ped_relation[match(paste(sibling_pairs[,2], sibling_pairs[,1]), zip)] <- "full-sib"
relation_frame$ped_relation[match(paste(halfsib_pairs[,1], halfsib_pairs[,2]), zip)] <- "half-sib"
relation_frame$ped_relation[match(paste(halfsib_pairs[,2], halfsib_pairs[,1]), zip)] <- "half-sib"
relation_frame$ped_relation[match(paste(grandp_pairs[,1], grandp_pairs[,2]), zip)] <- "grandp-child"
relation_frame$ped_relation[match(paste(grandp_pairs[,2], grandp_pairs[,1]), zip)] <- "grandp-child"

plot(coefficient ~ ibs0, data = relation_frame, type = 'n', xlab = "IBS0", ylab = "Coefficient of Relatedness", bty = "n", xlim = c(0,11500))
grid()
points(coefficient ~ ibs0, data = subset(relation_frame, is.na(ped_relation)), cex = 1.5)
points(coefficient ~ ibs0, data = subset(relation_frame, ped_relation == "parent-child"), pch = 22, bg = 'red', cex = 1.5)
points(coefficient ~ ibs0, data = subset(relation_frame, ped_relation == "full-sib"), pch = 23, bg = 'darkgreen', cex = 1.5)
points(coefficient ~ ibs0, data = subset(relation_frame, ped_relation == "half-sib"), pch = 24, bg = 'darkorange', cex = 1.5)
points(coefficient ~ ibs0, data = subset(relation_frame, ped_relation == "grandp-child"), pch = 25, bg = 'orchid', cex = 1.5)
abline(h = 0)

legend("topright", pch = c(22, 23, 24, 25), pt.bg = c('red', 'darkgreen', 'darkorange', 'orchid'),
       legend = c("parent-child", "full-sib", "half-sib", "grandp-child"))

relation_frame[which(relation_frame$ped_relation == "full-sib" & relation_frame$coefficient < 0.4),]

# 

plot(ibs2 ~ ibs0, data = relation_frame, type = 'n', xlab = "IBS0", ylab = "IBS2", bty = "n", xlim = c(0,11500))
grid()
points(ibs2 ~ ibs0, data = subset(relation_frame, is.na(ped_relation)), cex = 1.5, col = 'darkgrey')
points(ibs2 ~ ibs0, data = subset(relation_frame, ped_relation == "parent-child"), pch = 22, bg = 'red', cex = 1.5)
points(ibs2 ~ ibs0, data = subset(relation_frame, ped_relation == "full-sib"), pch = 23, bg = 'darkgreen', cex = 1.5)
points(ibs2 ~ ibs0, data = subset(relation_frame, ped_relation == "half-sib"), pch = 24, bg = 'darkorange', cex = 1.5)
points(ibs2 ~ ibs0, data = subset(relation_frame, ped_relation == "grandp-child"), pch = 25, bg = 'orchid', cex = 1.5)
abline(h = 0)

legend("topright", pch = c(22, 23, 24, 25), pt.bg = c('red', 'darkgreen', 'darkorange', 'orchid'),
       legend = c("parent-child", "full-sib", "half-sib", "grandp-child"))

query_relation <- function(q1, q2) {
  return(rbind(
    relation_frame[relation_frame$id1 == q2 & relation_frame$id2 == q1,],
    relation_frame[relation_frame$id1 == q1 & relation_frame$id2 == q2,]))
}

## Analyze candidate sex contigs

psc_bycontig <- read.table("bcftools_PSC_bycontig.decap.txt", header = FALSE)[,c(3,4,5,6,10,14)]
names(psc_bycontig) <- c("ID", "nrefhom", "nnonrefhom", "nhet", "mean_depth", "nmiss")
contigs_100kb <- read.table("BGI_contigs_100kb.txt", header = FALSE)$V1[-1]
contigs_100kb_lengths <- read.table("BGI_contigs_100kb.txt", header = FALSE)$V2[-1]
names(contigs_100kb_lengths) <- contigs_100kb

# 18 total individuals
psc_bycontig$contig <- rep(contigs_100kb, each = 18)
psc_bycontig$n_sites <- psc_bycontig$nrefhom + psc_bycontig$nnonrefhom + psc_bycontig$nhet + psc_bycontig$nmiss
psc_bycontig$heterozygosity <- psc_bycontig$nhet / psc_bycontig$n_sites

#

x <- aggregate(psc_bycontig[,-7], by = list(psc_bycontig$ID), FUN = sum)
x$heterozygosity <- x$nhet / (sum(as.numeric(contigs_100kb_lengths)) - x$nmiss)


# get rid of contigs that have no variants
psc_bycontig <- subset(psc_bycontig, !contig %in% unique(psc_bycontig$contig[psc_bycontig$n_sites < 1]))

sex_id <- factor(ped_input$V5, levels = 1:2, labels = c("M","F"))
names(sex_id) <- ped_input$V2

psc_bycontig$sex <- sex_id[as.character(psc_bycontig$ID)]
psc_split <- split(psc_bycontig, psc_bycontig$contig)

m_idxlist <- lapply(psc_split, "[[", "sex") %>% lapply(., "==", "M") %>% lapply(., which)
f_idxlist <- lapply(psc_split, "[[", "sex") %>% lapply(., "==", "F") %>% lapply(., which)

male_depthlist <- lapply(psc_split, "[[", "mean_depth") %>% mapply("[", ., m_idxlist, SIMPLIFY = FALSE)
female_depthlist <- lapply(psc_split, "[[", "mean_depth") %>% mapply("[", ., f_idxlist, SIMPLIFY = FALSE)

male_hetlist <- lapply(psc_split, "[[", "heterozygosity") %>% mapply("[", ., m_idxlist, SIMPLIFY = FALSE)
female_hetlist <- lapply(psc_split, "[[", "heterozygosity") %>% mapply("[", ., f_idxlist, SIMPLIFY = FALSE)

dp_pvals <- mapply(wilcox.test, male_depthlist, female_depthlist, SIMPLIFY = FALSE) %>% sapply(., "[", "p.value")
het_pvals <- mapply(wilcox.test, male_hetlist, female_hetlist, SIMPLIFY = FALSE) %>% sapply(., "[", "p.value")
sigdp_idx <- which(dp_pvals < 0.05)
sighet_idx <- which(het_pvals < 0.05)

dp_ratio <- mapply("/",
                   lapply(psc_split, "[[", "mean_depth") %>% mapply("[", ., f_idxlist, SIMPLIFY = FALSE) %>% lapply(., mean),
                   lapply(psc_split, "[[", "mean_depth") %>% mapply("[", ., m_idxlist, SIMPLIFY = FALSE) %>% lapply(., mean))

m_het <- lapply(psc_split, "[[", "heterozygosity") %>% mapply("[", ., m_idxlist, SIMPLIFY = FALSE) %>% sapply(., mean)

plot(dp_ratio, m_het, ylim = c(0,1), ylab = "mean heterozygosity", xlab = "depth ratio")

# Union of significant hits
sigunion_idx <- union(sighet_idx, sigdp_idx)
points(dp_ratio[sigunion_idx], m_het[sigunion_idx], col = 'red')

Y_contigs <- names(which(dp_ratio < 0.5))

sig_contigs <- names(male_depthlist[sigunion_idx])
psc_sig <- psc_bycontig[psc_bycontig$contig %in% sig_contigs,]
psc_sig <- aggregate(psc_sig[c("mean_depth", "heterozygosity")],
                     by = list(sex = psc_sig$sex, contig = psc_sig$contig), FUN = mean)
psc_nonsig <- psc_bycontig[!(psc_bycontig$contig %in% sig_contigs),]
psc_agg <- aggregate(psc_nonsig[c("mean_depth", "heterozygosity")],
                     by = list(sex = psc_nonsig$sex, contig = psc_nonsig$contig), FUN = mean)

ggplot(subset(psc_sig, contig %in% names(which(dp_ratio > 1.5 | dp_ratio < 0.5))),
       aes(x = sex, y = mean_depth, group = contig)) + 
  geom_line() + geom_point(size = 2, aes(color = sex)) + theme_minimal()

ggplot(subset(psc_sig, contig %in% names(which(dp_ratio < 0.5))),
       aes(x = sex, y = heterozygosity, group = contig)) + 
  geom_line() + geom_point(size = 2, aes(color = sex)) + theme_minimal()

ggplot(data = psc_nonsig, aes(x = mean_depth, y = heterozygosity, color = "nonsignificant")) + geom_point(size = 1.5) + xlim(0,100) +
   geom_point(data = psc_sig, aes(x = mean_depth, y = heterozygosity, color = sex), size = 1.5) + theme_bw() +
   scale_color_manual(values = c("#f8766d", "#619cff", "gray")) + labs(color = "sex")

X_contigs <- names(which(dp_ratio > 1.5))

omit_these <- subset(psc_sig, (sex == "M" & mean_depth < 25))$contig
omit_idx <- which(psc_sig$contig %in% omit_these)
sig_include <- psc_sig[-omit_idx,]
sig_exclude <- psc_sig[omit_idx,]

ggplot(data = psc_nonsig, aes(x = mean_depth, y = heterozygosity, color = "nonsignificant")) + geom_point(size = 1.5) + xlim(0,80) +
   geom_point(data = sig_include, aes(x = mean_depth, y = heterozygosity, color = "unfiltered"), size = 1.5) +
   geom_point(data = sig_exclude, aes(x = mean_depth, y = heterozygosity, color = sex), size = 1.5) +
   theme_bw() + scale_color_manual(values = c("#f8766d", "#619cff", "grey75", "grey50")) + labs(color = "sex")

candidate_contigs <- contigs_100kb[!contigs_100kb %in% c(X_contigs, Y_contigs)]
write.table(candidate_contigs, file = "candidate_contigs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

                                                              MARGIN = c(1,2), function(x) { gsub("\\|", "/", x) })