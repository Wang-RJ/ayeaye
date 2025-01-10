# Generate a data frame for callability from tables from shell

hapsize_table <- read.table("callability/callable_sites.txt")
auto_size <- hapsize_table$V1
names(auto_size) <- paste("trio", 1:19, sep = "")

callhomM_tables <- paste("callability/", list.files("callability/")[grep("hom_callabilityGQ", list.files("callability/"))], sep = "")
callhet_tables <- paste("callability/", list.files("callability/")[grep("het_callabilityGQ", list.files("callability/"))], sep = "")

GQ_levels <- paste("GQ", seq(20,80,10), sep = "")

# Import tables in list form
Chet_list <- lapply(callhet_tables, function(het_table) {
  het <- read.table(het_table, header = FALSE)
  names(het) <- c("denom", "transmit", "dp", "dpgq", "bam", "ab30", "ab35", "ab40", "ab45", "ab50")
  return(het)
})
names(Chet_list) <- GQ_levels

ChomM_list <- lapply(callhomM_tables, function(hom_table, AD = 1) {
  homM <- read.table(hom_table, header = FALSE)
  names(homM) <- c("denom", "transmit", "dp", "dpgq", "bamad0")
  return(homM)
})
names(ChomM_list) <- GQ_levels

# trios 16-19 are the product of father-daughter inbreeding, and Chet cannot be calculated
# on them as in the other trios due to the absence of 0/0-1/1 in the parents
#
# trios 13-15 are the father and sister-mother unions leading to mothers in trios 16-19
# we'll use their mean values of Chet as are they seem the most comparable

Chet_list <- lapply(Chet_list, function(table) {
  table[16:19,] <- data.frame(t(colMeans(table[13:15,]))); return(table)
  })

dncAD0_noGQAB <- candidates_DP %>% filter_gq(., 20) %>% filter_bam %>%
  filter_bam_dpvar(., lower_dp_bounds, upper_dp_bounds) %>% filter_rel
dncAD0_noGQAB <- dncAD0_noGQAB[!names(dncAD0_noGQAB) %in% c("trio2", "trio3")]

# Other questionable trios:
# trio1 has 1X2124 as mother, trio13 has 1x0356 as mother, trio8 has 15444 as proband
# though only trio2 and trio3 are from 1X3656

denovo_candidatesAD0 <- vector("list", 5)
names(denovo_candidatesAD0) <- as.character(seq(30, 50, 5))

for(AB_idx in as.character(seq(30, 50, 5))) {
  for(GQ in seq(20,80,10)) {
    denovo_candidatesAD0[[AB_idx]][[paste("GQ", GQ, sep = "")]] <- dncAD0_noGQAB %>% filter_gq(., GQ) %>% filter_allelicbalance(as.numeric(AB_idx)/100, 1) %>% bind_candidates
  }
}

callability_M <- list(AD0 = list())
for(AB_idx in as.character(seq(30, 50, 5))) {
  callability_M[["AD0"]][[AB_idx]] <- mapply(function(chet, chom) { chet[[paste("ab", AB_idx, sep = "")]] / chet[["dp"]] *
                                                                    (chom[["bamad0"]] / chom[["transmit"]])**2 },
                                             Chet_list, ChomM_list)
}

ratebytrio <- function(dn_list, callability) {
  lapply(GQ_levels, function(dnf_idx) {
    dnf <- dn_list[[dnf_idx]]
    
    n_mut <- sapply(split(dnf, dnf$trio), nrow)
    t_ord <- gsub("trio", "", names(n_mut)) %>% as.numeric
    n_mut <- n_mut[order(t_ord)]
    
    rate <- n_mut / (2 * auto_size[sort(t_ord)] * callability[sort(t_ord),dnf_idx])
    return(rate)
  })
}

for(AB_idx in as.character(seq(30,50,5))) {
  print(sapply(ratebytrio(denovo_candidatesAD0[[AB_idx]], callability_M[["AD0"]][[AB_idx]]), mean))
}

cbplot <- function(dn_frame, callability) {
  dn_frame <- lapply(dn_frame, function(df) { return(subset(df, !trio %in% c("trio2", "trio3"))) })
  callability_data <- data.frame(GQ = seq(20,80,10),
                                 mutations = unlist(lapply(dn_frame, nrow)),
                                 meanCallability = colSums(callability) / 17,
                                 meanCallableSites = t((auto_size %*% callability) / (17 * 1e9)),
                                 meanMutationRate = sapply(ratebytrio(dn_frame, callability), mean) * 1e8)
  callability_data$CI <- 1.96 * sqrt((callability_data$mutations / (2 * callability_data$meanCallableSites)) / (2 * callability_data$meanCallableSites)) / 100
  p1 <- ggplot(callability_data, aes(x = GQ, y = mutations)) + geom_line() + geom_point() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p2 <- ggplot(callability_data, aes(x = GQ, y = meanCallability)) + ylab("Callability") + geom_line() + geom_point() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p3 <- ggplot(callability_data, aes(x = GQ, y = meanCallableSites)) + ylab("Callable Sites") + geom_line() + geom_point() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p4 <- ggplot(callability_data, aes(x = GQ, y = meanMutationRate)) + ylab("Mutation Rate") + geom_point() + theme_bw() + geom_errorbar(aes(ymin = meanMutationRate - CI,
                                                                                                                                            ymax = meanMutationRate + CI),
                                                                                                                                        stat = "identity", width = 0.4)
  GGp1 <- ggplotGrob(p1)
  GGp2 <- ggplotGrob(p2)
  GGp3 <- ggplotGrob(p3)
  GGp4 <- ggplotGrob(p4)
  grid::grid.newpage()
  grid::grid.draw(rbind(GGp1, GGp2, GGp3, GGp4))
}

cbplot(denovo_candidatesAD0[["35"]], callability_M[["AD0"]][["35"]])

cbplot_ab <- function() {
  ab_seq = as.character(seq(30,50,5))
  callability_data <- data.frame(AB = as.numeric(ab_seq),
                                 mutations = sapply(ab_seq, function(x) { nrow(denovo_candidatesAD0[[x]][["GQ60"]]) }),
                                 meanCallability = colSums(sapply(ab_seq, function(x) { callability_M[["AD0"]][[x]][,"GQ60"] })) / 11,
                                 meanCallableSites = t((auto_size %*% sapply(ab_seq, function(x) { callability_M[["AD0"]][[x]][,"GQ60"] })) / (11* 1e9)),
                                 meanMutationRate = sapply(ab_seq, function(x) {
                                   sapply(ratebytrio(denovo_candidatesAD0[[x]], callability_M[["AD0"]][[x]])[6], mean) }) * 1e8)
  callability_data$CI <- 1.96 * sqrt((callability_data$mutations / (2 * callability_data$meanCallableSites)) / (2 * callability_data$meanCallableSites)) / 100
  p1 <- ggplot(callability_data, aes(x = AB, y = mutations)) + geom_line() + geom_point() + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p2 <- ggplot(callability_data, aes(x = AB, y = meanCallability)) + ylab("Callability") + geom_line() + geom_point() + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p3 <- ggplot(callability_data, aes(x = AB, y = meanCallableSites)) + ylab("Callable Sites") + geom_line() + geom_point() + theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p4 <- ggplot(callability_data, aes(x = AB, y = meanMutationRate)) + ylab("Mutation Rate") + geom_point() + geom_errorbar(aes(ymin = meanMutationRate - CI,
                                                                                                                               ymax = meanMutationRate + CI),
                                                                                                                           stat = "identity", width = 0.4) +
    theme_bw() #+ ylim(.5, 1.25)
  GGp1 <- ggplotGrob(p1)
  GGp2 <- ggplotGrob(p2)
  GGp3 <- ggplotGrob(p3)
  GGp4 <- ggplotGrob(p4)
  grid::grid.newpage()
  grid::grid.draw(rbind(GGp1, GGp2, GGp3, GGp4))
}

cbplot_ab()

# Rename in global environment for plots
denovo_candidates <- denovo_candidatesAD0[["35"]][["GQ60"]]
callability <- callability_M[["AD0"]][["35"]][-c(2:3),]
auto_size_bytrio <- auto_size[-c(2,3)]
