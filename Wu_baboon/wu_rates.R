# Quick calculation of mutation rates from baboons in Wu et al.

wu_mutations <- read.csv("Wu_DNMs.csv", check.names = FALSE)
wu_mutbyID <- table(wu_mutations$F1)

wu_trios <- read.csv("Wu_trios.csv", check.names = FALSE)
wu_trios <- subset(wu_trios, F1 %in% names(wu_mutbyID))

ids <- names(wu_mutbyID)
trio_idxorder <- match(ids, wu_trios$F1)

wu_rates <- wu_mutbyID / wu_trios[trio_idxorder,][["Callable Genome Size"]]
wu_rates <- rates / (1 - wu_trios[trio_idxorder,][["FNR"]])

wu_patage <- wu_trios[trio_idxorder,][["Father Age At Birth (y)"]]

plot(as.vector(wu_rates) ~ wu_patage, ylim = c(1e-8, 2.75e-8), xlim = c(4,16), xaxt = 'n', pch = 16, ylab = "mutation rate", xlab = "paternal age at conception")
points(candidate_count / (2e9 * 1-0.8) ~ paternalageatbirth_years, bg = 'red', pch = 21, cex = 1.25)
axis(1, at = seq(0,15,5))
