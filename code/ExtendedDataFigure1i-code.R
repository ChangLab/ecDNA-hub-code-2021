library(RColorBrewer)
library(ggplot2)
library(ggpubr)
source("data/custom_colors.R")

set.seed(1)
colo.df <- read.delim("data/COLO320DM_ecDNAcount_vs_bursting.txt")

# make list containings vectors of how many ecDNAs are transcribing for ecDNA cluster of given size
colo.dist <- list()
for (i in unique(colo.df$ecDNA_number)) {
    colo.dist[[as.character(i)]] <- c(
        rep(1, sum(colo.df$RNA_count[colo.df$ecDNA_number == i]))
        , rep(0, sum(colo.df$ecDNA_number[colo.df$ecDNA_number == i]) -
              sum(colo.df$RNA_count[colo.df$ecDNA_number == i])))
}

# resample cluster size-matched distribution to calculate bursting frequency for each ecDNA cluster
colo.df.resample <- data.frame()
for (i in colo.df$ecDNA_number) {
    dist.i <- colo.dist[[as.character(i)]]
    n <- length(dist.i)
    colo.df.resample <- rbind(colo.df.resample, c(i, sum(sample(dist.i, n, replace = TRUE)) / n))
}
colnames(colo.df.resample) <- c('ecDNA_number', 'bursting_frequency')

# classify hubs as singleton or multiple
colo.df.resample$single <- ifelse(colo.df.resample$ecDNA_number == 1,  "Single ecDNA", "ecDNA in hub")

#plotting
ggplot(colo.df.resample, aes(x = single, y = bursting_frequency, fill = single)) + 
    geom_violin(alpha = 0.6) + 
    geom_boxplot(width = 0.2) + 
    stat_compare_means(comparisons = list(c("Single ecDNA", "ecDNA in hub"))) +
    theme_classic() + 
    scale_fill_manual(values = my_custom_palettes$splash_of_salmon[c(4,5)]) +
    ylab("MYC transcriptional probability") + xlab("") + theme(legend.position = "none")