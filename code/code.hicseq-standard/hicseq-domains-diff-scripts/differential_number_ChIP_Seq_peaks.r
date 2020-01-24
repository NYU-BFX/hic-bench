# differential_number_ChIP_Seq_peaks.r script

# Load libraries
library(plyr)
library(reshape2)
library(ggplot2)

# Read inputs
args <- commandArgs()

chip1_active_file <- args[3]
chip2_active_file <- args[4]
chip1_inactive_file <- args[5]
chip2_inactive_file <- args[6]
chip1_unchanged_file <- args[7]
chip2_unchanged_file <- args[8]

obj1 <- args[9]
obj2 <- args[10]
bin.size <-  as.numeric(args[11])
out_prefix <- args[12]
analysis_type <- args[13]
num_peaks1 <- as.numeric(args[14])
num_peaks2 <- as.numeric(args[15])
num_common_TADs <- as.numeric(args[16])
domains_common_file <- args[17]
active_TADs_file <- args[18]
inactive_TADs_file <- args[19]
unchanged_TADs_file <- args[20]


chip1_active <- read.table(chip1_active_file, header = FALSE, sep = "\t")
chip2_active <- read.table(chip2_active_file, header = FALSE, sep = "\t")
chip1_inactive <- read.table(chip1_inactive_file, header = FALSE, sep = "\t")
chip2_inactive <- read.table(chip2_inactive_file, header = FALSE, sep = "\t")
chip1_unchanged <- read.table(chip1_unchanged_file, header = FALSE, sep = "\t")
chip2_unchanged <- read.table(chip2_unchanged_file, header = FALSE, sep = "\t")
domains_common <- read.table(domains_common_file, header = FALSE, sep = "\t")
active_TADs <- read.table(active_TADs_file, header = FALSE, sep = "\t")
inactive_TADs <- read.table(inactive_TADs_file, header = FALSE, sep = "\t")
unchanged_TADs <- read.table(unchanged_TADs_file, header = FALSE, sep = "\t")

colnames(domains_common) <- c("chr", "start", "end")
colnames(active_TADs) <- c("chr", "start", "end")
colnames(inactive_TADs) <- c("chr", "start", "end")
colnames(unchanged_TADs) <- c("chr", "start", "end")

# Keep the first 6 columns of the files
chip1_active <- chip1_active[,1:6]	
chip2_active <- chip2_active[,1:6]	
chip1_inactive <- chip1_inactive[,1:6]	
chip2_inactive <- chip2_inactive[,1:6]	
chip1_unchanged <- chip1_unchanged[,1:6]	
chip2_unchanged <- chip2_unchanged[,1:6]

colnames(chip1_active) <- c("chr", "start", "end", "chr_peak", "start_peak", "end_peak")
colnames(chip2_active) <- c("chr", "start", "end", "chr_peak", "start_peak", "end_peak")
colnames(chip1_inactive) <- c("chr", "start", "end", "chr_peak", "start_peak", "end_peak")
colnames(chip2_inactive) <- c("chr", "start", "end", "chr_peak", "start_peak", "end_peak")
colnames(chip1_unchanged) <- c("chr", "start", "end", "chr_peak", "start_peak", "end_peak")
colnames(chip2_unchanged) <- c("chr", "start", "end", "chr_peak", "start_peak", "end_peak")

# Count peaks per TAD for both objects

num_peaks_per_TAD_chip1_active <- count(chip1_active, c("chr", "start", "end"))
num_peaks_per_TAD_chip2_active <- count(chip2_active, c("chr", "start", "end"))
num_peaks_per_TAD_chip1_inactive <- count(chip1_inactive, c("chr", "start", "end"))
num_peaks_per_TAD_chip2_inactive <- count(chip2_inactive, c("chr", "start", "end"))
num_peaks_per_TAD_chip1_unchanged <- count(chip1_unchanged, c("chr", "start", "end"))
num_peaks_per_TAD_chip2_unchanged <- count(chip2_unchanged, c("chr", "start", "end"))

# Count average TAD length for differentially active TADs and all common TADs
domains_common$length <- domains_common$end - domains_common$start
mean_length_domains <- mean(domains_common$length)

active_TADs$length <- active_TADs$end - active_TADs$start
mean_length_active <- mean(active_TADs$length)

inactive_TADs$length <- inactive_TADs$end - inactive_TADs$start
mean_length_inactive <- mean(inactive_TADs$length)

unchanged_TADs$length <- unchanged_TADs$end - unchanged_TADs$start
mean_length_unchanged <- mean(unchanged_TADs$length)

# Normalize number of peaks per total number of peaks per sample overlapping the common TADs
# the total number of common TADs and the average TAD length
num_peaks_per_TAD_chip1_active$enrich <- (num_peaks_per_TAD_chip1_active$freq * mean_length_active) / ((num_peaks1 * mean_length_domains)/num_common_TADs)
num_peaks_per_TAD_chip2_active$enrich <- (num_peaks_per_TAD_chip2_active$freq * mean_length_active) / ((num_peaks2 * mean_length_domains)/num_common_TADs)

num_peaks_per_TAD_chip1_inactive$enrich <- (num_peaks_per_TAD_chip1_inactive$freq * mean_length_inactive) / ((num_peaks1 * mean_length_domains)/num_common_TADs)
num_peaks_per_TAD_chip2_inactive$enrich <- (num_peaks_per_TAD_chip2_inactive$freq * mean_length_inactive) / ((num_peaks2 * mean_length_domains)/num_common_TADs)

num_peaks_per_TAD_chip1_unchanged$enrich <- (num_peaks_per_TAD_chip1_unchanged$freq * mean_length_unchanged) / ((num_peaks1 * mean_length_domains)/num_common_TADs)
num_peaks_per_TAD_chip2_unchanged$enrich <- (num_peaks_per_TAD_chip2_unchanged$freq * mean_length_unchanged) / ((num_peaks2 * mean_length_domains)/num_common_TADs)

# Calculate p-values between the distributions for each TAD category
active_TADs_p_value <- t.test(as.numeric(num_peaks_per_TAD_chip1_active$enrich), as.numeric(num_peaks_per_TAD_chip2_active$enrich), paired = FALSE, alternative = "two.sided")
inactive_TADs_p_value <- t.test(as.numeric(num_peaks_per_TAD_chip1_inactive$enrich), as.numeric(num_peaks_per_TAD_chip2_inactive$enrich), paired = FALSE, alternative = "two.sided")
unchanged_TADs_p_value <- t.test(as.numeric(num_peaks_per_TAD_chip1_unchanged$enrich), as.numeric(num_peaks_per_TAD_chip2_unchanged$enrich), paired = FALSE, alternative = "two.sided")

# Calculate mean values
mean_active_chip1 <- mean(num_peaks_per_TAD_chip1_active$enrich)
mean_active_chip2 <- mean(num_peaks_per_TAD_chip2_active$enrich)
mean_inactive_chip1 <- mean(num_peaks_per_TAD_chip1_inactive$enrich)
mean_inactive_chip2 <- mean(num_peaks_per_TAD_chip2_inactive$enrich)
mean_unchanged_chip1 <- mean(num_peaks_per_TAD_chip1_unchanged$enrich)
mean_unchanged_chip2 <- mean(num_peaks_per_TAD_chip2_unchanged$enrich)

mean_p_values <- data.frame("active" = c(mean_active_chip1, mean_active_chip2, active_TADs_p_value$p.value),
                                "inactive" = c(mean_inactive_chip1, mean_inactive_chip2, inactive_TADs_p_value$p.value),
                                "unchanged" = c(mean_unchanged_chip1, mean_unchanged_chip2, unchanged_TADs_p_value$p.value))
rownames(mean_p_values) <- c("mean_chip1 ","mean_chip2", "pvalue")

write.table(mean_p_values, paste(out_prefix, analysis_type, "mean_values.tsv", sep = "_"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

colnames(num_peaks_per_TAD_chip1_active)[5] <- obj1
colnames(num_peaks_per_TAD_chip2_active)[5] <- obj2

colnames(num_peaks_per_TAD_chip1_inactive)[5] <- obj1
colnames(num_peaks_per_TAD_chip2_inactive)[5] <- obj2

colnames(num_peaks_per_TAD_chip1_unchanged)[5] <- obj1
colnames(num_peaks_per_TAD_chip2_unchanged)[5] <- obj2

# Merge for the same TAD category 
active_TADs_peaks <- merge(num_peaks_per_TAD_chip1_active, num_peaks_per_TAD_chip2_active, 
                           by = c("chr", "start", "end"), all = TRUE)
inactive_TADs_peaks <- merge(num_peaks_per_TAD_chip1_inactive, num_peaks_per_TAD_chip2_inactive, 
                           by = c("chr", "start", "end"), all = TRUE)
unchanged_TADs_peaks <- merge(num_peaks_per_TAD_chip1_unchanged, num_peaks_per_TAD_chip2_unchanged, 
                           by = c("chr", "start", "end"), all = TRUE)

# Format data for plotting
active_TADs_peaks_melted <- melt(active_TADs_peaks[,c(1,2,3,5,7)], id.vars =  c("chr", "start", "end"))
inactive_TADs_peaks_melted <- melt(inactive_TADs_peaks[,c(1,2,3,5,7)], id.vars =  c("chr", "start", "end"))
unchanged_TADs_peaks_melted <- melt(unchanged_TADs_peaks[,c(1,2,3,5,7)], id.vars =  c("chr", "start", "end"))

active_TADs_peaks_melted$category <-"active TADs"
inactive_TADs_peaks_melted$category <-"inactive TADs"
unchanged_TADs_peaks_melted$category <-"unchanged TADs"

all_TADs_peaks_melted <- rbind(active_TADs_peaks_melted, inactive_TADs_peaks_melted, unchanged_TADs_peaks_melted)

pdf(file = paste(out_prefix, analysis_type, obj1, obj2, "boxplot_number_of_peaks_per_TAD_norm.pdf", sep = "_"), width = 7, height = 5)
ggplot(all_TADs_peaks_melted, aes(x = variable, y = value, fill = variable)) + geom_boxplot() + facet_grid(category ~ .) +
  theme_bw() + scale_fill_discrete(name="Cell Type") + xlab(paste(analysis_type, "Enrichment per TAD", sep = " "))+ ylab("")+
  scale_fill_manual(values=c("dark red","blue")) + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") + coord_flip()
dev.off()

