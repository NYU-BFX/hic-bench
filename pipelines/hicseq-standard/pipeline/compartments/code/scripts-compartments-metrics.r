argv = commandArgs(trailingOnly = TRUE)
HKcounts = argv[1L]
TSScounts = argv[2L]
sampleName = argv[3L]

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

# Read Data
bed_cols = c("chr","start","end","pc1", "counts")
pc1.counts = fread(HKcounts, col.names = bed_cols)
setnames(pc1.counts, "counts", "HK.counts")
pc1.counts = pc1.counts[fread(TSScounts, col.names = bed_cols),
                        TSS.counts := i.counts, on = bed_cols[1:3]]
pc1.counts[, `:=`(group = sign(pc1), size = end - start)]
chromosomes = stringr::str_sort(unique(pc1.counts$chr), numeric = TRUE)
pc1.counts[, chr := factor(chr, chromosomes)]

# Aggregate Info
info.df = pc1.counts[, lapply(.SD, sum), c("chr", "group"),
                     .SDcols = c("HK.counts", "TSS.counts", "size")]
info.df[, `:=`(size.per.count.HK = size / HK.counts,
               size.per.count.TSS = size / TSS.counts)]

# Compute logRatio
ratio.df = info.df[, lapply(.SD, function(x) log(x[2]) - log(x[1])), chr,
                   .SDcols = colnames(info.df)[-c(1, 2)]]
ratio.df[is.nan(size.per.count.HK), size.per.count.HK := 0]
ratio.df[is.nan(size.per.count.TSS), size.per.count.TSS := 0]
ratio.df[, invertSign := fifelse(abs(size.per.count.HK) > abs(size.per.count.TSS),
                                 -sign(size.per.count.HK),
                                 -sign(size.per.count.TSS))]

# Invert Signs
# pc1.counts[ratio.df, pc1 := pc1 * i.invertSign, on = "chr"]
info.df[ratio.df, invertSign := i.invertSign, on = "chr"]

# Write metrics
ratio.df[, group := 2L]
info.df = rbind(info.df, ratio.df)[order(chr, group)]  # add log ratios to metrics
info.df[, `:=`(group = factor(group, c(-1, 1, 2), c(-1, 1, "log(ratio[1:-1])")),
               sampleName = sampleName,
               invertSign = invertSign < 0)]
setnames(info.df, c("chr", "group", "size"), c("chromosome", "signGroup", "size.bp"))
setnames(ratio.df, c("chr", "group", "size"), c("chromosome", "signGroup", "size.bp"))
info.df[, lapply(.SD, round, digits = 3), .SDcols = which(sapply(info.df, is.numeric))]
setcolorder(info.df, c("chromosome", "signGroup", "HK.counts", "TSS.counts", "size.bp",
                       "size.per.count.HK", "size.per.count.TSS", "invertSign", "sampleName"))
fwrite(info.df, "pc1_metrics_summary.txt", quote = FALSE, sep = "\t", scipen = 999)

# Figures
metrics = c("size.bp", "size.per.count.HK", "size.per.count.TSS")
df.gg3 = melt(ratio.df, id.vars = "chr", measure.vars = metrics,
              variable.name = "metric", value.name = "ratioMinusPlus")  # wide -> long
df.gg3[, metric := factor(metric, metrics)]  # respect order
p = ggplot(df.gg3,aes(chr, ratioMinusPlus)) +
  geom_point(aes(color = metric), size = 3) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
  labs(x = "", y = "log ratio [A:B]")
ggsave('metrics_AB.logratio.pdf', p, width = 12, useDingbats = FALSE)
