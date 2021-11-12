# FUNCTIONS ---------------------------------------------------------------

getEigenMat <- function(pca1.matrix, centrotelo) {
    # filter out centro- and telo-meres
    pca1.matrix <- makeGRangesFromDataFrame(pca1.matrix,
                                            keep.extra.columns = TRUE,
                                           starts.in.df.are.0based <- TRUE)
    pca1.matrix <- sort(subsetByOverlaps(pca1.matrix, centrotelo, invert = TRUE))
    pca1.matrix <- as.data.table(pca1.matrix)[, strand := NULL]
    setnames(pca1.matrix, c("seqnames", "width"), c("chr", "size"))

    # get matrix in long format
    df.sizes <- melt(pca1.matrix, id.vars = c("chr", "start", "end", "size"),
                  variable.name = "sampleName", value.name = "PC1")
    df.sizes[, `:=`(comp = fifelse(PC1 > 0, "A", "B"),
                  size = size / 1e6)]
    # aggregate metrics per sample & compartment
    df.sizes <- df.sizes[!is.na(PC1) & PC1 != 0,
                         .(size = sum(size), meanPC1 = mean(PC1), medianPC1 = median(PC1)),
                         c("sampleName", "comp")]
    df.sizes[, percSize := size / sum(size), sampleName]

    return(list("mat"=pca1.matrix, "df.sizes"=df.sizes))
}

getBarPlot <- function(data, mapping, grid = TRUE){
    if (grid) {
        # a bit ugly but... get df.sizes from global env
        sz <- 3
        mapping <- lapply(mapping, rlang::quo_get_expr)
        mapping <- sapply(mapping, rlang::as_string)
        df.sizes.gg <- df.sizes[sampleName %in% mapping]
    } else {
        sz <- 5
        df.sizes.gg <- data
    }
    df.sizes.gg[, sizeLabel := sprintf("%.0f Mb\n%s (%s)", size, comp,
                                 scales::percent(percSize, .1))]

    p <- ggplot(data = df.sizes.gg, aes(sampleName, size)) +
        geom_col(aes(fill = comp), color = "black", width = 0.6) +
        geom_text(aes(label = sizeLabel),colour="white", size = sz,
              position = position_stack(vjust = 0.5))+
        scale_fill_manual(values = c("red", "blue")) +
        scale_x_discrete(labels = labeller_vec)+
        labs(x = "", y = "compartments size (Mb)", fill = "Compartment")

    if (grid) {
        p <- p + theme(axis.text.x = element_blank(),
                     legend.position = "none")
    } else {
        # flip for easier readability
        p <- p + coord_flip() + theme_minimal(base_size = 14)
    }
    p
}

getDensPlot <- function(data, mapping, grid = TRUE) {
    # get vars
    mapping <- lapply(mapping, rlang::quo_get_expr)
    mapping <- sapply(mapping, rlang::as_string)
    mapping <- mapping[c("x", "y")]
    # wide -> long
    temp <- as.data.table(data[, mapping])
    setnames(temp, labeller_vec[mapping])
    temp <- melt(temp, measure.vars = 1:2, variable.name = "group", value.name = "PC1")
    # plot
    ggplot(data = temp, aes(x=PC1, fill=group)) +
        geom_density(show.legend = TRUE, alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dotdash", color = "darkred") +
        scale_fill_manual(values = c("green", "purple"))+
        labs(fill = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"))
}

getHexPlot <- function(data, mapping, grid = TRUE, nbins = 70){
    sz <- ifelse(grid, 3, 5)
    mat_range <- range(data, na.rm = TRUE)
    # get vector to compute correlation
    xv <- rlang::eval_tidy(mapping$x, data)
    yv <- rlang::eval_tidy(mapping$y, data)

    ggplot(data, mapping) +
        stat_binhex(bins = nbins) +
        geom_abline(slope = 1, intercept = 0, linetype = 2) +
        annotate(geom = "text", label = sprintf("r = %.3g", cor(xv, yv)),
                 x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, size = sz) +
        scale_fill_gradientn(colours=c("blue", "orange", "red"), trans="log10") +
        labs(fill = "log10 count\n(eigenvector-1)\n") +
        coord_cartesian(xlim = mat_range, ylim = mat_range) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"))
}

getEigenPCA <- function(mat, df.color) {
    pca <- prcomp(t(mat))
    names <- colnames(mat)
    autoplot(pca, label = FALSE, size=5, col = df.color$sampleColor) +
        geom_point(size=2) +
        ggrepel::geom_text_repel(aes(label = df.color$sampleName), segment.size = 0.5,
                        force = 3, size = 4, point.padding = 1, segment.alpha = 1) +
        theme(panel.background = element_rect(fill = "white"),
              panel.border = element_rect(colour = "black",fill=NA,size = 2),
              axis.text.x = element_text(size = 20,colour = "black"),
              axis.text.y = element_text(size = 20,colour = "black"),
              axis.title.x = element_text(size = 20,colour = "black"),
              axis.title.y = element_text(size = 20,colour = "black"),
              legend.position = "none")
}

getEigenHmap <- function(mat, bins.fraction = 0.3){
    ix <- order(matrixStats::rowSds(mat), decreasing = TRUE)
    mat <- mat[ix, ]
    #tot.bins = number of top most variable pc1 bins to plot
    tot.bins <- round(nrow(mat) * bins.fraction)

    gplots::heatmap.2(mat[1:tot.bins, , drop = FALSE],
              main = "",
              col = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
              dendrogram = "both",
              scale = "none",
              hclustfun = function(x) hclust(x, method = "ward.D2"),
              distfun = function(x) dist(x, method = "manhattan"),
              labRow = FALSE,
              cexRow=1,
              cexCol=1,
              margins=c(12,7),
              trace="none",
              srtCol=45,
              key = TRUE,
              keysize = 1.5,
              density.info = "none")
    title(paste(scales::percent(bins.fraction, 1),
                "most variable genomic bins\n(Eigenvector-1)"), cex.main = 0.7)
}

getEigenPairs <- function(mat) {
    ggpairs(mat,
            lower = list(continuous = getHexPlot),
            diag = list(continuous = getBarPlot),
            upper = list(continuous = getDensPlot),
            legend = c(2, 1),
            progress = FALSE) +
        theme_bw(base_size=14)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text=element_text(size=8, color="black"),
              axis.line = element_line(color = "black"),
              legend.direction = "vertical",
              legend.position = "left",
              strip.background = element_rect(fill="grey"))
}

getMetricsPlot <- function(DT, metrics) {
        DT <- melt(DT, id.vars = c("chromosome", "sampleName"), measure.vars = metrics,
                  variable.name <- "metric", value.name = "logratio.AB")
        DT <- DT[is.finite(logratio.AB) & !is.na(logratio.AB)]
        DT[, metric := factor(metric, metrics)]
        p <- ggplot(DT, aes(metric, -logratio.AB)) +
            geom_boxplot(aes(fill = metric)) +
            ggbeeswarm::geom_beeswarm(alpha = .25) +
            geom_hline(yintercept = 0, linetype = 4, color = "darkred") +
            facet_wrap(~chromosome, nrow = 1L) +
            labs(x = "", y = "log(A/B)") +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_blank(),
                  axis.line = element_line(colour = "black"))
}

# classify switched pc1 bins (pairwise-comparison)
classifyPc1Bins <- function(p, mat, delta.cut) {
    df.mat <- mat[, c("chr", "start", "end", "size", p), with = FALSE]
    setnames(df.mat, p, c("X", "Y")) # easier notation
    df.mat[, `:=`(pc1.diff = Y - X,
                  delta = (Y - X) / abs(Y),  # this is asymmetric
                  compX = fifelse(X > 0, "A", "B"),
                  compY = fifelse(Y > 0, "A", "B"))
          ]
    df.mat[, `:=`(switch = paste0(compX, fifelse(abs(delta) < delta.cut, compX, compY)),
                  compX = NULL, compY = NULL)  # no longer needed
          ]
    df.mat[X == 0 | Y == 0, switch := NA_character_]
    df.mat
}

fix_sample_names <- function(labelsMat) {
    if (any(grepl("^\\d", labelsMat))) labelsMat <- paste0("S", labelsMat)  #to avoid ggplot error when first character of label is numeric
    labelsMat <- gsub(labelsMat, pattern = "-", replacement = "_")
    labelsMat <- gsub(labelsMat, pattern = ":", replacement = ".")
}

read_pca_matrix <- function (path) {
    DT <- fread(path)
    DT[, c("chr", "start", "end") := tstrsplit(V1, split = "[:-]", type.convert = TRUE)][, V1 := NULL]
    setcolorder(DT, c("chr", "start", "end"))
    labelsMat <- setdiff(names(DT), bed.cols)
    setnames(DT, labelsMat, fix_sample_names(labelsMat))
    DT
}

read_figure_labels <- function(sample_labels, labelsMat) {
    sampleLabels <- data.table(sampleLabels = scan(sample_labels, what = "", quiet = TRUE))  # sampleLabels format = groupName:sampleName
    sampleLabels[, c("groupName", "sampleName") := tstrsplit(sampleLabels, split = ":", fixed = TRUE)]
    sampleLabels[groupName == "", groupName := sampleName]
    sampleLabels[, groupName := forcats::fct_inorder(groupName)]  # respect user order
    sampleLabels[, sampleName := fix_sample_names(sampleName)]
    # check if labels match with matrix columns
    if (all(labelsMat %in% sampleLabels$sampleName)) {
        sampleLabels = sampleLabels[match(sampleName, labelsMat)]
    } else if (nrow(sampleLabels) != length(labelMats)) {
        stop("Number of labels does not match number of samples")
    } # now by assumption sampleLabels order == labelsMat
    droplevels(sampleLabels$groupName)
}

# load libraries
suppressPackageStartupMessages({
    library("data.table")
    library("RColorBrewer")
    library("ggplot2")
    library("GGally")
    library("ggfortify")  # autoplot
    library("GenomicRanges")
    library("optparse")
})


# PARSE ARGUMENTS ---------------------------------------------------------

# process command-line arguments
cmdline_args <- commandArgs(trailingOnly=T);

# TODO: rewrite argument parsing. Does not make much sense
option_list <- list(
  make_option(c("-m","--pca1-matrix"), help="A matrix containing the samples PC1 values for each sample [default \"%default\"]."),
  make_option(c("-o","--output-dir"), help="Output directory (required) [default \"%default\"]."),
  make_option(c("-c","--centrotelo-file"), help="A bed file containing centromere and telomere coordinates [default \"%default\"]."),
  make_option(c("-L","--sample-labels"), default="", help="A file containing sample labels in its first column [default \"%default\"]."),
  make_option(c("-r","--metrics-file"), default="", help="A text file containing the AB logratio metrics for each sample [default \"%default\"]."),
  make_option(c("-d","--delta-cutoff"), default=1.2, help="A delta value cutoff for compartment switch call [default \"%default\"]."),
  make_option(c("--show-text"), action="store_true",default=TRUE, help="Display sample label text on PCA plot."),
  make_option(c("--use-short-names"), action="store_true",default=FALSE, help="Use only second part of the name (after the colon)."),
  make_option(c("--plain"), action="store_true",default=FALSE, help="Create only PCA and heatmap plots.")
)
usage <- "scripts-compartments-stats.r [OPTIONS] MATRIX";

# get command line options & input arguments
arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
opt <- arguments$options;
files <- opt[["pca1-matrix"]]
metrics_file <- opt[["metrics-file"]]
out_dir <- opt[["output-dir"]]
delta.cut <- as.numeric(opt[["delta-cutoff"]])
sample_labels <- opt[["sample-labels"]]
# show_text <- opt[["show-text"]]
use_short_names <- opt[["use-short-names"]]
plain <- opt[["plain"]]
centrotelo_file <- opt[["centrotelo-file"]]


# DATA PROCESSING ----------------------------------------------------------

# create output directories
if (!dir.exists(out_dir)) dir.create(out_dir) else warning("output directory already exists, it will be overwritten!", call. = FALSE)


### LOAD DATA ###

bed.cols <- c("chr", "start", "end")
centrotelo <- makeGRangesFromDataFrame(fread(centrotelo_file, col.names = bed.cols))
pca1.matrix <- read_pca_matrix(files)
sampleNames <- setdiff(colnames(pca1.matrix), bed.cols)

## PROCESS DATA ###

comparisons <- t(combn(sampleNames, 2L))
matAndSizes <- getEigenMat(pca1.matrix, centrotelo) # filter pca1 matrix and get compartment sizes
mat <- as.matrix(matAndSizes$mat[, 5:ncol(matAndSizes$mat)])
df.sizes <- matAndSizes$df.sizes

# Compartment Switches
df.comparisons <- apply(comparisons, 1L, classifyPc1Bins,
                        mat = matAndSizes$mat, delta.cut = delta.cut)
names(df.comparisons) <- apply(comparisons, 1, paste, collapse = ":")
df.switch <- rbindlist(df.comparisons, idcol = "comparisonID")
df.switch <- df.switch[!is.na(switch),
                      .(size = sum(size)/1e6, N = .N),
                      c("comparisonID", "switch")]
df.switch[, switch := factor(switch, c("AB", "AA", "BB", "BA"))]

# set color groups
sampleLabels <- if (use_short_names && file.exists(sample_labels)) read_figure_labels(sample_labels, sampleNames) else factor(sampleNames)
groupColor <- c(brewer.pal(7,"Set1"), brewer.pal(7,"Set2"), brewer.pal(7,"Set3"))[1:nlevels(sampleLabels)]
groupColor <- setNames(groupColor, levels(sampleLabels))
sampleLabels <- data.table(sampleName = sampleNames,
                           sampleLabel = sampleLabels,
                           sampleColor = groupColor[as.character(sampleLabels)])
labeller_vec <- sampleLabels[, setNames(as.character(sampleLabel), sampleName)]


# MAIN --------------------------------------------------------------------

if (ncol(mat) > 2) {
    p <- getEigenPCA(mat, df.color = sampleLabels)
    ggsave(file.path(out_dir, "eigenPCA.pdf"), p, useDingbats = FALSE)

    pdf(file.path(out_dir, "eigenHMAP.pdf"), width = 10, useDingbats = FALSE)
    getEigenHmap(mat)
    dev.off()
}

if (plain == TRUE) quit()

mat <- as.data.frame(mat)

# make hexbin-bar-density grid plot (all vs all)
p <- getEigenPairs(mat)
ggsave(file.path(out_dir, "eigenPairsGrid.pdf"), p, height = 14, width =16, useDingbats=FALSE)

# get separate hexbin and density plots (all vs all)
hexbin_dir <- file.path(out_dir, "hexbin_plots")
density_dir <- file.path(out_dir,"density_plots")
if (!dir.exists(hexbin_dir)) dir.create(hexbin_dir) else warning("hexbin_plots directory already exists, it will be overwritten!", call. = FALSE)
if (!dir.exists(density_dir)) dir.create(density_dir) else warning("density_plots directory already exists, it will be overwritten!", call. = FALSE)

for (k in 1:nrow(comparisons)) {
    x <- sym(comparisons[k, 1])
    y <- sym(comparisons[k, 2])
    suffix <- paste0(comparisons[k, 2], "_vs_", comparisons[k, 1], ".pdf")

    p <- getHexPlot(mat, aes(!!x, !!y), grid=FALSE)
    ggsave(file.path(hexbin_dir, paste0("eigenHexbinPlot_", suffix)), useDingbats = FALSE)

    p <- getDensPlot(mat, aes(!!x, !!y), grid=FALSE)
    ggsave(file.path(density_dir, paste0("eigenDensityPlot", suffix)), useDingbats = FALSE)
}

# get all barplots in one plot
p <- getBarPlot(df.sizes, NULL, grid = FALSE)
ggsave(file.path(out_dir, "size_barplots.pdf"), p, width = 9,height = 5, useDingbats=FALSE)

#plot bars for pc1 classified bins
p <- ggplot(df.switch, aes(size, comparisonID)) +
    geom_col(aes(fill = switch), position = "stack", width = .6, size = .5) +
    scale_y_discrete(labels = function(x) sub(":", " vs ", x)) +
    scale_fill_manual(values=c("blue","darkred","darkblue","red"))+
    labs(x = "size (Mb)", y = "",  fill = "", title = "Compartment Switches",
         subtitle = sprintf("Criteria: PC1 sign flips & abs(relative delta cutoff) is > %g", delta.cut)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size=9),
          plot.subtitle = element_text(hjust = 0.5, size=9))
ggsave(file.path(out_dir, "switch_bins_barplot.pdf"), p, width = 12, useDingbats = FALSE)

# write df.switch
df.switch[, switch := paste0(switch, ".number")]
df.switch <- dcast(df.switch, comparisonID ~ switch, value.var = "N")
df.switch[, c("sampleName1", "sampleName2") := tstrsplit(comparisonID, split = ":", fixed = T)]
setcolorder(df.switch, c("sampleName1", "sampleName2", "AA.number", "AB.number", "BB.number", "BA.number", "comparisonID"))
fwrite(df.switch, file.path(out_dir, "pc1.switch_summary.tsv"), sep="\t")

# write every comparison individually
for (k in 1:nrow(comparisons)) {
    DT <- df.comparisons[[k]]
    DT[, coord := paste(chr, start, end, sep = ":")]
    DT[, c("chr", "start", "end", "size") := NULL]
    setcolorder(DT, "coord")
    setnames(DT, c("X", "Y"), comparisons[k, ])
    fname <- sprintf("pc1.switch_%s_vs_%s.tsv", comparisons[k, 1], comparisons[k, 2])
    fwrite(DT, file.path(out_dir, fname), sep="\t")
}

# Plot Metrics
if (file.exists(metrics_file)){
    metrics <- c("size.per.count.HK","size.per.count.TSS", "size.bp")
    DT <- fread(metrics_file)[chromosome !=  "chrY" & signGroup == "log(ratio[1:-1])"]
    p <- getMetricsPlot(DT, metrics)
    ggsave(file.path(out_dir, "metrics_logratio.AB.pdf"), p,
           width = 14, height = 7, useDingbats=FALSE)
}

