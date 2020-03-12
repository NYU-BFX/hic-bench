#!/usr/bin/env Rscript
# Libraries ---------------------------------------------------------------

library(argparser)
library(data.table)
library(RColorBrewer)

norm_score <- function(x) x / sum(x, na.rm = TRUE) * 1e6
rotate <- function(x) t(apply(x, 2L, rev))

# Other Defaults
cols = colorRampPalette(colors = c("darkblue", "blue4", "blue3", "white", "red3", "red4", "darkred"), space = "Lab")(100)
cols_red = colorRampPalette(colors = brewer.pal(n = 9L, name = 'YlOrRd'), space = "Lab")(100)


# Parse Arguments ---------------------------------------------------------

p <- arg_parser("R-script to generate 'cute' contact matrices")
p <- add_argument(p, "mat1", help = "First matrix")
p <- add_argument(p, "mat2", help = "Second matrix")
p <- add_argument(p, "left", type = "integer", help = "Starting position")
p <- add_argument(p, "right", type = "integer", help = "Ending position")
p <- add_argument(p, "sample1", help = "Name for 1st sample")
p <- add_argument(p, "sample2", help = "Name for 2nd sample")
p <- add_argument(p, "output", help = "Output prefix")
p <- add_argument(p, "--normalize", help = "Normalize the matrices", flag = TRUE)
p <- add_argument(p, "--binsize", default = 1e4L, help = "Bin Size")
p <- add_argument(p, "--pseudo", default = 0.1, help = "Pseudo count to add before log")
p <- add_argument(p, "--max_logFC", default = 2, help = "Maximum value for logFC to color")
p <- add_argument(p, "--height", default = 2048L, help = "Image height in px")
p <- add_argument(p, "--width", default = 2048L, help = "Image width in px")
p <- add_argument(p, "--res", default = 300, help = "Image resolution")

argv <- parse_args(p)
#argv <- parse_args(p, strsplit("--normalize /gpfs/home/sn1110/carrolllab/0719_diagnosed_relapsed/1019_hic_5pairs_shallow/__05a-matrix-ic/results/matrix-ic.by_sample.cutoff_0/matrix-filtered.by_sample.res_20kb/filter.by_sample.mapq_20/align.by_sample.bowtie2/D_D_S7/matrix.chr4.tsv  /gpfs/home/sn1110/carrolllab/0719_diagnosed_relapsed/1019_hic_5pairs_shallow/__05a-matrix-ic/results/matrix-ic.by_sample.cutoff_0/matrix-filtered.by_sample.res_20kb/filter.by_sample.mapq_20/align.by_sample.bowtie2/D_R_S8/matrix.chr4.tsv  73400000  75800000  Diagnosis  Relapse  cxcl8_73.5-75.8mb", "\\s+")[[1]])

mat1_file <- argv$mat1
mat2_file <- argv$mat2
pos_left <- argv$left
pos_right <- argv$right
sample_name_1 <- argv$sample1
sample_name_2 <- argv$sample2
out_prefix <- argv$output
pseudo_count <- argv$pseudo
normalize <- argv$normalize  # WARNING: Default behavior changed
bin.size <- argv$binsize
max_abs <- argv$max_logFC
fig_height = argv$height
fig_width = argv$width
fig_res = argv$res


# Main --------------------------------------------------------------------

i <- pos_left / bin.size
j <- pos_right / bin.size

# Read Inputs
mat1 <- as.matrix(fread(mat1_file),  rownames = 1L)
mat2 <- as.matrix(fread(mat2_file),  rownames = 1L)

if (normalize) {
    mat1 <- norm_score(mat1)
    mat2 <- norm_score(mat2)
}

mat1 = mat1[i:j, i:j, drop = FALSE]
mat2 = mat2[i:j, i:j, drop = FALSE]
mixed_matrix <- mat1
mixed_matrix[lower.tri(mixed_matrix)] <- mat2[lower.tri(mat2)]

# Rotate matrix for plotting
mat1 = rotate(mat1)
mat2 = rotate(mat2)
mixed_matrix = rotate(mixed_matrix)

# Plot logFC matrices
# WARNING: In the original script the log2 was taken over log10. Put at the end if you want to retain behavior
logFC <- log2(mat2 + pseudo_count) - log2(mat1 + pseudo_count)
max_abs <- pmin(max_abs, max(abs(logFC), na.rm = TRUE))
png(paste0(out_prefix, "_logFC_", sample_name_1, "-vs-", sample_name_2, ".png"),
    height = fig_height, width = fig_width, res = fig_res)
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(logFC, col = cols, zlim = c(-max_abs, max_abs), xaxs = "i", axes = FALSE)
dev.off()

# Plot log10 of 1st matrix
mat1 <- log10(mat1 + pseudo_count)
png(paste0(out_prefix, "_", sample_name_1, ".png"),
    height = fig_height, width = fig_width, res = fig_res)
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(mat1,col = cols_red, xaxs = "i", axes = FALSE)
box(lwd = 10)
dev.off()

# Plot log10 of 2nd matrix
mat2 <- log10(mat2 + pseudo_count)
png(paste0(out_prefix,  "_", sample_name_2, ".png"),
    height = fig_height, width = fig_width, res = fig_res)
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(mat2, col = cols_red, xaxs = "i", axes = FALSE)
box(lwd = 10)
dev.off()

# Merge 2 matrices (1st on top, 2nd at the bottom)
mixed_matrix <- log10(mixed_matrix + pseudo_count)
png(paste0(out_prefix, "_", sample_name_1, "-top_", sample_name_2, "-bottom.png"),
    height = fig_height, width = fig_width, res = fig_res, pointsize = 45)
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(mixed_matrix, col = cols_red, xaxs = "i", axes = FALSE)
box(lwd = 10)
dev.off()

