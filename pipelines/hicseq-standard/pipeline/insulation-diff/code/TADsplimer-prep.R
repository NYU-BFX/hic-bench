#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

mat.internal = args[1]
inpdir.internal = args[2]
object.internal = args[3]
outdir.internal = args[4]
out.prefix.internal = args[5]

matrix.to.start = read.table(paste0(inpdir.internal, "/", mat.internal), sep="\t", header=T)
matrix.to.start = as.matrix(matrix.to.start[, -1])

if (!dir.exists(paste0(outdir.internal, "/", object.internal))) {
  dir.create(paste0(outdir.internal, "/", object.internal))
}
fn.out = paste0(outdir.internal, "/", object.internal, "/", out.prefix.internal, "_", mat.internal)
write.table(matrix.to.start, file=fn.out, sep="\t", quote=F, row.names=F, col.names=F)
