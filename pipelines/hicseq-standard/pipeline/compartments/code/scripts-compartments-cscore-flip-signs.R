### USAGE ###
# EXAMPLE:
# Rscript ./code/scripts-compartments-cscore-flip-signs.R $outdir/__jdata/each.${chr}/${chr}_cscore.bedgraph /gpfs/data/tsirigoslab/hicbench-repository/data-repo/genomes/mm10/HK_genes.bed

args = commandArgs(trailingOnly=T)
bedgraph.input = args[1]
hkgene.input = args[2]
# samplePath = argv[1L] #bedgraph file where to apply the chromosome pc1 sign inversion/s
# chromsList = argv[2L] #comma separeted list of the chromosomes to be flipped

# message(bedgraph.input)
# message(file.exists(bedgraph.input))
# message(getwd())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))

chr.bedgraph.header = fread(bedgraph.input, sep="\n", nrows=1, header=F)

# For chromosomes X, Y, M, there is a possibility that CscoreTools did not produce *any* output.
# There is a possibility that chr.bedgraph has no data, and the downstream product cannot be generated.

# Error handling
chr.bedgraph.nr = fread(cmd=paste0("wc -l ", bedgraph.input))$V1
if (chr.bedgraph.nr == 1) {
  message("Because of lack of data, no part of sign-changing script was run. This happened probably because it was chrX, Y, or M.")
} else {
  chr.bedgraph = fread(bedgraph.input, skip=1, header=F)
  chr.bedgraph.gr = GRanges(chr.bedgraph$V1, IRanges(chr.bedgraph$V2 + 1, chr.bedgraph$V3), strand="*")
  chr.bedgraph.gr$score = chr.bedgraph$V4

  hkgene.gr = import.bed(hkgene.input)

  # print(summary(chr.bedgraph.gr$score))

  chr.subsetbyoverlaps.bedgraph.gr = subsetByOverlaps(chr.bedgraph.gr, hkgene.gr)

  # print(summary(chr.subsetbyoverlaps.bedgraph.gr$score))
  subsetbyoverlaps.gr.length = length(chr.subsetbyoverlaps.bedgraph.gr)

  if (subsetbyoverlaps.gr.length == 0) {
    message("Intersection with housekeeping gene resulted in zero regions. No use of KS test here.")
  } else {
    ks.test.result.correct = ks.test(chr.bedgraph.gr$score, chr.subsetbyoverlaps.bedgraph.gr$score, alternative="greater")
    ks.test.result.reverse = ks.test(chr.bedgraph.gr$score, chr.subsetbyoverlaps.bedgraph.gr$score, alternative="less")
    pdf(paste0(bedgraph.input, ".quality.pdf"), width=4, height=4)

    # First. Black is all
    plot.ecdf(chr.bedgraph.gr$score, col="black", lwd=1.5, verticals=T, pch="", col.01line=F, main=basename(bedgraph.input), xlab="Raw Cscore", ylab="Cumulative distribution")
    # Second. Red is enriched for HK genes
    plot.ecdf(chr.subsetbyoverlaps.bedgraph.gr$score, col="red", lwd=1.5, verticals=T, pch="", col.01line=F, add=T)

    legend("topleft", c("All", "Enriched for HK genes", paste0("KS test D: ", round(ks.test.result.correct$statistic, 2), ", p: ", sprintf("%.3f", ks.test.result.correct$p.value))), fill=c("black", "red", NA), bty="n")

    dev.off()

    # The condition defined for "needs flipping": "correct" is false and "reverse" is true.
    # D cutoff is +0.20 and p cutoff is 0.01 (pretty lenient).
    needs.to.be.flipped = (!(ks.test.result.correct$statistic >= +0.20 & ks.test.result.correct$p.value <= 0.01)) & (ks.test.result.reverse$statistic >= +0.20 & ks.test.result.reverse$p.value <= 0.01)

    # Condition for ambiguous matching: hope this is very rare, and that's what I observe by far.
    # One: both not significant
    ambiguous.condition.one = (ks.test.result.correct$statistic < +0.20 & ks.test.result.correct$p.value > 0.01) & (ks.test.result.reverse$statistic < +0.20 & ks.test.result.reverse$p.value > 0.01)
    if (ambiguous.condition.one) {
      print("The KS test did not produce significant difference on either direction. Please check if the enrichment is actually happening!")
    }

    # Two: both significant widely.
    ambiguous.condition.two = (ks.test.result.correct$statistic >= +0.05 & ks.test.result.correct$p.value <= 0.01) & (ks.test.result.reverse$statistic >= +0.05 & ks.test.result.reverse$p.value <= 0.01)
    if (ambiguous.condition.two) {
      print("The KS test on either direction is producing significant difference. This means KS test is inadequate to determine the signs. Please use your best judgment whether to flip the sign or not.")
    }


    # If statement works only if flipping is needed.
    if (needs.to.be.flipped) {
      # Rename the old file
      file.rename(bedgraph.input, paste0(bedgraph.input, ".old"))

      # times -1, score
      chr.bedgraph.gr$score = -1 * chr.bedgraph.gr$score
      chr.subsetbyoverlaps.bedgraph.gr$score = -1 * chr.subsetbyoverlaps.bedgraph.gr$score

      # draw cdf plot once again
      ks.test.result.correct = ks.test(chr.bedgraph.gr$score, chr.subsetbyoverlaps.bedgraph.gr$score, alternative="greater")

      pdf(paste0(bedgraph.input, ".quality.Sign.flipped.pdf"), width=4, height=4)

      plot.ecdf(chr.bedgraph.gr$score, col="black", lwd=1.5, verticals=T, pch="", col.01line=F, main=basename(bedgraph.input), xlab="Raw Cscore", ylab="Cumulative distribution")
      plot.ecdf(chr.subsetbyoverlaps.bedgraph.gr$score, col="red", lwd=1.5, verticals=T, pch="", col.01line=F, add=T)

      legend("topleft", c("All", "Enriched for HK genes", paste0("KS test D: ", round(ks.test.result.correct$statistic, 2), ", p: ", sprintf("%.3f", ks.test.result.correct$p.value))), fill=c("black", "red", NA), bty="n")

      dev.off()

      options(scipen=999)
      # Write cscore bedgraph file
      fwrite(chr.bedgraph.header, file=bedgraph.input, quote=F, col.names=F)
      fwrite(data.frame(seqnames(chr.bedgraph.gr), start(chr.bedgraph.gr) - 1, end(chr.bedgraph.gr), chr.bedgraph.gr$score), file=bedgraph.input, append=T, sep="\t", quote=F, col.names=F, row.names=F)

    }

  }
}
