library(ggplot2)
library(reshape2)
library(scales)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=F)

colnames(x)[1] = c("cnv")
x = melt(x, id.vars=c("cnv"))
x$cn = round(x$value)
if (sum(x$cn > 9)) { x[x$cn > 9,]$cn = 9; }

x$cn = factor(x$cn, levels=0:9)
nsamples = length(unique(x$variable))
nbins = sqrt(nsamples)
if (nbins < 30) { nbins = 30; }

# Plot CNVs
for (CNV in unique(x$cnv)) {
  print(CNV)
  df = x[x$cnv == CNV,]
  p = ggplot(data=df, aes(x=value))
  for(i in 0:9) {
    p = p + geom_histogram(data=subset(df, cn == i), aes(fill=cn), bins=nbins)
  }
  p = p + xlab("Copy-number")
  p = p + ylab("Count")
  p = p + scale_x_continuous(breaks=0:10, labels=comma)
#  p = p + scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"), drop=F)
  p = p + scale_fill_manual(values=c("#ff7f00", "#1f78b4","#33a02c","#e31a1c","#6a3d9a", "#fdbf6f", "#a6cee3", "#b2df8a", "#fb9a99", "#cab2d6"), drop=F)
  p = p + ggtitle(CNV)
  ggsave(p, file=paste0(CNV, ".png"), width=24, height=6)
}
