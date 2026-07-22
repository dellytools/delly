library(ggplot2)
library(reshape2)
library(scales)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=F, na.strings=c(".", "NA"))

colnames(x)[1] = c("cnv")
x = melt(x, id.vars=c("cnv"))
x = x[!is.na(x$value) & x$value != -1,]   # drop missing / invalid RDCN (-1 sentinel)
x$cn = round(x$value)
if (sum(x$cn > 9)) { x[x$cn > 9,]$cn = 9; }

x$cn = factor(x$cn, levels=0:9)
nsamples = length(unique(x$variable))
nbins = 2 * ceiling(sqrt(nsamples))
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
  p = p + scale_x_continuous(breaks=0:10, labels=comma, limits=c(0,6))
  p = p + scale_fill_manual(values=c("#ff7f00", "#1f78b4","#33a02c","#e31a1c","#6a3d9a", "#fdbf6f", "#a6cee3", "#b2df8a", "#fb9a99", "#cab2d6"), drop=F)
  p = p + ggtitle(CNV)
  ggsave(p, file=paste0(CNV, ".png"), width=24, height=6)
}
