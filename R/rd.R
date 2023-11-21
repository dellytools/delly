library(ggplot2)
library(scales)
library(gtable)
library(grid)

chrNamesLong = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22", "chrX")
chrNamesShort = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
minCN = 0
maxCN = 8
seg = data.frame()
if (length(args)>1) {
   seg = read.table(args[2], header=F, sep="\t")
   colnames(seg) = c("chr", "start", "end", "id", "cn")
}

# Fix chromosome ordering
if (sum(x$chr %in% chrNamesLong) > sum(x$chr %in% chrNamesShort)) { chrs = chrNamesLong; } else { chrs = chrNamesShort; }
x = x[x$chr %in% chrs,]
x$chr = factor(x$chr, levels=chrs)
if (nrow(seg) > 0) {
 seg = seg[seg$chr %in% chrs,]
 seg$chr = factor(seg$chr, levels=chrs)
}

# Whole genome
p = ggplot(data=x, aes(x=start, y=x[,6]))
p = p + geom_point(pch=21, color="black", fill="black", size=0.5)
p = p + xlab("Chromosome")
p = p + ylab("Copy-number")
p = p + scale_x_continuous(labels=comma)
if (nrow(seg)) { p = p + geom_segment(data=seg, aes(x=start, y=cn, xend=end, yend=cn), color="#31a354", size=1.2); }
p = p + facet_grid(. ~ chr, scales="free_x", space="free_x")
p = p + scale_y_continuous(labels=comma, breaks = c(minCN:maxCN), limits=c(minCN, maxCN))
p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
p = p + ggtitle(args[1])
ggsave(p, file="plot.wholegenome.png", width=24, height=6)
print(warnings())

# By chromosome
for(chrname in unique(x$chr)) {
 print(chrname)
 sub = x[x$chr == chrname,]
 sl = seg[seg$chr == chrname,]
 p = ggplot(data=sub, aes(x=start, y=sub[,6]))
 p = p + geom_point(pch=21, color="black", fill="black", size=0.5)
 p = p + ylab("Copy-number") + xlab(chrname)
 p = p + scale_x_continuous(labels=comma, breaks = scales::pretty_breaks(n=20))
 p = p + scale_y_continuous(labels=comma, breaks = c(minCN:maxCN), limits=c(minCN, maxCN))
 if (nrow(sl)) { p = p + geom_segment(data=sl, aes(x=start, y=cn, xend=end, yend=cn), color="#31a354", size=1.2); }
 p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
 ggsave(p, file=paste0("plot.", chrname, ".png"), width=24, height=6)
 print(warnings())
}

