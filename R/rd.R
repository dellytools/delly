library(DNAcopy)
library(ggplot2)
library(scales)
library(gtable)
library(grid)

chrNamesLong = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22", "chrX")
#chrNamesLong = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22", "chrX", "chrY")
chrNamesShort = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
#chrNamesShort = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

# Params
minCN = 0
sdUNDO = 1.5

# Parse coverage table
args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)

# Fix chromosome ordering
if (sum(x$chr %in% chrNamesLong) > sum(x$chr %in% chrNamesShort)) { chrs = chrNamesLong; } else { chrs = chrNamesShort; }
x = x[x$chr %in% chrs,]
x$chr = factor(x$chr, levels=chrs)

# Log-ratio values
x$pos = (x$start + x$end) / 2
baseCN = median(x[,6])
x$logratio = log2(x[,6] / baseCN)

# Segmentation
if (length(args)>1) {
 # Provided on the command-line
 seg = read.table(args[2], header=F, sep="\t")
 colnames(seg) = c("chr", "start", "end", "id", "cn")
} else {
 # Segment logR
 cnaData = CNA(x$logratio, maploc=x$pos, chrom=x$chr, sampleid="Sample", presorted=T)
 cnaSegments = segment(smooth.CNA(cnaData), undo.splits="sdundo", undo.SD=sdUNDO)
 seg = data.frame(segments.summary(cnaSegments))
 colnames(seg) = c("ID", "chr", "start", "end", "num.mark", "seg.mean", "seg.sd", "seg.median", "seg.mad")
 seg$cn = 2^seg$seg.median * baseCN
 seg$chr = factor(seg$chr, levels=chrs)
}

# Whole genome
maxCN = as.integer(max(x[,6])+1)
#maxCN = 8
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
 maxCN = as.integer(max(sub[,6])+1)
 p = ggplot(data=sub, aes(x=start, y=sub[,6]))
 p = p + geom_point(pch=21, color="black", fill="black", size=0.5)
 p = p + ylab("Copy-number") + xlab(chrname)
 p = p + scale_x_continuous(labels=comma, breaks = scales::pretty_breaks(n=20))
 p = p + scale_y_continuous(labels=comma, breaks = c(minCN:maxCN), limits=c(minCN, maxCN))
 if (nrow(sl)) { p = p + geom_segment(data=sl, aes(x=start, y=cn, xend=end, yend=cn), color="#31a354", size=1.2); }
 p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
 p = p + ggtitle(args[1])
 ggsave(p, file=paste0("plot.", chrname, ".png"), width=24, height=6)
 print(warnings())
}
