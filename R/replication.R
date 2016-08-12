library(ggplot2)
library(reshape2)
library(scales)

# Read data
args=commandArgs(trailingOnly=TRUE)
cov = read.table(args[1], header=T)
chr=NA
if (length(args) > 1) {
   s=unlist(strsplit(args[2], ":"))
   chr=s[1]
   minStart = 0
   maxStart = max(cov[cov$chr==chr,]$start) + 1
   if (length(s) > 1) {
      v=unlist(strsplit(s[2], "-"))
      minStart = as.integer(v[1])
      maxStart = as.integer(v[2])
   }
}
cov=cov[cov[,5]!=0,]
s=boxplot(cov[,5], plot=F)$stats
lq=s[2,1]
uq=s[4,1]
iqrcov=cov[cov[,5]>=lq & cov[,5]<=uq,5]
mu=mean(iqrcov)
sdev=sd(cov[,5])
cov=cov[cov[,5]>=(mu - 3*sdev) & cov[,5]<=(mu+3*sdev),]
cov$zscores=(cov[,5] - mu)/sdev
lab=colnames(cov)[5]

if (is.na(chr)) { chrs = unique(cov$chr); } else { chrs = c(chr); }

# Plot individual chromosomes
plot = T
for (cr in chrs) {
    chrCov = cov[cov$chr==cr,]
    if (!is.na(chr)) chrCov = chrCov[(chrCov$start>=minStart) & (chrCov$start < maxStart),]
    chrCov$spl = fitted(with(chrCov, smooth.spline(start, zscores, spar=0.3)))

    # Find contiguous blocks
    medDist = median(chrCov$start[-1] - chrCov$start[1:nrow(chrCov)-1])
    runs = rle((chrCov$start[-1] - chrCov$start[1:nrow(chrCov)-1]) == medDist)
    idx = 1
    rd = data.frame()
    for(i in 1:length(runs$values)) {
    	  if ((runs$values[i]==T) && (runs$length[i]>20)) {
	     subset = chrCov[idx:(idx+runs$length[i]),]
	     subset$grp = idx
	     rd = rbind(rd, subset)
	  }
	  idx = idx + runs$length[i]
    }

    # Plot segments
    if (plot) {
        p1=ggplot(data=rd, aes(x=start, y=zscores)) + geom_point(alpha=1/4, size=1)
    	p1=p1 + xlab(paste0(lab, " - ", cr)) + ylab("Read depth")
	p1=p1 + scale_y_continuous(labels=comma) 
	p1=p1 + scale_x_continuous(labels=comma)
	p1=p1 + geom_line(data=rd, aes(x=start, y=spl, group=grp), color="darkblue", size=1.5)
	ggsave(paste0(lab, ".", cr, ".png"), width=16, height=6)
    } else {
      	colnames(rd)[colnames(rd)=="spl"] = paste0(lab, "spl");
        write.table(rd[,c("chr", "start", paste0(lab, "spl"))], file=paste0(lab, ".", cr, ".spl.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    }
}
print(warnings())
quit()

