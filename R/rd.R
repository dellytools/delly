library(ggplot2)
library(reshape2)
library(scales)

gg_color_hue = function(n) {
  hues=seq(15,375,length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

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
md=median(cov[,5])
cov$cn=(cov[,5]/md)*2
lab=colnames(cov)[5]

if (is.na(chr)) { chrs = unique(cov$chr); } else { chrs = c(chr); }

# Plot individual chromosomes
for (cr in chrs) {
    chrCov = cov[cov$chr==cr,]
    if (!is.na(chr)) chrCov = chrCov[(chrCov$start>=minStart) & (chrCov$start < maxStart),]
    if (max(chrCov$cn) > 6) { chrCov[chrCov$cn>6,]$cn = 6; }
    
    p1=ggplot(data=chrCov, aes(x=start, y=cn)) + geom_point(alpha=1/4, size=0.5)
    p1=p1 + xlab(paste0(lab, " - ", cr)) + ylab("CN estimate")
    p1=p1 + scale_y_continuous(breaks=c(0,1,2,3,4,5,6), labels=c("0", "1", "2", "3", "4", "5", ">6"), limits=c(0,6))
    p1=p1 + scale_x_continuous(labels=comma)
    ggsave(paste0(lab, ".", cr, ".png"), width=16, height=6)
}
print(warnings())
quit()
