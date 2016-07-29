library(ggplot2)
library(reshape2)
library(scales)

# Read data
args=commandArgs(trailingOnly=TRUE)
cov = read.table(args[1], header=T)
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


# Plot individual chromosomes
for (cr in unique(cov$chr)) {
    chrCov=cov[cov$chr==cr,]
    p1=ggplot(data=chrCov, aes(x=start, y=zscores)) + geom_point(alpha=1/4, size=1)
    p1=p1 + stat_smooth() + scale_y_continuous(labels=comma) 
    p1=p1 + xlab(paste0(lab, " - ", cr)) + ylab("Read depth")
    p1=p1 + scale_x_continuous(labels=comma)
    ggsave(paste0(lab, ".", cr, ".png"), width=16, height=6)
}
print(warnings())
quit()

