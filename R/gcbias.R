library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
x = read.table("gc.table", header=T)
x$gc = x$gcsum / (nrow(x)-1)
x$fractionSample = x$fractionSample * 100
x$fractionReference = x$fractionReference * 100
df = melt(x[,c("gc","fractionSample","fractionReference")], id.vars=c("gc"))

# Whole genome
p = ggplot(data=df, aes(x=gc, y=value))
p = p + geom_bar(aes(color=variable, fill=variable), stat="identity")
p = p + xlab("GC content")
p = p + ylab("Obs / Exp")
p = p + ylim(0, max(max(x$obsexp), max(x$fractionSample + x$fractionReference)))
p = p + geom_line(data=x, aes(x=gc, y=obsexp), color="black")
ggsave(p, file="gcbias.png", width=12, height=6)
print(warnings())
