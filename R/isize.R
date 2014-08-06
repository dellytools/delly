library(ggplot2)
library(grid)

args=commandArgs(trailingOnly=TRUE)
orient=c("F+", "F-", "R+", "R-")
ins=read.table(args[1], header=T)
count=c()
titleVal="Insert Size Histogram ("
for (o in orient) {
    s=sum(ins[ins$Orient==o,]$Count)
    count=c(count, s)
    lab=paste("#", o, ":", s, sep="")
    titleVal=paste(titleVal, lab, sep="  ")
}
titleVal=paste(titleVal, ")", sep=" ")
countTotal=sum(count)
	
png(paste(args[1], "png", sep="."), height=800, width=800)
grid.newpage()
pushViewport(viewport(layout=grid.layout(4,2)))
p1=ggplot(data=ins) + geom_histogram(aes(x=Size, y=Count, colour=Orientation, fill=Orientation), stat="identity") + ggtitle(titleVal)
p2=ggplot(data=ins) + geom_histogram(aes(x=Size, y=Count, colour=Orientation, fill=Orientation), stat="identity") + facet_wrap(~ Orientation, nrow=2, ncol=2)
print(p1, vp = viewport(layout.pos.row=1:2, layout.pos.col=1:2))
print(p2, vp = viewport(layout.pos.row=3:4, layout.pos.col=1:2))
z=dev.off()
