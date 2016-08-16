library(ggplot2)
library(reshape2)
library(scales)
library(optparse)
library(VariantAnnotation)

gg_color_hue = function(n) {
  hues=seq(15,375,length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


optl = list(
     	    make_option(c("-f", "--file"), type="character", default=NULL, help="input coverage file"),
	    make_option(c("-r", "--region"), type="character", default=NULL, help="region chr:start-end"),
	    make_option(c("-v", "--vcf"), type="character", default=NULL, help="SV vcf file"),
	    make_option(c("-s", "--sample"), type="character", default=NULL, help="sample ID")
	    )

opt_parser = OptionParser(option_list=optl)
opt = parse_args(opt_parser)

if (is.null(opt$file)) {
   print_help(opt_parser)
   stop("Input coverage file is missing!", call.=FALSE)
}

  

cov = read.table(opt$file, header=T)
col = 5
if (!is.null(opt$sample) && (opt$sample %in% colnames(cov))) {  col = which(colnames(cov) == opt$sample); }
chr = NULL
if (!is.null(opt$region)) {
   s=unlist(strsplit(opt$region, ":"))
   chr=s[1]
   minStart = 0
   maxStart = max(cov[cov$chr==chr,]$start) + 1
   if (length(s) > 1) {
      v=unlist(strsplit(s[2], "-"))
      minStart = as.integer(v[1])
      maxStart = as.integer(v[2])
   }
}
sv = NULL
if (!is.null(opt$vcf)) {
  vcf = readVcf(opt$vcf, "hg19")
  svChr = as.vector(seqnames(rowRanges(vcf)))
  svStart = start(rowRanges(vcf))
  svEnd = info(vcf)$END
  svId = names(rowRanges(vcf))
  svType = info(vcf)$SVTYPE
  sv = data.frame(svChr, svStart, svEnd, svId, svType)
  sv[sv$svType=="TRA",]$svEnd = sv[sv$svType=="TRA",]$svStart
  sv$middle = as.integer( (sv$svStart + sv$svEnd) / 2 )
  sv$yv = 4
  sv = melt(sv, id.vars=c("svChr", "svType", "svId", "yv"))
  sv[sv$variable=="middle",]$yv = sv[sv$variable=="middle",]$yv + 1.5
}


  

cov=cov[cov[,col]!=0,]
md=median(cov[,col])
cov$cn=(cov[,col]/md)*2
lab=colnames(cov)[col]

if (is.null(chr)) { chrs = unique(cov$chr); } else { chrs = c(chr); }

# Plot individual chromosomes
for (cr in chrs) {
    chrCov = cov[cov$chr==cr,]
    if (!is.null(chr)) chrCov = chrCov[(chrCov$start>=minStart) & (chrCov$start < maxStart),]
    maxcn = ceiling(quantile(chrCov$cn, 0.999))
    if (maxcn < 6) maxcn = 6
    if (max(chrCov$cn) > maxcn) { chrCov[chrCov$cn>maxcn,]$cn = maxcn; }
    
    p1=ggplot(data=chrCov, aes(x=start, y=cn)) + geom_point(alpha=1/4, size=0.5)
    p1=p1 + xlab(paste0(lab, " - ", cr)) + ylab("CN estimate")
    p1=p1 + scale_y_continuous(breaks=0:maxcn, labels=as.character(0:maxcn), limits=c(0,maxcn))
    p1=p1 + scale_x_continuous(labels=comma)
    if (!is.null(sv)) {
       svChr = sv[sv$svChr==cr,]
       if (nrow(svChr)) {
              p1=p1 + geom_line(data=svChr, aes(x=value, y=yv, group=svId, colour=svType));
       }
    }
    ggsave(paste0(lab, ".", cr, ".png"), width=16, height=6)
}
print(warnings())
quit()
