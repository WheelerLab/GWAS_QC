args <- commandArgs(trailingOnly=T)
"%&%" <- function(a, b) paste(a, b, sep="")

header <- args[1]
pcsfile <- header %&% '.evec'
#popfile <- header %&% '.pop'

pcs <- read.table(pcsfile)
#pop <- read.table(popfile)
#colnames(pop) <- 'pop'
#new <- cbind(pcs,pop)

pdf(file=header %&% '.PCA.plot.pdf')
plot(pcs$V2,pcs$V3)
#legend('bottomright',legend=unique(new$pop),col=unique(new$pop),pch=1)
dev.off()

