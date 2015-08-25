args <- commandArgs(trailingOnly=T)
"%&%" <- function(a, b) paste(a, b, sep="")

header <- args[1]
pcsfile <- header %&% '.evec'
popfile <- header %&% '.pop'

pcs <- read.table(pcsfile)
pop <- read.table(popfile)
colnames(pop) <- 'pop'
new <- cbind(pcs,pop)

pdf(file=header %&% '.PCA.plot.pdf')
plot(new$V2,new$V3,col=new$pop)
legend('bottomright',legend=unique(new$pop),col=unique(new$pop),pch=1)
dev.off()

###choose individuals from GWAS	for homogeneous	GWAS###

ceu <- new[new$pop=='CEU',]
gwas <- new[new$pop == 'GWAS',]

uV2 <- mean(ceu$V2) + 20*sd(ceu$V2)
lV2 <- mean(ceu$V2) - 20*sd(ceu$V2)

uV3 <- mean(ceu$V3) + 20*sd(ceu$V3)
lV3 <- mean(ceu$V3) - 20*sd(ceu$V3)

inclusion <- gwas[gwas$V2 >= lV2,]
inclusion <- inclusion[inclusion$V2 <= uV2,]
inclusion <- inclusion[inclusion$V3 >= lV3,]
inclusion <- inclusion[inclusion$V3 <= uV3,]

samples <- inclusion[,1]

write.table(samples,file=header %&% ".euro.GWAS.PCAs",quote=F,row.names=F,col.names=F)
