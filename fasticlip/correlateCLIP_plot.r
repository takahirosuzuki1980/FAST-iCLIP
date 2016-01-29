args=commandArgs(TRUE)
fn=args[1]

print("Reading in file")
data=read.delim(fn, header=FALSE, nrow=5)
nc = ncol(data)
data=read.delim(fn, header=FALSE, colClasses=c(rep("character",nc-2),rep("numeric",2)))
x1 = log10(data[,(nc-1)]+1)
x2 = log10(data[,nc]+1)
r.pear=cor(x1, x2)
r.spe=cor(x1, x2, method="spe")
print(paste("Pearson correlation:", r.pear, sep=' '))
print(paste("Spearman correlation:", r.spe, sep = ' '))

png("plot.png")
plot(x1, x2, pch=19, cex=0.3, xlab="log10 RT stops in Dataset 1", ylab="log10 RT stops in Dataset 2", main="RT stop scatterplot")
dev.off()