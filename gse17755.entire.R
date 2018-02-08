# source("http://bioconductor.org/biocLite.R")
# biocLite(c("GO.db","GOstats","Rgraphviz", "hgu95av2.db"))

#load("E:/work/UCB/project/gene expression data/GSE17755/R data/current.RData")
# load("current.RData")
# save.image("current.RData")
rm(list=ls())
options(stringsAsFactors=FALSE) 

# library("GEOquery")
library("limma")
# library("Biobase")
# library("affy")
#library("impute")
#library("arrayQuality")
#library("arrayQualityMetrics") 

setwd("E:/UCB/project/gene expression data/GSE17755/gse17755")
ds.file.list <- list.files()

#The row names and column names of the data matrix
ds.cnames <- unlist(lapply(ds.file.list, function(x) strsplit(x, split="_")[[1]][[1]]))
ds.cnames
ds.gsm <- read.table(ds.file.list[1], sep="\t", fill=T, header=T, quote="")
ds.rnames <- ds.gsm[,2]

#read in the disease patient samples
for(i in 1:length(ds.file.list)){
  ds.gsm <- read.table(ds.file.list[i], sep="\t", fill=T, header=T, quote="")
	if(i==1){
	  ds.red.df <- matrix(NA, nrow=nrow(ds.gsm), ncol=length(ds.file.list), dimnames=list(ds.rnames, ds.cnames))
		ds.green.df <- matrix(NA, nrow=nrow(ds.gsm), ncol=length(ds.file.list), dimnames=list(ds.rnames, ds.cnames))
		ds.red.bg.df <- matrix(NA, nrow=nrow(ds.gsm), ncol=length(ds.file.list), dimnames=list(ds.rnames, ds.cnames))
		ds.green.bg.df <- matrix(NA, nrow=nrow(ds.gsm), ncol=length(ds.file.list), dimnames=list(ds.rnames, ds.cnames))
	}
  
  if (length(which(ds.gsm[,2]==ds.rnames))!=length(ds.rnames)){
    print("probe names are differnet")
  } 
  ds.red.df[, i] <- ds.gsm[, 25]
  ds.green.df[, i] <- ds.gsm[, 26]
  ds.red.bg.df[, i] <- ds.gsm[, 27]
  ds.green.bg.df[, i] <- ds.gsm[, 28]
	print(i)
}

setwd("E:/UCB/project/gene expression data/GSE17755/healthy individual")
hi.file.list=list.files()

hi.cnames <- unlist(lapply(hi.file.list, function(x) strsplit(x, split="_")[[1]][[1]]))
hi.cnames

hi.gsm <- read.table(hi.file.list[1], sep="\t", fill=T, header=T, quote="")
hi.rnames <- hi.gsm[,2]

#read in the healthy control samples
for(i in 1:length(hi.file.list)){
  hi.gsm <- read.table(hi.file.list[i], sep="\t", fill=T, header=T, quote="")
  if(i==1){
    hi.red.df <- matrix(NA, nrow=nrow(hi.gsm), ncol=length(hi.file.list), dimnames=list(hi.rnames, hi.cnames))
    hi.green.df <- matrix(NA, nrow=nrow(hi.gsm), ncol=length(hi.file.list), dimnames=list(hi.rnames, hi.cnames))
    hi.red.bg.df <- matrix(NA, nrow=nrow(hi.gsm), ncol=length(hi.file.list), dimnames=list(hi.rnames, hi.cnames))
    hi.green.bg.df <- matrix(NA, nrow=nrow(hi.gsm), ncol=length(hi.file.list), dimnames=list(hi.rnames, hi.cnames))
  }
  
  if (length(which(hi.gsm[,2]==hi.rnames))!=length(hi.rnames)){
    print("probe names are differnet")
  } 
  
  hi.red.df[, i] <- hi.gsm[, 25]
  hi.green.df[, i] <- hi.gsm[, 26]
  hi.red.bg.df[, i] <- hi.gsm[, 27]
  hi.green.bg.df[, i] <- hi.gsm[, 28]
  print(i)
}

#order the probe names to ensure that thte probe names are identical
order.ind=unlist(lapply(hi.rnames, match, ds.rnames))

order.ds.red.df <- ds.red.df[order.ind, ]
order.ds.green.df <- ds.green.df[order.ind, ]
order.ds.red.bg.df <- ds.red.bg.df[order.ind, ]
order.ds.green.bg.df <- ds.green.bg.df[order.ind, ]

#red.log.df <- log(red.df, 2)
#green.log.df <- log(green.df, 2)
#
#red.bg.df <- log(red.bg.df, 2)
#green.bg.df <- log(green.bg.df, 2)
summary(hi.red.df)
summary(hi.green.df)

summary(ds.red.df)
summary(ds.green.df)

cnames=c(ds.cnames, hi.cnames)
rnames=rownames(hi.red.df)
  
red.df=cbind(order.ds.red.df, hi.red.df)
green.df=cbind(order.ds.green.df, hi.green.df)
red.bg.df=cbind(order.ds.red.bg.df, hi.red.bg.df)
green.bg.df=cbind(order.ds.green.bg.df, hi.green.bg.df)

#impute the missing values
# green.bg.impute <- round(impute.knn(green.bg.df)$data, 3)
# red.bg.impute <- round(impute.knn(red.bg.df)$data, 3)

green.bg.impute=na.omit(green.bg.df)
red.bg.impute=na.omit(red.bg.df)

green.df.impute=na.omit(green.df)
red.df.impute=na.omit(red.df)

# green.impute=na.omit(green.df)
# red.impute=na.omit(red.df)

#inspect the index of the missing value
# which(is.na(green.df), arr.ind=TRUE)
# which(is.na(red.df), arr.ind=TRUE)
# which(is.na(green.bg.df), arr.ind=TRUE)
# which(is.na(red.bg.df), arr.ind=TRUE)

# dim(red.bg.impute)
# dim(green.bg.impute)
# dim(green.impute)
# dim(red.impute)
# dim(green.df)

#green.df.impute <- round(impute.knn(green.df)$data, 3)
#red.df.impute <- round(impute.knn(red.df)$data, 3)

#filter the probes with low expression intensity
#calculate the gene level expression mean and sd
#green.bg.sd <- apply(green.bg.df, 1, function(x) sd(as.numeric(x)))
#red.bg.sd <- apply(red.bg.df, 1, function(x) sd(as.numeric(x)))
#green.bg.mean <- apply(green.bg.df, 1, function(x) mean(as.numeric(x)))
#red.bg.mean <- apply(red.bg.df, 1, function(x) mean(as.numeric(x)))
#
#green.thre <- green.bg.mean+2*green.bg.sd
#red.thre <- red.bg.mean+2*red.bg.sd
#green.filter <- green.df
#red.filter <- red.df

#up to here, not enough space.
# rm(list=c("red.filter", "green.filter", "red.bg.df", "green.bg.df")

#filter the probes with low expression intensity
green.bg.sd <- apply(green.bg.impute, 2, function(x) sd(as.numeric(x)))
red.bg.sd <- apply(red.bg.impute, 2, function(x) sd(as.numeric(x)))
green.bg.mean <- apply(green.bg.impute, 2, function(x) mean(as.numeric(x)))
red.bg.mean <- apply(red.bg.impute, 2, function(x) mean(as.numeric(x)))
green.thre <- green.bg.mean+2*green.bg.sd
red.thre <- red.bg.mean+2*red.bg.sd
green.filter <- green.df
red.filter <- red.df

for(i in 1:length(cnames)){
	green.filter[which(green.df[, i]<green.thre[i]), i] <- NA
	red.filter[which(red.df[, i]<red.thre[i]), i] <- NA
}

red.filter[1:10, 1:10]

gc()

rm(list=c("order.ds.green.bg.df", "order.ds.green.df", "order.ds.red.bg.df", "order.ds.red.df", "order.ind"))
rm(list=c("hi.cnames", "hi.file.list", "hi.green.bg.df", "hi.green.df", "hi.gsm", "hi.red.bg.df", "hi.red.df", "hi.rnames"))
rm(list=c("ds.cnames", "ds.file.list", "ds.green.bg.df", "ds.green.df", "ds.gsm", "ds.red.bg.df", "ds.red.df", "ds.rnames"))
   
#assemble the Red and Green matrix into RG_list
RG <- list(R=red.filter, G=green.filter, Rb=red.bg.df, Gb=green.bg.df)
RG_list <- as(RG, "RGList")

setwd("E:/UCB/project/gene expression data/GSE17755/")

RG_list$genes <- readGAL("RA.gal", fill=T)
RG_list$genes[1:30, ]
RG_list$printer <- getLayout(RG_list$genes)

jpeg(file="H:/Zhe Liu/project/gene expression data/GSE17755/entire/chip_intensity.jpg")
imageplot(log2(RG_list$R[,1]), RG_list$printer, low="white", high="red")
dev.off()

jpeg(file="H:/Zhe Liu/project/gene expression data/GSE17755/entire/chip_background.jpg")
imageplot(log2(RG_list$Rb[,1]), RG_list$printer, low="white", high="red")
dev.off()

boxplot(data.frame(cbind(red.filter[, 1:3], green.filter[, 1:3])),
		main="Background",col=rep(c("green","red"),each=3))

jpeg(file="H:/Zhe Liu/project/gene expression data/GSE17755/entire/RG_before_norm.jpg")
plotDensities(RG_list, log=T, arrays=NULL, singlechannels=NULL, groups=NULL, col=NULL)
dev.off()

rm(list=c("green.bg.df", "green.bg.impute", "green.bg.mean", "green.bg.sd", "green.df", "green.df.impute", "green.filter", "green.thre", "red.bg.df", "red.bg.impute", "red.bg.mean", "red.bg.sd", "red.df", "red.df.impute", "red.filter", "red.thre"))

#within array normalization
RG_norm_within=normalizeWithinArrays(RG_list, layout=RG_list$printer, method="printtiploess")#, layout=layout_param

jpeg(file="H:/Zhe Liu/project/gene expression data/GSE17755/entire/RG_norm_within.jpg")
plotDensities(RG_norm_within, log=T, arrays=NULL, singlechannels=NULL, groups=NULL, col=NULL)
dev.off()

#for(i in 1:ncol(MA_norm_within[[1]])){
#	plotMA(list(M=MA_norm_within[[1]][, i], A=MA_norm_within[[2]][, i]))
#}

#between array normalization
RG_norm_btw=normalizeBetweenArrays(RG_norm_within, method="quantile")

jpeg(file="H:/Zhe Liu/project/gene expression data/GSE17755/entire/RG_norm_btw.jpg")
plotDensities(RG_norm_btw, log=T, arrays=NULL, singlechannels=NULL, groups=NULL, col=NULL)
dev.off()

#dump the environment
setwd("H:/Zhe Liu/project/gene expression data/GSE17755/")
save.image(file = "entire current.RData")
# setwd("E:/UCB/project/gene expression data/GSE17755/report tip4/")

#impute the missing values
M.matrix <- na.omit(RG_norm_btw$M)
A.matrix <- na.omit(RG_norm_btw$A)

#MA plot
for (i in 1:length(cnames)){
	ma.dir <- paste("H:/Zhe Liu/project/gene expression data/GSE17755/entire/maplot/maplot", i, ".jpg", sep="")
	jpeg(file=ma.dir)
	ma.plot(A=A.matrix[, i], M=M.matrix[, i], cex=1)
	dev.off()
	print(i)
}

#box plot of the M matrix
# rownames(M.matrix) <- NULL
jpeg(file="H:/Zhe Liu/project/gene expression data/GSE17755/entire/ratio_boxplot.jpg")
boxplot(data.frame((M.matrix)), main="Intensity ratio boxplot",col=rep(c("green"),each=111))
dev.off()

# #filter the genes with more than 5% of the values are NA
# f1 <- function(x){ 
# 	if(length(which(is.na(x)))<(0.05*length(x))){
# 		return(x)
# 	}else{
#     return(NULL)
# 	}
# }
# 
# #inspect the number of missing values
# dim(na.omit(RG_norm_btw$M))
# 
# M.matrix.filter <- apply(RG_norm_btw$M, 1, f1)
# A.matrix.filter <- apply(RG_norm_btw$A, 1, f1)
# 
# probe.names=(names(M.matrix.filter)[which(unlist(lapply(M.matrix.filter,length))!=0)])
# 
# lapply(M.matrix.filter,length)
# 
# M.matrix.cbind <- do.call("rbind", M.matrix.filter)
# A.matrix.cbind <- do.call("rbind", A.matrix.filter)
# 
# rownames(M.matrix.cbind)=probe.names
# rownames(A.matrix.cbind)=probe.names
# 
# #differentially expressed gene analysis
# MA.list.filter=as(list(M=M.matrix.cbind, A=A.matrix.cbind), "MAList")
# 
# design=c(rep(1, ncol(M.matrix.cbind)))
# fit=lmFit(MA.list.filter, design)
# fit=eBayes(fit)
# topGene=topTable(fit, adjust="fdr",number=nrow(M.matrix.cbind))
# # length(which(rownames(d) %in% topGene[which(topGene[,6]<0.05), 1]))
# # topGene[5000, ]
# #use the top gene with FDR value less than 0.05
# topGene.fdr=topGene[which(topGene[,6]<0.05), ]
# M.clustering=M.matrix.cbind[which(rownames(M.matrix.cbind) %in% topGene.fdr[,1]), ]
# dim(M.clustering)
# 
# M.exprs <- round(impute.knn(M.matrix.cbind)$data, 3)
# A.exprs <- round(impute.knn(A.matrix.cbind)$data, 3)

#differentially expressed gene analysis
MA.list.filter=as(list(M=M.matrix, A=A.matrix), "MAList")

design=cbind(rep(1, ncol(M.matrix)), c(rep(1, 111), rep(0, 45)))

fit=lmFit(MA.list.filter, design)
fit=eBayes(fit)
fit$p.value[1:10, ]
M.matrix[1:10, 1:3]

topGene=topTable(fit, adjust="fdr",number=nrow(M.matrix))
# length(which(rownames(d) %in% topGene[which(topGene[,6]<0.05), 1]))
#use the top gene with FDR value less than 0.05
topGene.fdr=topGene[which(topGene[,6]<0.001), ]
dim(topGene.fdr)
M.clustering=M.matrix[which(rownames(M.matrix) %in% topGene.fdr[,1]), ]
dim(M.clustering)

rank.M.matrix=M.matrix[unlist(lapply(topGene[,1], match, rownames(M.matrix))), ]

#save the object for consensus clustering

setwd("H:/Zhe Liu/project/gene expression data/GSE17755/gene level clustering/clustering5")
save(M.clustering, topGene.fdr, file="M.exprs.entire.RData")

#M.matrix <- round(impute.knn(RG_norm_btw$M)$data, 3)
#A.matrix <- round(impute.knn(RG_norm_btw$A)$data, 3)

##quality assessment of the data
#eset <- ExpressionSet(assayData=as.matrix(M.matrix))
#eset
#
#setwd("E:/work/UCB/project/gene expression data/GSE17755/report tip4/")
#arrayQualityMetrics(eset, outdir = "within_btw_norm", force=T, do.logtransform=F, spatial=F)

#differentially expression gene analysis
#code to analyze the gene expression data
exp.data<- "E:/UCB/project/gene expression data/GSE17755/GSE17755_family.soft"
gds <- getGEO(filename=exp.data, GSEMatrix=T)
class(gds)  ## list -- ExpressionSet class (eSet)

gsmplatforms <- lapply(GSMList(gds),function(x) {Meta(x)$platform})
gsmplatforms

Table(GSMList(gds)[[1]])[1:5,]
# and get the column descriptions
Columns(GSMList(gds)[[1]])[1:5,]

fac.ds=do.call(rbind, lapply(GSMList(gds)[1:112], function(x) x@header))
fac.hi=do.call(rbind, lapply(GSMList(gds)[190:234], function(x) x@header))

fac <- rbind(fac.ds, fac.hi)
str(fac)
dim(fac)

ad.df <- as.data.frame(matrix(NA, nrow=nrow(fac), ncol=6))

for(i in 1:nrow(fac)){
	ad.ma1 <- do.call(rbind, lapply(fac[i, 2][[1]],  function(x) unlist(strsplit(x, split=": "))))
	ad.ma2 <- do.call(rbind, lapply(fac[i, 3][[1]],  function(x) unlist(strsplit(x, split=": "))))
	rownames(ad.ma1) <- paste("ch1", ad.ma1[,1], sep="@")
	ad.ma1 <- ad.ma1[,2]
	rownames(ad.ma2) <- paste("ch2", ad.ma2[,1], sep="@")
	ad.ma2 <- ad.ma2[,2]
	ad.ma <- cbind(t(ad.ma1), t(ad.ma2))
	ad.df[i, ] <- ad.ma
}

colnames(ad.df) <- colnames(ad.ma)
pdata <- cbind(fac, ad.df)
pheno <- as(pdata[-70, ],"AnnotatedDataFrame")

#load("E:/work/UCB/project/gene expression data/GSE17755/R data/current.RData")

#perform quality control
#arrayQualityMetrics(expressionset = gse17755, outdir = "QAnorm", force = TRUE, 
#    do.logtransform = F, spatial=F)

#contruct gene ID mapping matrix
feature <- ((GPLList(gds))[[1]]@dataTable@table)
feature <- data.frame(lapply(feature, as.character), stringsAsFactors=FALSE)

feature.id <- as.character(feature[, 1])
feature.matrix <- as.data.frame(matrix(NA, nrow=nrow(M.matrix), ncol=13))

rnames <- rownames(M.matrix)
rnames[1:10]
mat.len <- length(rnames)
feature.len <- length(feature.id)

name.list=as.list(rnames)
names(name.list)=seq(1:mat.len)

ptm <- proc.time()
for(i in 1:mat.len){
	if(length(which(feature.id %in% rnames[i]))!=0){
		feature.matrix[i, ]=feature[which(feature.id %in% rnames[i]), ]
	}
  print(i)
}

proc.time() - ptm

feature.matrix2 <- as.matrix(feature.matrix)
rownames(feature.matrix2) <- rownames(M.matrix)
colnames(feature.matrix2) <- colnames(feature)
feature.df <- as(as.data.frame(feature.matrix2),"AnnotatedDataFrame")

#construct expression set object
MA17755 <- as(list(M=M.matrix, A=A.matrix), "MAList")
gse17755 <- ExpressionSet(assayData=M.matrix, phenoData=pheno, featureData=feature.df)

#differentially expressed gene analysis 
pheno <- pData(gse17755)

# f <- (c(rep(0, nrow(pheno))))
# 
# fit = lmFit(MA17755)
# fit2 = eBayes(fit)
# 
# results = topTable(fit2, adjust = "BH", number=1000000)
# topgenes = results[results[, "P.Value"] < 0.001, ]
# topgenes[1:10, ncol(topgenes)]

featureNames(gse17755) <- feature.matrix2[, 5]
length(unique(feature.matrix2[, 5]))

#take the probe with largest Median Absolute Deviation of the log ratio matrix
entrez.factor <- as.factor(feature.matrix2[, 5])

M.mad <- apply(M.matrix, 1, mad)
M.mad[2]
M.index <- tapply(M.mad, entrez.factor, which.max)
M.index[2]

index <- tapply(1:nrow(M.matrix), entrez.factor, function(x) x)
index[2]

collapse.exprs <- t(mapply(function(x, y) M.matrix[x[y], ], index[2:length(index)], M.index[2:length(M.index)]))
collapse.exprs[1:3, 1:3]
dim(collapse.exprs)

write.table(collapse.exprs, file = "gse17755.entire.csv", quote = F, sep = "\t")
#dump the datasets and environment
setwd("H:/Zhe Liu/project/gene expression data/GSE17755/")
save(gse17755, file="gse17755.entire.RData")
save.image(file = "gse177755.entire.image.RData")
# load("gse177755.image.RData")

#up to here!!!!!
#the following code is to generate the files as the input of GSEA software.
rank.rnames=rownames(rank.M.matrix)

rank.rnames[1:10]
mat.len <- length(rank.rnames)
feature.len <- length(feature.id)

name.list=as.list(rank.rnames)
names(name.list)=seq(1:mat.len)

rank.feature.matrix <- data.frame(matrix("", nrow=nrow(M.matrix), ncol=13))

ptm <- proc.time()
for(i in 1:mat.len){
  if(length(which(feature.id %in% rank.rnames[i]))!=0){
      rank.feature.matrix[i, ]=feature[which(feature.id %in% rank.rnames[i]), ]
  }
}
proc.time() - ptm

rank.feature.matrix2 <- rank.feature.matrix
rownames(rank.feature.matrix2) <- rownames(rank.M.matrix)
colnames(rank.feature.matrix2) <- colnames(feature)
rank.gene=cbind(rank.feature.matrix2[, 5], topGene[, 4])

entrez.factor <- as.factor(feature.matrix2[, 5])
rank.gene.entrez <- tapply(rank.gene[,2], entrez.factor, which.max)
index <- tapply(1:nrow(rank.gene), entrez.factor, function(x) x)
index[2]

rank.index <- tapply(rank.gene[,2], entrez.factor, which.max)
rank.index[2]

collapse.rank.gene.list <- t(mapply(function(x, y) rank.gene[y[x], ], rank.index, index))
collapse.rank.gene.list[1:3, ]
collapse.rank.gene.list=cbind(collapse.rank.gene.list[,1], round(as.numeric(collapse.rank.gene.list[,2]),3))

write.table(na.omit(collapse.rank.gene.list), file="preranked gene list.entire.txt", row.names=F, col.names=F, quote=F, sep="\t")
chip.data=cbind(feature[, 1], feature[, 6], feature[, -c(1,6)])
chip.data[which(chip.data=="", arr.ind=T)]=NA
write.table(chip.data, file="gse17755.entire.chip", row.names=F, col.names=T, quote=F, sep="\t")



