options(stringsAsFactors = FALSE) 

library("impute")

setwd("C:/work/microarray/DevGSEA/dataset/expr data/raw data")

f1<-read.table("Table_NKI_295_1.txt",sep="\t",header=TRUE)
f2<-read.table("Table_NKI_295_2.txt",sep="\t",header=TRUE)
f3<-read.table("Table_NKI_295_3.txt",sep="\t",header=TRUE)
f4<-read.table("Table_NKI_295_4.txt",sep="\t",header=TRUE)
f5<-read.table("Table_NKI_295_5.txt",sep="\t",header=TRUE)
f6<-read.table("Table_NKI_295_6.txt",sep="\t",header=TRUE)

combined_matrix<-cbind(f1,f2[3:ncol(f2)],f3[3:ncol(f3)],f4[3:ncol(f4)],f5[3:ncol(f5)],f6[3:ncol(f6)])
combined_matrix=combined_matrix[-1,]
dim(combined_matrix)
col_names=gsub("Sample.", "", colnames(combined_matrix)[grep("Sample", colnames(combined_matrix), )])

#substract those flag is 0 which means that they are control probes or bad splot
flag<-combined_matrix[, seq(1:295)*5+2]
flag=apply(flag, c(1,2), as.numeric)
which(rowsum(flag)!=295)

#there is errors here but they will not affect the results
sel_flag<- apply(flag,1,function(x) any(x)!=0)
array_sel<-combined_matrix[sel_flag,]
ratio_array=array_sel[,c(1,2,seq(1:295)*5-2)]

rownames(ratio_array)=ratio_array[, 1]
ratio_array=ratio_array[,-c(1,2)]
colnames(ratio_array)=col_names
dim(na.fail(ratio_array))

#log 2 transformation of the ratio data
log_ratio_array=apply(ratio_array, c(1,2), function(x) log(10^as.numeric(x), 2))

#impute the missing value
log_ratio_array<- round(impute.knn(log_ratio_array)[[1]], 3)
write.table(log_ratio_array, "log_ratio_imputed.txt", quote=FALSE, row.names=T, col.names=T, sep="\t")

#get the intensity data
intensity_array=array_sel[,c(1,2,seq(1:295)*5+1)]
rownames(intensity_array)=intensity_array[, 1]
intensity_array=intensity_array[,-c(1,2)]
colnames(intensity_array)=col_names
dim(na.fail(intensity_array))

#log 2 transformation of the intensity data
log_intensity_array=apply(intensity_array, c(1,2), function(x) log(10^as.numeric(x), 2))

#impute the missing value
log_intensity_array<- round(impute.knn(log_intensity_array)[[1]], 3)
write.table(log_intensity_array, "log_intensity_imputed.txt", quote=FALSE, row.names=T, col.names=T, sep="\t")

################################################################################
library("limma")
library("arrayQualityMetrics")

#load in the data for pre-processing
setwd("C:/work/microarray/DevGSEA/dataset/expr data/raw data")
log_ratio_imputed=read.table("log_ratio_imputed.txt", sep="\t", header=T, check.names=FALSE)
log_intensity_imputed=read.table("log_intensity_imputed.txt", sep="\t", header=T, check.names=FALSE)

#rownames(log_ratio_imputed)=log_ratio_imputed[,1]
#rownames(log_intensity_imputed)=log_intensity_imputed[,1]
#
#log_ratio_imputed[1:10,1:10]
#log_intensity_imputed[1:10,1:10]
#
#log_ratio_imputed=log_ratio_imputed[,-1]
#log_intensity_imputed=log_intensity_imputed[,-1]

#####map entrez gene ID and gene name
mapping<-read.table("ID map.txt", sep="\t", header=TRUE)
index_list<-list()
cache <- new.env(hash=TRUE, parent=emptyenv())

for (i in 1:nrow(mapping)){
	assign(as.character(mapping[i,1]), i, envir=cache)
}

index_list<- mget(as.character(rownames(log_ratio_imputed)), cache, ifnotfound=NA)
index<- do.call("rbind", index_list)

#select the non-NA value of the log ratio matrix
log_ratio_imputed_map<- cbind(mapping[index[,1],2:3], log_ratio_imputed)
log_ratio_imputed_sel<-log_ratio_imputed_map[!is.na(log_ratio_imputed_map[,1]),]

entrez_factor=as.factor(log_ratio_imputed_sel[, 2])

#take the median value of the log ratio matrix
output_log_ratio_list=list()

for(i in 3:ncol(log_ratio_imputed_sel)){
	output_log_ratio_list[[i-2]]=tapply(log_ratio_imputed_sel[, i], entrez_factor, median)
}

output_log_ratio_df=do.call("cbind", output_log_ratio_list)
colnames(output_log_ratio_df)=colnames(log_ratio_imputed)
dim(output_log_ratio_df)
output_log_ratio_df[1:5,1:5]

#select the non-NA value of the log intensity matrix
log_intensity_imputed_map<- cbind(mapping[index[,1],2:3], log_intensity_imputed)
log_intensity_imputed_sel<-log_intensity_imputed_map[!is.na(log_intensity_imputed_map[,1]),]

#take the median value of the log intensity matrix
output_log_intensity_list=list()

for(i in 3:ncol(log_intensity_imputed_sel)){
	output_log_intensity_list[[i-2]]=tapply(log_intensity_imputed_sel[, i], entrez_factor, median)
}

output_log_intensity_df=do.call("cbind", output_log_intensity_list)

colnames(output_log_intensity_df)=colnames(log_intensity_imputed)
dim(output_log_intensity_df)
output_log_intensity_df[1:5,1:5]

#quality assessment of the data
pheno=read.AnnotatedDataFrame("descrp.txt", head=T, sep="\t", fill=T)
eset1<- new("ExpressionSet", exprs=as.matrix(output_log_intensity_df), phenoData=pheno)
eset1

fac=colnames(pData(eset1))[grep("EVENTmeta", colnames(pData(eset1)))][[1]]

arrayQualityMetrics(expressionset = eset1, outdir = "QAraw_new3", force=T, do.logtransform=F, intgroup=fac, spatial=F) 

#normalization
MA=list(M=output_log_ratio_df, A=output_log_intensity_df)
MA_list=as(MA, "MAList")

#quality assessment of the data by limma QA function
for(i in 1:ncol(MA_list[[1]])){
	plotMA(list(M=MA_list[[1]][, i], A=MA_list[[2]][, i]))
}

pdf(file="C:/work/microarray/DevGSEA/dataset/expr data/raw data/MA_before_norm_new.pdf")
plotDensities(MA_list, log=T, arrays=NULL, singlechannels=NULL, groups=NULL, col=NULL)
dev.off()

layout_param<- list(ngrid.r=1, ngrid.c=1,
		nspot.r=nrow(MA[[1]]), nspot.c=1)

#within array normalization
MA_norm_within=normalizeWithinArrays(MA_list, layout=layout_param, method="loess")

for(i in 1:ncol(MA_norm_within[[1]])){
	plotMA(list(M=MA_norm_within[[1]][, i], A=MA_norm_within[[2]][, i]))
}

pdf(file="C:/work/microarray/DevGSEA/dataset/expr data/raw data/MA_norm_within_new.pdf")
plotDensities(MA_norm_within, log=T, arrays=NULL, singlechannels=NULL, groups=NULL, col=NULL)
dev.off()

#between array normalization
MA_norm_btw=normalizeBetweenArrays(MA_norm_within, method="quantile")

for(i in 1:ncol(MA_norm_btw[[1]])){
	plotMA(list(M=MA_norm_btw[[1]][, i], A=MA_norm_btw[[2]][, i]))
}

pdf(file="C:/work/microarray/DevGSEA/dataset/expr data/raw data/MA_norm_btw_new.pdf")
plotDensities(MA_norm_btw, log=T, arrays=NULL, singlechannels=NULL, groups=NULL, col=NULL)
dev.off()

#quality assessment of the data
pheno=read.AnnotatedDataFrame("descrp.txt", head=T, sep="\t", fill=T)
eset2 <- new("ExpressionSet", exprs=as.matrix(MA_norm_btw[[2]]), phenoData=pheno)
eset2

fac=colnames(pData(eset2))[grep("EVENTmeta", colnames(pData(eset2)))][[1]]

arrayQualityMetrics(expressionset = eset2, outdir = "QAraw_new4", force=FALSE, do.logtransform=F, intgroup=fac, spatial=F) 
output_file=round(exprs(eset2), 3)

write.table(output_file, "new_exprs_entrez_id.txt", quote=FALSE, row.name=T, col.name=T, sep="\t")

#substract ER positive/negative data

descrp=read.table("descrp.txt", sep="\t", head=T)

er_pos_exprs=output_file[, which(descrp[,10]==1)]
er_neg_exprs=output_file[, which(descrp[,10]==0)]

#output_pos_data<- round(t(apply(er_pos_exprs, 1, z_score)), 3)
#output_neg_data<- round(t(apply(er_neg_exprs, 1, z_score)), 3)

#apply(output_pos_data, c(1), mean)
#apply(output_pos_data, c(1), sd)
#
#apply(output_neg_data, c(1), mean)
#apply(output_neg_data, c(1), sd)

#output_pos_file=cbind(array_id_sel[sub_ind, 1], output_pos_data)
#output_neg_file=cbind(array_id_sel[sub_ind, 1], output_neg_data)

#dim(output_pos_file)
#dim(output_neg_file)

er_pos_df=descrp[which(descrp[,10]==1), ]
er_neg_df=descrp[which(descrp[,10]==0), ]

write.table(er_pos_exprs, "new_exprs_entrez_er_pos.txt", quote=FALSE, row.name=T, col.name=T, sep="\t")
write.table(rbind(er_pos_df[, 4]), "new_exprs_entrez_er_pos_cls.txt", quote=FALSE, row.name=F, col.name=F, sep="\t")

write.table(er_neg_exprs, "new_exprs_entrez_er_neg.txt", quote=FALSE, row.name=T, col.name=T, sep="\t")
write.table(rbind(er_neg_df[, 4]), "new_exprs_entrez_er_neg_cls.txt", quote=FALSE, row.name=F, col.name=F, sep="\t")

#dump the data for further usage
dump_list=c("er_pos_exprs", "er_neg_exprs", "er_pos_df", "er_neg_df")
dump(dump_list, file = "dump_exprs.R", append = FALSE, 
		control = "all", envir = parent.frame(), evaluate = TRUE)

###########################
#output as gsea input file
###########################

#entire dataset
setwd("C:/work/microarray/DevGSEA/input data/New data")

write("#1.2", "new_exprs_entrez_id.gct", append=F, sep="\t")
write(paste(nrow(output_file), ncol(output_file), sep="\t"), "new_exprs_entrez_id.gct", append=T, sep="\t")
write.table(cbind(rep(NA, nrow(output_file)), output_file), "new_exprs_entrez_id.gct", quote=F, append=T, sep="\t")

#####entire dataset class file
write(paste(ncol(output_file), 2, 1, sep=" "), "new_exprs_entrez_id.cls", append=F, sep="\t")
write("#free metasis", "new_exprs_entrez_id.cls", append=T, sep="\t")

cls=gsub("0", "free", descrp[, 4])
cls=gsub("1", "metasis", cls)
cat(cls, file="new_exprs_entrez_id.cls", sep=" ", fill=FALSE, labels=NULL, append=T)






#compute the differentailly expressioin, choose the most 
#significantly expressed genes
#header<- read.table("colname.txt",sep="\t")
#header[1:10]
#
#t_test_fun=function(x){
#	free_ind=which(header=="free")
#	meta_ind=which(header=="metasis")
#	
#	return(t.test(x[which(header=="free")], x[which(header=="metasis")])$p.value)
#}
#
#t_stat=apply(norm_array, 1, t_test_fun)
#
######obtain the minimium t-test p-value transcript to represent gene
#entrez.factor=as.factor(array_id_sel[, 2])
#
#sub_t_stat=t_stat[!is.na(array_id[,1])]
#names(sub_t_stat)=array_id_sel[, 2]
#
#output.list=list()
#output.list=tapply(sub_t_stat, entrez.factor, which.min)
#
#sub_ind=c()
#
#for(i in 1:length(output.list)){
#	sub_ind=c(sub_ind, which(names(output.list[i])==names(sub_t_stat))[output.list[i]])
#}
#
#output_file=array_id_sel[sub_ind,-c(1,2)]
#rownames(output_file)=array_id_sel[sub_ind,2]

#output_norm_file<- round(t(apply(output_file, 1, z_score)),3)
#output_norm_file<- cbind(array_id_sel[sub_ind, 1], output_norm_file)
#
#apply((output_norm_file[, -c(1)]), c(1), function(x) mean(as.numeric(x)))
#apply((output_norm_file[, -c(1)]), c(1), function(x) sd(as.numeric(x)))
#
#write.table(output_norm_file, "new_exprs_entrez_id.txt", quote=FALSE, row.name=T, col.name=T, sep="\t")
#
