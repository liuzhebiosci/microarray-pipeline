setwd("/Users/zheliu/Dropbox/work/Zhe Liu's meeting/thesis")
library("ArrayExpress")
AEset = ArrayExpress("E-MEXP-886")

save(AEset, file = "AEset.RData")
class(AEset)
AEset
colnames(pData(AEset))
fac = colnames(pData(AEset))[grep("Factor", colnames(pData(AEset)))]
fac
library("arrayQualityMetrics")
arrayQualityMetrics(expressionset = AEset, outdir = "QAraw", force = FALSE, do.logtransform = TRUE, intgroup = fac)
library("affy")
rAEset = rma(AEset)
arrayQualityMetrics(expressionset = rAEset, outdir = "QAnorm", force = FALSE, intgroup = fac)
cAEset = rma(AEset[, -1])
save.image(file = "introduction figure.RData")

#differentially expressed genes
#library("limma")
#groups = pData(cAEset)[, fac]
#groups[groups == "wildfitype"] = "WT" groups[groups == "ataxin 1 -/-"] = "KO" f = factor(groups)
#design = model.matrix(~f)
#colnames(design) = c("WTvsRef", "KOvsWT") design
#fit = lmFit(cAEset, design)
#fit2 = eBayes(fit)
#results = topTable(fit2, coef = "KOvsWT", adjust = "BH", number = nrow(cAEset))
#topgenes = results[results[, "P.Value"] < 0.001, ]
#colnames(results)
#m = exprs(cAEset[topgenes[, "ID"], ])

# functional analysis
#colnames(m) = 1:9
#colours = c("lightgreen", "mediumpurple")
#col = colours[f]
#heatmap(m, ColSideColors = col, margin = c(5, 8))
#library("GSEABase")
#annotation(cAEset)
#library("moe430a.db")
#gsc = GeneSetCollection(cAEset, setType = KEGGCollection()) Am = incidence(gsc)
#dim(Am)
#nsF = cAEset[colnames(Am), ] library("genefilter")
#rtt = rowttests(nsF, fac)
#rttStat = rtt$statistic
#Am2 = Am[selectedRows, ]
#tAadj = tA/sqrt(rowSums(Am2))
#names(tA) = names(tAadj) = rownames(Am2)
#library(KEGG.db)
#smPW = which(tAadj == min(tAadj))
#pwName = KEGGPATHID2NAME[[names(smPW)]] pwName
#KEGG2heatmap(names(smPW), nsF, "moe430a", col = colorRampPalette(c("white", + "darkblue"))(256), ColSideColors = col, margin = c(5, 8))
#sessionInfo()

