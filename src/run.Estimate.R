library(estimate)
library(data.table)
library(ggplot2)
library(ggthemes)

root.dir <- "~/OneDrive/projects/PDX/"

source("~/OneDrive/projects/PDX/src/load.PDX.R")
rnaseq <- get.PDX()
meta <- load.meta()

meta <- meta[project=="characterisation" & patient=="AB861"]
rnaseq <- subset.PDX(rnaseq, rnaseq$meta[,Sample.name] %in% meta[,Sample.name])

rnaseq <- merge.PDX(rnaseq)
rnaseq$meta <- merge(rnaseq$meta,meta,sort=F)
rnaseq$meta <- rnaseq$meta[,which(unlist(lapply(rnaseq$meta, function(x)!all(is.na(x))))),with=F]
rnaseq$meta[,merged.lanes:=abbreviate(merged.lanes)]

x <- rnaseq$human
rownames(x) <- mapIds(org.Hs.eg.db,keys=rownames(x),keytype = "ENSEMBL",column = "ENTREZID", multiVals = "first")
x <- x[!is.na(rownames(x)),]
x <- x[!duplicated(rownames(x)),]

write.table(x,file = "input.tsv",quote = F,sep="\t")
filterCommonGenes(input.f = "input.tsv", output.f = "output.gct", id = "EntrezID")
estimateScore("output.gct","score.gct",platform="illumina")
scores <- fread(paste(readLines("score.gct")[-(1:2)],collapse = "\n"))
file.remove(c("input.tsv","output.gct","score.gct"))

scores <- melt(scores[,-1],variable.name = "Sample.name")
scores$Sample.name <- as.character(scores$Sample.name)
scores <- merge(scores,rnaseq$meta[,.(Sample.name=gsub("-",".",Sample.name),tissue.grafted,tissue.sampled)],by="Sample.name")

ggplot(scores) + aes(x=Sample.name,y=value,fill=tissue.sampled) + geom_col(position="dodge") + facet_grid(Description~tissue.grafted,scales = "free_x",space = "free_x") + theme_tufte() + theme(axis.text.x=element_text(angle=90,hjust=1))
ggplot(scores[!is.na(tissue.grafted)]) + aes(x=tissue.grafted,y=value,fill=tissue.sampled) + geom_boxplot() + facet_grid(Description~.,scales = "free") + theme_tufte()
