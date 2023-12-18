source("~/OneDrive/projects/PDX/src/load.PDX.R")
rnaseq <- get.PDX()
meta <- load.meta()
meta <- meta[project=="collaborators"]
rnaseq <- subset.PDX(rnaseq, rnaseq$meta[,Sample.name] %in% meta[,Sample.name])
meta <- meta[,which(unlist(lapply(meta, function(x)!all(is.na(x))))),with=F]
rnaseq <- subset.PDX(rnaseq,rnaseq$meta[,Sample.name!="17_R1  251  CNTL 1R (7752)"])
rnaseq <- merge.PDX(rnaseq)
rnaseq$meta <- merge(rnaseq$meta,meta,sort=F)
rnaseq$meta[,merged.lanes:=NULL]
rnaseq$meta[,project:=NULL]
x <- rnaseq$human
write.csv(x,file="~/Desktop/collabs.counts.csv")

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

db <- org.Hs.eg.db
tx <- TxDb.Hsapiens.UCSC.hg19.knownGene

map <- select(db, keys = rownames(x), columns =c("ENTREZID"), keytype="ENSEMBL")
map <- map[!is.na(map[,2]),]
map <- map[!(duplicated.default(map$ENSEMBL)),]
map$start <- mapIds(tx,column = "TXSTART",keys = map[,2],keytype = "GENEID",multiVals = function(x){min(x)})
map$end <- mapIds(tx,column = "TXEND",keys = map[,2],keytype = "GENEID",multiVals = function(x){max(x)})
map$len <- map$end - map$start
map <- map[!is.na(map[,5]),]

x <- x[rownames(x) %in% map$ENSEMBL,]
map <- map[match(map$ENSEMBL,rownames(x)),]

FPKM <- t(t(x / colSums(x)) / map$len) * 1E9
write.csv(FPKM,file="~/Desktop/collabs.FPKM.csv")