library(data.table)

micro <- list()
micro$counts <- fread("data/ExpressionSamples.txt")

#Process microarray data
gene.list <- micro$counts[,Gene]
micro$counts[,Gene:=NULL]
micro$counts <- as.matrix(micro$counts)
rownames(micro$counts) <- gene.list

#Parse microarray IDs and generate meta
micro$meta <- data.table(Sample.name=colnames(micro$counts))
micro$meta[,patient:=tstrsplit(Sample.name,"-")[1]]
micro$meta[,primary:=substr(patient,nchar(patient),nchar(patient))!='M']
micro$meta[primary==FALSE,patient:=substr(patient,1,nchar(patient)-1)]
micro$meta[,model:=tstrsplit(Sample.name,"-")[2]]
micro$meta[,model:=tstrsplit(model,"R")[1]]
micro$meta[,xenograft:=substr(model,1,1)=='X']
micro$meta[,cells:=lapply(tstrsplit(model,"C")[2],function(x)!is.na(x))]