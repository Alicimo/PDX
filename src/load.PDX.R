pkgs <- c("data.table,limma,edgeR,Biobase,biomaRt")
invisible(lapply(strsplit(pkgs,",")[[1]],require,ch=T))


get.PDX <- function(purge=TRUE,normalise="voom",tissue.subset="human"){
  rnaseq <- load.data()
  #if(purge){ rnaseq <- purge.small(rnaseq) }
  #if(normalise=="voom"){ rnaseq <- normalise.voom(rnaseq) }
  #if(normalise=="RPKM"){ rnaseq <- normalise.RPKM(rnaseq) }
  rnaseq
}
  
load.data <- function(){
  rnaseq <- list()
  
  rnaseq$meta <- fread("~/OneDrive/projects/PDX/data/featureCounts.meta.tsv")
  
  counts <- fread("~/OneDrive/projects/PDX/data/featureCounts.counts.tsv")
  gene.list <- counts[,Geneid]
  gene.info <- counts[,1:6]
  counts <- as.matrix(counts[,-(1:6)])
  rownames(counts) <- gene.list
  
  if(any(rnaseq$meta[,fname] != colnames(counts))){stop()}
  
  rnaseq$mouse <- list(
    data=counts[substr(rownames(counts),1,4)=="ENSM",],
    gene.meta=gene.info[substr(rownames(counts),1,4)=="ENSM",]
  )
  rnaseq$human <- list(
    data=counts[substr(rownames(counts),1,4)=="ENSG",],
    gene.meta=gene.info[substr(rownames(counts),1,4)=="ENSG",]
  )

  rnaseq$meta[,mouse.library.size:=colSums(rnaseq$mouse$data)]
  rnaseq$meta[,human.library.size:=colSums(rnaseq$human$data)]
  
  #Return data
  rnaseq
}

load.meta <- function(){
  rbindlist(lapply(list.files("~/OneDrive/projects/PDX/data/meta/","*.csv",full.names = T),function(fname){
    x <- fread(fname,check.names = T,sep=",")
    x$project <- gsub(".csv","",basename(fname))
    x
  }),fill=T)
}

subset.meta <- function(meta,x){
  meta <- meta[x]
  meta <- meta[,which(unlist(lapply(meta, function(x)!all(is.na(x))))),with=F]
  meta <- meta[,which(unlist(lapply(meta, function(x)length(unique(x))!=1))),with=F]
  meta
}

unknown.meta <- function(){
  data <- get.PDX()
  meta <- load.meta()
  unknown <- data$meta[!(Sample.name %in% meta[project!="unknown",Sample.name]),.(SLX=paste(unique(Pool),collapse = ",")),Sample.name]
  fwrite(unknown,"~/OneDrive/projects/PDX/data/meta/unknown.csv")
}

normalise.voom <- function(rnaseq){
  #ßrnaseq$counts <- rnaseq$counts[rowMedians(rnaseq$counts) > 0,]
  rnaseq$counts <- rnaseq$counts[rowSums(cpm(rnaseq$counts) > 5) >= 10,]
  x <- DGEList(rnaseq$counts)
  plotMDS(x,labels=rnaseq$meta[,Sample.name])
  x <- calcNormFactors(x, method="TMM")
  d <- model.matrix(~0+interaction(rnaseq$meta$Flowcell,rnaseq$meta$Lane))
  rnaseq$counts <- as.matrix(voom(x,design=d,plot=TRUE)$E)
  rnaseq
}

normalise.RPKM <- function(rnaseq){
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  x <- getBM(attributes=c("ensembl_gene_id","start_position","end_position"),filter="ensembl_gene_id",values=rownames(rnaseq$counts),mart=ensembl)
  x$length <- x$end_position - x$start_position
  rnaseq$counts <- rnaseq$counts[rownames(rnaseq$counts) %in% x$ensembl_gene_id,]
  x <- x$length[match(rownames(rnaseq$counts),x$ensembl_gene_id)]
  rnaseq$counts <- rpkm(rnaseq$counts,x,log=T)
  rnaseq
}

subset.PDX <- function(rnaseq,x){
  PDX.subset <- list()
  PDX.subset$human <- list(
    gene.meta = rnaseq$human$gene.meta,
    data = rnaseq$human$data[,x]
  )
  PDX.subset$mouse <- list(
    gene.meta = rnaseq$mouse$gene.meta,
    data = rnaseq$mouse$data[,x]
  )
  PDX.subset$meta <- rnaseq$meta[x]
  if( dim(PDX.subset$meta)[1] != dim(PDX.subset$mouse$data)[2] ){stop()} 
  PDX.subset
}

merge.PDX <- function(rnaseq){
  merged <- list()
  merged$meta <- list()
  merged$human <- list()
  merged$mouse <- list()
  
  for (sname in rnaseq$meta[,unique(Sample.name)]){
    x <- rnaseq$meta[,Sample.name==sname]
    merged$meta <-c(merged$meta,
                    list(data.table(
                      Sample.name=sname,
                      merged.lanes=paste(sort(rnaseq$meta[x,paste(Pool,Flowcell,Lane,sep=":")]),collapse=";")
                    )))
    merged$human <- c(merged$human,list(rowSums(rnaseq$human$data[,x,drop=F])))
    merged$mouse <- c(merged$mouse,list(rowSums(rnaseq$mouse$data[,x,drop=F])))
  }
  
  merged$meta <- rbindlist(merged$meta)
  merged$human <- do.call(cbind,merged$human)
  merged$mouse <- do.call(cbind,merged$mouse)
  
  rownames(merged$human) <- rownames(rnaseq$human$data)
  colnames(merged$human) <- merged$meta[,Sample.name]
  rownames(merged$mouse) <- rownames(rnaseq$mouse$data)
  colnames(merged$mouse) <- merged$meta[,Sample.name]
  
  merged$meta[,mouse.library.size:=colSums(merged$mouse)]
  merged$meta[,human.library.size:=colSums(merged$human)]
  
  merged$mouse <- list(
    gene.meta=rnaseq$mouse$gene.meta,
    data=merged$mouse
      )
  merged$human <- list(
    gene.meta=rnaseq$human$gene.meta,
    data=merged$human
  )
  
  
  merged
}