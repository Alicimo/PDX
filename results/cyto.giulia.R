pkgs <- c("data.table,limma,edgeR,GSVA,org.Hs.eg.db,org.Mm.eg.db,ggplot2,ggthemes,pheatmap,biomaRt")
invisible(lapply(strsplit(pkgs,",")[[1]],require,ch=T))

rnaseq <- readRDS("rnaseq.RDS")
cytos <- fread("uniquecytokinelist.csv")

ensembl.human = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
ensembl.mouse = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

norm.count.matrix <- function(y,zero.count.shift=0.5,lib.sizes){
  y <- t(1E6 * t(y) / (lib.sizes * calcNormFactors(y,lib.sizes)))
  y <- log(y+zero.count.shift)
  return(y)
}

keyMap <- getBM(
  attributes=c('hgnc_symbol','ensembl_gene_id'),
  filters = 'hgnc_symbol', 
  values = cytos$Human_gene, 
  mart = ensembl.human
)

human.cytos <- unique(merge(cytos[,.(Cytokine,Human_gene)],keyMap,by.x="Human_gene",by.y="hgnc_symbol",all.x=T,sort=F))
human.cytos <- human.cytos[-which(ensembl_gene_id %in% c('ENSG00000228978','ENSG00000204490','ENSG00000228849','ENSG00000206439','ENSG00000228321','ENSG00000230108','ENSG00000223952','ENSG00000271503'))]

keyMap <- getBM(
  attributes=c('mgi_symbol','ensembl_gene_id'),
  filters = 'mgi_symbol', 
  values = cytos$Mouse_gene, 
  mart = ensembl.mouse
)

mouse.cytos <- unique(merge(cytos[,.(Cytokine,Mouse_gene)],keyMap,by.x="Mouse_gene",by.y="mgi_symbol",all.x=T,sort=F))
mouse.cytos <- mouse.cytos[-which(ensembl_gene_id=="ENSMUSG00000115982")]

## Boxplot

x <- norm.count.matrix(rnaseq$human$data,0.5,colSums(rnaseq$human$data))
x <- x[match(human.cytos$ensembl_gene_id,rownames(x)),]
rownames(x) <-human.cytos$Cytokine
x <- cbind(rnaseq$meta,t(x))
x <- melt(x,id.vars = 1:ncol(rnaseq$meta))
x[,compartment:="human"]

y <- norm.count.matrix(rnaseq$mouse$data,0.5,colSums(rnaseq$mouse$data))
y <- y[match(mouse.cytos$ensembl_gene_id,rownames(y)),]
rownames(y) <- mouse.cytos$Cytokine
y <- cbind(rnaseq$meta,t(y))
y <- melt(y,id.vars = 1:ncol(rnaseq$meta))
y[,compartment:="mouse"]

z <- rbind(x,y)

(p <- ggplot(z) + aes(x=variable,y = value,fill=sample.type) + geom_boxplot() + facet_grid(compartment~.) + theme_bw() + theme(axis.text.x = element_text(angle=90,hjust=1)) + labs(x="",y="Expression",fill="") + ggtitle("Cytokines"))

## Heatmap

x <- norm.count.matrix(rnaseq$human$data,0.5,colSums(rnaseq$human$data))
x <- x[match(human.cytos$ensembl_gene_id,rownames(x)),]
rownames(x) <- human.cytos$Cytokine
x <- x[!is.na(rowSums(x)),]
x <- x[hclust(dist(x))$order,]

y <- norm.count.matrix(rnaseq$mouse$data,0.5,colSums(rnaseq$mouse$data))
y <- y[match(mouse.cytos$ensembl_gene_id,rownames(y)),]
rownames(y) <- paste0("m.",mouse.cytos$Cytokine)
y <- y[!is.na(rowSums(y)),]
y <- y[hclust(dist(y))$order,]

ann <- data.frame(rnaseq$meta[,.(sample.type,merged.lanes)])
rownames(ann) <- rnaseq$meta$Sample.name

pheatmap(rbind(x,y),annotation_col = ann,cluster_rows = F,gaps_row = rep(nrow(x),3))