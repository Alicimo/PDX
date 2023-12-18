x <- DGEList(data)
x <- calcNormFactors(x,method="TMM")
v <- voom(x,model.matrix(~0+rnaseq$meta$project))
v <- v$E

gl <- fread("~/Desktop/network_genes.txt",header = F)
gl.converted <- data.table(biomaRt::select(org.Hs.eg.db,keys=gl$V1,keytype = "SYMBOL",columns = c("ENSEMBL","ENTREZID","SYMBOL","GENENAME")))

v.gl <- v[match(gl.converted$ENSEMBL, rownames(v)),]

i <- !is.na(rownames(v.gl))
v.gl <- v.gl[i,]
gl.converted <- gl.converted[i,]

write.csv(v.gl,"~/Desktop/network_genes_expression.csv",row.names = T)