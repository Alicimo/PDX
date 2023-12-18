library(data.table)
meta <- fread("~/Desktop/PDX.rnaseq.char.meta.csv")
rnaseq <- fread("~/Desktop/PDX.rnaseq.char.human.TPMs_corrected.csv")

x <- list()
for(p.ID in sort(meta[,length(unique(sample.type)),patient][V1>1,patient])){
  for(priID in meta[patient==p.ID & sample.type=="primary",unique(tissue.sampled)]){
    for(xenID in meta[patient==p.ID & sample.type=="xenograft",unique(tissue.sampled)]){
      for(pasID in meta[patient==p.ID & sample.type=="xenograft" & tissue.sampled==xenID,unique(passage)]){
      
        i <- meta[,.I[patient==p.ID & sample.type=="primary" & tissue.sampled==priID]]
        j <- meta[,.I[patient==p.ID & sample.type=="xenograft" & tissue.sampled==xenID & passage==pasID]]
        
        primary <- rowMeans(rnaseq[,i+1,with=F])
        xeno <- rowMeans(rnaseq[,j+1,with=F])
  
        x[[length(x)+1]] <- data.table(patient=p.ID,primary=priID,xenograft=xenID,passage=paste0("x",pasID),cor=cor(primary,xeno))
      }
    }
  }
}
x <- rbindlist(x)

