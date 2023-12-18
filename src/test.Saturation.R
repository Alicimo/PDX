library(estimate)
library(data.table)
library(ggplot2)
library(ggthemes)

root.dir <- "~/OneDrive/projects/PDX/"

source(paste0(root.dir,"src/load.PDX.R"))
rnaseq <- get.PDX()
meta <- load.meta()

meta <- meta[project=="characterisation" & patient=="AB861"]
rnaseq <- subset.PDX(rnaseq, rnaseq$meta[,Sample.name] %in% meta[,Sample.name])

rnaseq <- merge.PDX(rnaseq)
rnaseq$meta <- merge(rnaseq$meta,meta,sort=F)
rnaseq$meta <- rnaseq$meta[,which(unlist(lapply(rnaseq$meta, function(x)!all(is.na(x))))),with=F]
rnaseq$meta[,merged.lanes:=abbreviate(merged.lanes)]

human.saturation <- apply(rnaseq$human,2,function(x){
  x <- factor(sample(rep(1:length(x),times=x)))
  N <- length(x)
  y <- lapply(0:floor(N/1E6), function(i){
    table(x[(i*1E6+1):min((i+1)*1E6,N)])
  })
  y <- colSums(t(apply(Reduce(cbind,y),1,cumsum)) > 5)
  data.table(reads=c(seq(1E6,N,1E6),N),genes.detected=y)
})
human.saturation <- rbindlist(lapply(1:length(human.saturation),function(i){melt(as.data.table(human.saturation[i]),id.vars = 1)}))
names(human.saturation) <- c("N.reads","Sample.name","Genes.detected")
human.saturation[,Sample.name:=gsub(".genes.detected","",Sample.name)]

ggplot(human.saturation) + aes(x=N.reads,y=Genes.detected,group=Sample.name) + geom_line() + geom_point() + theme_tufte(20)
ggplot(human.saturation[,head(.SD,-1),Sample.name][,.(tail(N.reads,-1),tail(Genes.detected,-1)/head(Genes.detected,-1)-1),Sample.name]) + aes(x=V1,y=V2,group=Sample.name) + geom_line() + geom_point() + theme_tufte(20)

##### mouse

mouse.saturation <- apply(rnaseq$mouse,2,function(x){
  x <- factor(sample(rep(1:length(x),times=x)))
  N <- length(x)
  y <- lapply(0:floor(N/1E5), function(i){
    table(x[(i*1E5+1):min((i+1)*1E5,N)])
  })
  if(length(y)>1){
    y <- colSums(t(apply(Reduce(cbind,y),1,cumsum)) > 5)
  } else {
    y <- sum(y[[1]] > 5)
  }
  data.table(reads=c(seq_len(length(y)-1)*1E5,N),genes.detected=y)
})
mouse.saturation <- rbindlist(lapply(1:length(mouse.saturation),function(i){melt(as.data.table(mouse.saturation[i]),id.vars = 1)}))
names(mouse.saturation) <- c("N.reads","Sample.name","Genes.detected")
mouse.saturation[,Sample.name:=gsub(".genes.detected","",Sample.name)]

ggplot(mouse.saturation) + aes(x=N.reads,y=Genes.detected,group=Sample.name) + geom_line() + geom_point() + theme_tufte(20)
ggplot(mouse.saturation[,head(.SD,-1),Sample.name][,.(tail(N.reads,-1),tail(Genes.detected,-1)/head(Genes.detected,-1)-1),Sample.name]) + aes(x=V1,y=V2,group=Sample.name) + geom_line() + geom_point() + theme_tufte(20)


mouse.model <- mouse.saturation[,lm(Genes.detected~poly(log(N.reads),2))]
mouse.extra.lanes.frac <- rnaseq$meta[merged.lanes=="SLX-14085" & xenograft==TRUE,as.list(predict(mouse.model,data.frame(N.reads=c(mouse.library.size,mouse.library.size*(4/3),mouse.library.size*(5/3),mouse.library.size*2)))),Sample.name][,100*.SD[,2:4]/min(.SD)-100,Sample.name]
names(mouse.extra.lanes.frac) <- c("Sample.name","1.extra","2.extra","3.extra")
mouse.extra.lanes.genes <- rnaseq$meta[merged.lanes=="SLX-14085" & xenograft==TRUE,as.list(predict(mouse.model,data.frame(N.reads=c(mouse.library.size,mouse.library.size*(4/3),mouse.library.size*(5/3),mouse.library.size*2)))),Sample.name][,floor(.SD[,2:4] - min(.SD)),Sample.name]
names(mouse.extra.lanes.genes) <- c("Sample.name","1.extra","2.extra","3.extra")
