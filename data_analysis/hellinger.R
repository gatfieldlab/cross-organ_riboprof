####      HELLINGER DISTANCES OF RELATIVE TRANSCRIPT EXPRESSION     ####
##### analysis,using ONLY TRANSLATED/PROTEIN CODING TRANSCRIPTS (so transcripts present in k.annot(prepared.cds) and not all transcripts detected in cufflinks(that also contain processed transcripts and NMD))
# input is output of get_fmi.py
fmi.k.input <- read.table('fmi_kidney.txt', header=FALSE)
fmi.k <- fmi.k.input[which(fmi.k.input$V2 %in% k.annot[which(k.annot$tr.type == 'KNOWN'),'trID']),]
fmi.tot <- aggregate(V3~V1, data=fmi.k, FUN='sum')
fmi.k <- merge(fmi.k, fmi.tot, by='V1')
fmi.k$V4 <- fmi.k$V3.x/fmi.k$V3.y
fmi.k <- fmi.k[-4]
colnames(fmi.k) <- c('gene','transcript','fmi','rel_exp')

fmi.l.input <- read.table('fmi_liver.txt', header=FALSE)
fmi.l <- fmi.l.input[which(fmi.l.input$V2 %in% l.annot[which(l.annot$tr.type == 'KNOWN'),'trID']),]
fmi.tot <- aggregate(V3~V1, data=fmi.l, FUN='sum')
fmi.l <- merge(fmi.l, fmi.tot, by='V1')
fmi.l$V4 <- fmi.l$V3.x/fmi.l$V3.y
fmi.l <- fmi.l[-4]
colnames(fmi.l) <- c('gene','transcript','fmi','rel_exp')


fmi.kl <- merge(fmi.k, fmi.l, by=c('gene','transcript'), all=TRUE)
fmi.kl <- fmi.kl[-c(3,5)]
fmi.kl <- cbind(fmi.kl[,1:2], apply(fmi.kl[,3:4], c(1,2), function(x) if(is.na(x)) 0 else x))
colnames(fmi.kl) <- c('gid','trid','rel_exp_kidney','rel_exp_liver')
rel.exp.k <- aggregate( rel_exp_kidney ~ gid, fmi.kl, as.vector)
rel.exp.l <- aggregate( rel_exp_liver ~ gid, fmi.kl, as.vector)
#hell.gene.list <- intersect(rel.exp.k$gid, rel.exp.l$gid) #instead I will do the hell.dist matrix on kl.expressed.ids.
rownames(rel.exp.k) <- rel.exp.k$gid
rel.exp.k <- rel.exp.k[-1]
rownames(rel.exp.l) <- rel.exp.l$gid
rel.exp.l <- rel.exp.l[-1]

hell.dist <- data.frame(hell_dist = rep(NA, length(kl.expressed.ids)), row.names = kl.expressed.ids)
for(gid in kl.expressed.ids){
  hell.dist[gid,] <- sqrt(1- sum(sqrt(rel.exp.k[gid,][[1]]*rel.exp.l[gid,][[1]])))
}
