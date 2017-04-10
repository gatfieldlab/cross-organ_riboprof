####     CIRCADIAN ANALYSIS KIDNEY     ####
writeLines('   ... creating DESeq count data-sets and estimating dispersions')
cd.k.rp <- get_cds(round(k.rp[k.rp_filter, k.rp.cds]), design, norm.factors.k.rp)
cd.k.tr <- get_cds(round(k.tr[k.tr_filter, k.tr.cds]), design, norm.factors.k.tr)
theta.k.rp <- get_dispersions(cd.k.rp)
theta.k.tr <- get_dispersions(cd.k.tr)


writeLines('Performing circadian analysis, this will take time ...')

tpoints <- rep(seq(0,22,2),each=2)

writeLines('Cyclic model selection based on AIC')
for (type in c('k.rp', 'k.tr')) {
  norm.data <- get(paste('norm', 'cds', type, sep='.'))
  theta <- get(paste('theta', type, sep='.'))
  writeLines(paste('   ... calculating', type))
  d <- apply(cbind(norm.data, theta),1, function(x) {
    fits <- fitModel(x[1:length(tpoints)], tpoints, weights=rep(1,24)) #theta=x[length(x)])
    weights.aic <- waic(delta.aic.all(fits))  #originally: waic(delta.aic.all2(fits)) Mod Aug2015
    #weights.aic <- c(weights.aic[1:2], sum(weights.aic[3:4]))
    pars.sigmoid <- get_sigmoid_pars(fits$sigmoid)
    pars.sinusodial <- get_params(fits$sinusodial, 24)
    pars.null <- coef(fits$null)[[1]]
    c(pars.sigmoid, pars.sinusodial, pars.null, weights.aic, fits$weight_mode)
  })
  d <- data.frame(t(d))
  for (i in c(1:19)) { d[,i] <- as.numeric(as.vector(d[,i]))}
  colnames(d) <- c(paste('sigmoid', c('fc', 'p', 'l', 'm', 'B', 'min', 'max', 'ratio', 'phase'), sep='.'), 
                   paste('sinusoidal', c('mean', 'amp', 'min', 'max', 'ratio', 'phase'), sep='.'),
                   'null.intercept', paste('waic', c('sigmoid', 'sinusoidal', 'null'), sep='.'), 'weight_mode') 
  assign(paste('model', type, sep='.'), d)
}

writeLines('Detecting circadian transcripts with given cutoffs from models')
sig.cutoff <- 3.5
ratio.cutoff <- 1.5
for (type in c('k.rp', 'k.tr')) {
  cur.model <- get(paste('model', type, sep='.'))
  d <- apply(cur.model, 1, function(x) {
    weights.aic <- as.numeric(x[17:19])
    if (weights.aic[3] < 0.05) {
      if (!any(is.na(weights.aic[1:2])) & round(weights.aic[1]/weights.aic[2],2) >= sig.cutoff) {
        decision <- 'sigmoid'
        ratio <- as.numeric(x[8])
        phase <- as.numeric(x[9])
      } else {
        decision <- 'sinusoidal'
        ratio <- as.numeric(x[14])
        phase <- as.numeric(x[15])
      }
      cur.cyclic <- round(ratio,2) >= ratio.cutoff 
    } else {
      decision <- 'null'
      cur.cyclic <- F
      ratio <- 1.0
      phase <- NA
    }
    c(cur.cyclic, decision, ratio, phase)
  })
  d <- data.frame(t(d))
  for (i in 3:4) { d[,i] <- as.numeric(as.vector(d[,i]))}
  colnames(d) <- c('cyclic', 'decision', 'ratio', 'phase')
  assign(paste('cyclic', 'model', type, sep='.'), d)
}

k.cyclic.model <- merge(cyclic.model.k.rp, cyclic.model.k.tr, by=0, all=TRUE, suffixes=c('.k.rp','.k.tr'))
rownames(k.cyclic.model) <- k.cyclic.model$Row.names
k.cyclic.model <- k.cyclic.model[-1]
k.cyclic.model$cyclic.k.rp <- as.logical(k.cyclic.model$cyclic.k.rp)
k.cyclic.model$cyclic.k.tr <- as.logical(k.cyclic.model$cyclic.k.tr)
no.na.k.cyclic.model <- apply(k.cyclic.model[c(1,5)], c(1,2), function(x) if (is.na(x)) FALSE else x)
na.k.cyclic.model <- data.frame(cyclic.k.rp = k.cyclic.model$cyclic.k.rp, cyclic.k.tr = k.cyclic.model$cyclic.k.tr)
rownames(na.k.cyclic.model) <- rownames(k.cyclic.model)
no.cutoff.k.cyclic.model <- data.frame(cyclic.k.rp = k.cyclic.model$decision.k.rp != 'null', cyclic.k.tr = k.cyclic.model$decision.k.tr != 'null')  ##### ! maybe decision doesnt need .k. !!
rownames(no.cutoff.k.cyclic.model) <- rownames(k.cyclic.model)
k.cyclic.model.sets <- get_sets(no.na.k.cyclic.model)
k.cyclic.model.na.sets <- get_sets(na.k.cyclic.model)
k.cyclic.model.no.cutoff.sets <- get_sets(no.cutoff.k.cyclic.model)

# easy access to gene IDs in sets
for (sett in c('sets', 'na.sets', 'no.cutoff.sets')) {
  csets <- get(paste('k.cyclic', 'model', sett, sep='.'))
  for (lettr in levels(csets)) {
    assign(paste('set',lettr,sett,sep='.'), names(csets)[which(csets == lettr)])
  }
}

#   source('JTK_CYCLE.R')
#   jtkdist(12, as.vector(table(design)))
#   jtk.init(12,2)
#   jtk.efficiency <- get_jtk(eff.noinf)
###### EFFICIENCIES FITTING TO HARMONIC MODEL
#cyclic.eff <- as.data.frame(t(apply(eff.log.noinf, 1, function(x) { get_params(harmonic_regression(x, tpoints, 24), 24, F, T)})))
k.cyclic.eff.fits <- apply(cbind(k.eff.log.noinf, k.eff.weights), 1, function(x) { harmonic_regression(x[1:24], tpoints, 24, x[25:48])})
#cyclic.eff.pvals <- sapply(cyclic.eff.fits, lm.fpval)
#names(cyclic.eff.pvals) <- rownames(eff.log.noinf)
k.cyclic.eff.pvals <- sapply(k.cyclic.eff.fits, function(x) pchisq(2*(logLik(x)-logLik(update(x, .~1))),2,lower.tail = FALSE))
#cyclic.eff.altfdr <- p.adjust(cyclic.eff.altp, method = 'BH')
k.cyclic.eff.fdr <- p.adjust(k.cyclic.eff.pvals, method = 'BH')
k.cyclic.eff.pars <- as.data.frame(t(sapply(k.cyclic.eff.fits, get_params, 24, F, T)))
k.cyclic.eff.decision <- k.cyclic.eff.fdr < 0.05 & k.cyclic.eff.pars$ratio > 1.5


# DAILY PROFILE PLOTS #
for(identifier in names(k.cyclic.model.na.sets)) proplot.kidney(identifier,subfolder='plots',setfolder=TRUE)

#ROSE HISTOGRAMS
ggsave(filename='Rose_diagram_kidney_tr_BD.pdf', rose24h(cyclic.model.k.tr[cyclic.model.k.tr$cyclic=='TRUE','phase'], '#cc6600', 'All circadian at totalRNA level',c(0,150)))
ggsave(filename='Rose_diagram_kidney_tr_setD.pdf', rose24h(cyclic.model.k.tr[k.na.set.d, 'phase'], '#cc6600', 'SetD at totalRNA level',c(0,50)))
ggsave(filename='Rose_diagram_kidney_tr_setB.pdf', rose24h(cyclic.model.k.tr[k.na.set.b, 'phase'], '#cc6600', 'SetB at totalRNA level',c(0,150)))

ggsave(filename='Rose_diagram_kidney_rp_CD.pdf', rose24h(cyclic.model.k.rp[cyclic.model.k.rp$cyclic=='TRUE','phase'], '#0066cc', 'All circadian at translation level',c(0,150)))
ggsave(filename='Rose_diagram_kidney_rp_setD.pdf', rose24h(cyclic.model.k.rp[k.na.set.d, 'phase'], '#0066cc', 'SetD at translation level',c(0,50)))
ggsave(filename='Rose_diagram_kidney_rp_setC_babel_between_66.pdf', rose24h(cyclic.model.k.rp[babel.k.set.c.bet, 'phase'], '#0066cc', 'SetC (babel) at translation level',c(0,10)))

# PHASE Heatmaps kidney
for (set in LETTERS[2:4]) {
  my.heatmap.k(set=set, cyclic.sets = k.cyclic.model.sets)
}


####     BABEL KIDNEY     ####
writeLines('Performing babel analysis, this will take a really long time ...')
k.bab.rna <- k.tr[k.expressed.ids,k.tr.cds]
k.bab.rp <- k.rp[k.expressed.ids,k.rp.cds]
colnames(k.bab.rna) <- paste(design, c('1','2'), sep='.')
colnames(k.bab.rp) <- paste(design, c('1','2'), sep='.') # save(list = c('k.bab.rna', 'k.bab.rp'), file = 'import_babel.RData'). Open R on the server and load(list = c('k.bab.rna', 'k.bab.rp'), file = 'import_babel.RData')
## run on server:
options(mc.cores = 32)
design <- rep(paste('ZT', formatC(seq(0,22,2), width=2, format='d', flag='0')
library(babel)
k.bab <- babel(k.bab.rna, k.bab.rp, group=design, nreps=10000000)
## save k.bab and load it back on PC
k.sels <- list()
k.res <- matrix(nrow=12, ncol=2)
for (i in 1:12) {
  k.sel <- which(k.bab$combined[[i]]['FDR'] < 0.05)
  k.sels[[i]] <- k.sel
  direction <- summary(as.factor(k.bab$combined[[i]][k.sel,'Direction'])) # ..combined[[1]] ---> ..combined[[i]] ((??))
  k.res[i,] <- direction
}
common <- k.sels[[1]]
for (i in 2:12) {
  common <- intersect(common, k.sels[[i]])
}
k.res.common <- matrix(c(69,49), nrow=12, ncol=2, byrow=T)
k.res.unique <- k.res - k.res.common
k.res.all <- cbind(k.res.common[,1], k.res.unique[,1], k.res.common[,2], k.res.unique[,2])

## Hierarchical clustering
k.sels.comp <- matrix(ncol=12,nrow=12)
for (i in 1:12) {
  for (j in i:12) {
    unionset <- union(k.sels[[i]], k.sels[[j]])
    commonset <- intersect(k.sels[[i]], k.sels[[j]])
    k.sels.comp[j,i] <- length(commonset)/length(unionset)
  }
}
k.sels.fit <- hclust(as.dist(1-k.sels.comp), method='ward')

pdf('babel.within.occupancy.clustering_kidney.pdf', useDingbats=F)
layout(c(1,1))
par(mar=c(5,4,4,2)+0.1)
plot(k.sels.fit, labels=paste('ZT',formatC(seq(0,22,2), width=2, format='d', flag='0'), sep=''), xlab='Jaccard distance', main='Hierarchical clustering of time-points based on similarity between\nsets of genes with significant alterations in ribosome occupancy',sub='')
dev.off()

## Within occupancy histogram plot
pdf('babel.within.occupancy.histogram_kidney.pdf', useDingbats=F)
myax <- seq(0,5000,100)
layout(matrix(c(1,2),1,2))
par(mar=c(2,2,2,0))
order_as_clustering <- TRUE
if (order_as_clustering) { # ORDER SAME AS IN CLUSTERING
  barplot(t(-k.res.all[rev(k.sels.fit$order),1:2]), horiz=T, axes=F, xlim=c(-500,0))
  axis(side=1, at=-myax, labels=myax)
  par(mar=c(2,0.5,2,1.5))
  barplot(t(k.res.all[rev(k.sels.fit$order),3:4]), horiz=T, axes=F, names.arg=paste('ZT',seq(0,22,2)[rev(k.sels.fit$order)],sep=''), axis.lty=0, las=1, xlim=c(0,500), xlab='High ribosome occupancy')
  axis(side=1, at=myax, labels=myax)
} else { # ORDER ZT0 -> ZT22
  barplot(t(-k.res.all[12:1,1:2]), horiz=T, axes=F, xlim=c(-500,0))
  axis(side=1, at=-myax, labels=myax)
  par(mar=c(2,0.5,2,1.5))
  barplot(t(k.res.all[12:1,3:4]), horiz=T, axes=F, names.arg=paste('ZT',seq(22,0,-2),sep=''), axis.lty=0, las=1, xlim=c(0,500), xlab='High ribosome occupancy')
  axis(side=1, at=myax, labels=myax)
}
dev.off()
layout(c(1,1))
par(mar=c(5,4,4,2)+0.1)

k.within.sig <- c()
for (j in 1:12) k.within.sig <- union(k.within.sig, k.sels[[j]])
k.within.sig.ids <- k.expressed.ids[k.within.sig]

# Betwwen time-points analysis (TE changes over timepoints)
k.between.sig <- c()
for (i in 1:length(k.bab$between)) {
  k.between.sig <- union(k.between.sig, which(k.bab$between[[i]]['FDR'] < 0.05))
}
k.between.sig.ids <- k.expressed.ids[k.between.sig]
k.all.sig <- union(k.within.sig, k.between.sig)
k.all.sig.ids <- union(k.within.sig.ids, k.between.sig.ids)

# Heatmaps
for (type in c('k.within', 'k.between', 'k.all')) {
  for (set in LETTERS[2:4]) {
    my.heatmap.k(set=set, cyclic.sets=k.cyclic.model.sets[get(paste(type,'sig','ids', sep='.'))], extra_label=paste('sig', type, sep='_'))
  }
}

# PCA analysis and 3D plot
log2.cds.k.tr <- data.frame(t(apply(norm.cds.k.tr,1,function(x) log2(x/mean(x)))))
log2.cds.k.rp <- data.frame(t(apply(norm.cds.k.rp,1,function(x) log2(x/mean(x)))))
log2.k.diff <- log2.cds.k.rp[k.expressed.ids,] - log2.cds.k.tr[k.expressed.ids,]
log2.k.diff.noinf <- no.inf(log2.k.diff)
log2.k.mean.diff <- data.frame(t(apply(log2.k.diff.noinf, 1, function(x) colMeans(rbind(x[seq(1,23,2)], x[seq(2,24,2)]), na.rm=T))))
pca.k <- princomp(log2.k.mean.diff[k.all.sig.ids,])

library('scatterplot3d')
s3d <- scatterplot3d(pca.k$loadings[,1:3], pch=21, bg=c(rep('white',6), rep('black',6)), type='p', grid=F, box=T, xlim=c(-0.5, 0.5), ylim=c(-0.6,0.6), zlim=c(-0.6,0.6))
s3d$plane3d(0,0,0, lty='dotted', lwd=0.7)
tpos <- s3d$xyz.convert(pca.k$loadings[,1:3])
pred2d <- s3d$xyz.convert(pca.k$loadings[,1], pca.k$loadings[,2],rep(0, nrow(pca.k$loadings)))
segments(tpos$x, tpos$y, pred2d$x, pred2d$y, lty=4)
text(tpos, paste('ZT',formatC(seq(0,22,2), width=2, format='d', flag='0'), sep=''), xpd=T, pos=2)

# Tables
k.bab.direction <- k.bab$combined[[1]]$Direction[k.all.sig]
k.bab.agree <- rep(TRUE, length(k.bab.direction))
for (i in 1:12) {
  k.bab.agree <- k.bab.agree & k.bab$combined[[i]]$Direction[k.all.sig] == k.bab.direction
}
k.bab.direction[!k.bab.agree] <- 0
k.bab.direction <- as.factor(k.bab.direction)
levels(k.bab.direction) <- c('under', 'mixed', 'over')
babel.k.cyclic.sets <- data.frame(set=k.cyclic.model.sets[k.all.sig.ids], between=k.all.sig.ids %in% k.between.sig.ids, within=k.all.sig.ids %in% k.within.sig.ids, direction=k.bab.direction)
write.csv(babel.k.cyclic.sets, file='babel.cyclic.sets_kidney.csv')

