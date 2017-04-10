##### DETECTION OF TRANSLATED UORFS  IN KIDNEY ###
uorf.cols <-c('gid','trid', 'start', 'stop', 'cov', 'frame0', 'frame1', 'frame2') # cov (coverage) is how many positions have reads (at least one read)

uorf.k.rp <- read.table('uorfs_kidney_rp2.txt', header=F, col.names=uorf.cols, stringsAsFactors = F)
uorf.k.tr <- read.table('uorfs_kidney_tr2.txt', header=F, col.names=uorf.cols, stringsAsFactors = F)

uorf.k.tr <- cbind(uorf.k.tr, cov.f = uorf.k.tr$cov/(uorf.k.tr$stop - uorf.k.tr$start))
uorf.k.tr <- cbind(uorf.k.tr, count = rowSums(uorf.k.tr[,6:8]))
uorf.k.tr.pval <- apply(uorf.k.tr, 1, function(x) {
  obs <- as.numeric(x[6:8])
  p.val <- if (sum(obs) > 0) chisq.test(obs)$p.value else NA
  p.val})
uorf.k.tr <- cbind(uorf.k.tr, pval=uorf.k.tr.pval, fdr=p.adjust(uorf.k.tr.pval,'fdr'))
uorf.k.tr <- cbind(uorf.k.tr, frame=uorf.k.tr$start %% 3)
uorf.k.tr <- cbind(uorf.k.tr, frame.obs = unlist(apply(uorf.k.tr[,6:8],1,function(x) {
  m <- which(x == max(x))
  if (length(m) > 1) m <- NA
  m-1})))

uorf.k.rp <- cbind(uorf.k.rp, cov.f = uorf.k.rp$cov/(uorf.k.rp$stop -uorf.k.rp$start))
uorf.k.rp <- cbind(uorf.k.rp, count = rowSums(uorf.k.rp[,6:8]))
uorf.k.rp.pval <- apply(uorf.k.rp, 1, function(x) {
  obs <- as.numeric(x[6:8])
  p.val <- if (sum(obs) > 0) chisq.test(obs)$p.value else NA
  p.val})
uorf.k.rp <- cbind(uorf.k.rp, pval=uorf.k.rp.pval, fdr=p.adjust(uorf.k.rp.pval,'fdr'))
uorf.k.rp <- cbind(uorf.k.rp, frame=uorf.k.rp$start %% 3)
uorf.k.rp <- cbind(uorf.k.rp, frame.obs = unlist(apply(uorf.k.rp[,6:8],1,function(x) {
  m <- which(x == max(x))
  if (length(m) > 1) m <- NA
  m-1})))

uorf.filter.k <- intersect(which(uorf.k.rp$cov.f > 0.1 & uorf.k.rp$fdr < 0.05 & uorf.k.rp$frame == uorf.k.rp$frame.obs), 
                           which(uorf.k.tr$cov.f > 0.1 & uorf.k.tr$fdr < 0.05 & uorf.k.tr$frame == uorf.k.tr$frame.obs))
k.uorf.translated <- uorf.k.rp$gid[setdiff(which(uorf.k.rp$cov.f > 0.1 & uorf.k.rp$fdr < 0.05 & uorf.k.rp$frame == uorf.k.rp$frame.obs), uorf.filter.k)]
k.uorf.translated = intersect(k.uorf.translated,k.expressed.ids) #1593
k.uorf.translated_trid <- unique(uorf.k.rp$trid[which(uorf.k.rp$gid %in% k.uorf.translated)])
