# DEFINE the KIDNEY  COUNT SAMPLES to be used (files containing counts), and composite sizes.
# same is done for liver.
k.rp.samples <- read.table('rp_count_files_vFinal.txt', stringsAsFactors=F, col.names=c('sample','library','rep','path','name','ext','use'), header=F)
k.rp.samples <- k.rp.samples[order(k.rp.samples$sample),]
k.tr.samples <- read.table('tr_count_files_vFinal.txt', stringsAsFactors=F, col.names=c('sample','library','rep','path','name','ext','use'), header=F)
k.tr.samples <- k.tr.samples[order(k.tr.samples$sample),]
#k.annot is the file containing the expressed transcripts in kidney with info about cds start/end, etc; output of prepare_gtf_v3.py Size is total length(5'+cds+3')
k.annot <- read.table('prepared_cds_v38.75_filtered_kidney_vFinal.txt', header=F, sep='\t', stringsAsFactors=F, 
                      col.names=c('geneID', 'status', 'trID', 'tr.type', 'protID', 'size', 'cds.start', 'cds.end'))
k.composite <- k.annot[k.annot$tr.type == 'composite',] #composite contains only the combined info of transcripts belonging to one gene. It contains combined (composite) size.
k.composite.sizes <- with(k.composite, cbind(x5utr=cds.start, cds=cds.end-cds.start, x3utr=size-cds.end, utr=NA, not.in.db=NA))
row.names(k.composite.sizes) <- k.composite$geneID 
k.annot$status <- factor(k.annot$status)
k.annot$tr.type <- factor(k.annot$tr.type)
# Size of raw count files (= nº genes) 13769 as of Feb2015; 13830(MArch2015,final)
k.tsize  <- dim(k.composite)[1]

# READ COUNT FILES FOR KIDNEY. Output is 'k.rp' and 'k.tr' object, containing number of reads (in 5',cds,3') per gene(rows) per timepoint(columns)
# same is done for liver.
for (dset in c('k.rp', 'k.tr')) {
  # Read counts of selected samples
  samples <- subset(get(paste(dset, 'samples', sep='.')), use=='Y')
  counts <- read_counts(samples, 'all', k.tsize) 
  # Extract region info
  regions <- unlist(lapply(strsplit(colnames(counts), split='.', fixed=T),function(x) paste(x[2:length(x)], collapse='.')))
  # Merge separate runs of same sample
  ttr <- t(counts)
  ttr.agr <- aggregate(x=ttr, by=list(rep(paste(samples$library, samples$rep, sep='_'), each=5), regions), sum)
  ttr.names <- apply(ttr.agr, 1, function(x) paste(x[1:2], collapse='.'))
  counts <- t(ttr.agr[,3:ncol(ttr.agr)])
  colnames(counts) <- ttr.names
  assign(paste(dset), counts)
  # Split data according to regions ?? necessary ??
  for (region in c('cds', 'x5utr', 'x3utr')) {
    assign(paste(dset, region, sep='.'), grep(region, colnames(counts)))
  }
  rm(ttr, counts, samples, ttr.agr, regions)
}

# Low-level filtering based on raw counts.CDS only based filtering 
#Filters a row of counts (gene counts) by min number(min.fraction) of columns with a min count (min.count;default=10). So we keep genes with at least 10 counts (on the cds) in at least 1/4 of samples (6 samples with at least 10 counts)
k.rp_filter <- get_count_filter(k.rp[,k.rp.cds], min.fraction=1/4) 
k.tr_filter <- get_count_filter(k.tr[,k.tr.cds], min.fraction=1/4)
k.rp_filter.ids <- rownames(k.rp)[k.rp_filter] 
k.tr_filter.ids <- rownames(k.tr)[k.tr_filter] 
k.expressed.ids <- intersect(k.rp_filter.ids, k.tr_filter.ids) #pass the filter in RP and TR in kidney.  12423

# Upper-Quartile based normalization
# Based on only CDS counts
design <- rep(paste('ZT', formatC(seq(0,22,2), width=2, format='d', flag='0'), sep=''), each=2)

norm.factors.k.rp <- get_norm_factors(k.rp[k.rp_filter, k.rp.cds], design)
norm.factors.k.tr <- get_norm_factors(k.tr[k.tr_filter, k.tr.cds], design)
norm.factors.l.rp <- get_norm_factors(l.rp[l.rp_filter, l.rp.cds], design)
norm.factors.l.tr <- get_norm_factors(l.tr[l.tr_filter, l.tr.cds], design)
norm.cds.k.rp <- t(apply(k.rp[k.rp_filter,k.rp.cds],1,'/',norm.factors.k.rp))
norm.cds.k.tr <- t(apply(k.tr[k.tr_filter,k.tr.cds],1,'/',norm.factors.k.tr))
norm.cds.l.rp <- t(apply(l.rp[l.rp_filter,l.rp.cds],1,'/',norm.factors.l.rp))
norm.cds.l.tr <- t(apply(l.tr[l.tr_filter,l.tr.cds],1,'/',norm.factors.l.tr))

norm.cds.k.rp.mean <- 0.5 * (norm.cds.k.rp[,seq(1,23,2)] + norm.cds.k.rp[,seq(2,24,2)])
norm.cds.k.tr.mean <- 0.5 * (norm.cds.k.tr[,seq(1,23,2)] + norm.cds.k.tr[,seq(2,24,2)])
norm.cds.l.rp.mean <- 0.5 * (norm.cds.l.rp[,seq(1,23,2)] + norm.cds.l.rp[,seq(2,24,2)])
norm.cds.l.tr.mean <- 0.5 * (norm.cds.l.tr[,seq(1,23,2)] + norm.cds.l.tr[,seq(2,24,2)])

#RPKM #reads / (length ⨯ total #reads in sample)
writeLines('   ... calculating RPKM values')
libsize.k.rp <- get_norm_factors(k.rp[k.rp_filter, k.rp.cds], design, 'efflibsize') 
libsize.k.tr <- get_norm_factors(k.tr[k.tr_filter, k.tr.cds], design, 'efflibsize') 
libsize.l.rp <- get_norm_factors(l.rp[l.rp_filter, l.rp.cds], design, 'efflibsize') 
libsize.l.tr <- get_norm_factors(l.tr[l.tr_filter, l.tr.cds], design, 'efflibsize') 
rpkm.cds.k.rp <- norm.cds.k.rp / k.composite.sizes[k.rp_filter.ids, 'cds'] / libsize.k.rp * 1e+09 
rpkm.cds.k.tr <- norm.cds.k.tr / k.composite.sizes[k.tr_filter.ids, 'cds'] / libsize.k.tr * 1e+09
rpkm.cds.l.rp <- norm.cds.l.rp / l.composite.sizes[l.rp_filter.ids, 'cds'] / libsize.l.rp * 1e+09 
rpkm.cds.l.tr <- norm.cds.l.tr / l.composite.sizes[l.tr_filter.ids, 'cds'] / libsize.l.tr * 1e+09
rpkm.cds.k.rp.mean <- 0.5 * (rpkm.cds.k.rp[,seq(1,23,2)] + rpkm.cds.k.rp[,seq(2,24,2)])
rpkm.cds.k.tr.mean <- 0.5 * (rpkm.cds.k.tr[,seq(1,23,2)] + rpkm.cds.k.tr[,seq(2,24,2)])
rpkm.cds.l.rp.mean <- 0.5 * (rpkm.cds.l.rp[,seq(1,23,2)] + rpkm.cds.l.rp[,seq(2,24,2)])
rpkm.cds.l.tr.mean <- 0.5 * (rpkm.cds.l.tr[,seq(1,23,2)] + rpkm.cds.l.tr[,seq(2,24,2)])

#TRANSLATION EFFICIENCIES
#kidney (same is done for liver)
k.eff <- rpkm.cds.k.rp[k.expressed.ids,] / rpkm.cds.k.tr[k.expressed.ids,]
k.eff.noinf <- as.data.frame(no.inf(k.eff))
k.eff.log <- log2(k.eff) #log efficiencies
k.eff.log.noinf <- no.inf(k.eff.log)
k.eff.log.noinf.reps <- array(c(k.eff.log.noinf[,seq(1,23,2)], k.eff.log.noinf[,seq(2,24,2)]), dim=c(nrow(k.eff.log.noinf), 12, 2))
k.eff.log.noinf.mean <- apply(k.eff.log.noinf.reps, c(1,2), mean, na.rm=T) #average efficiency of the two replicates
rownames(k.eff.log.noinf.mean) <- rownames(k.eff.log.noinf)
k.eff.weights <- sqrt(rpkm.cds.k.rp[k.expressed.ids,] * rpkm.cds.k.tr[k.expressed.ids,]) ## weights based on magnitude -- sqrt(TR * RPF). used when fitting efficiencies to harmonic regression in circ.analysis
