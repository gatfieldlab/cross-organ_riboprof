# Functions with no dependency

### Required packages -- Please all packages here
my_packages = c('edgeR','DESeq', 'ggplot2', 'MASS', 'minpack.lm', 'VennDiagram', 'gtable', 'gdata',
                'gtools', 'RColorBrewer', 'babel', 'ash')
dependency <- unlist(lapply(my_packages, require, character.only=T))
dependency <- data.frame(is.ok = dependency, row.names=my_packages)

#source('/home/barpat/Projects/riboprof/analysis/motion_graphs/motion_graphs.R')

# Filters a row of counts by min number(fraction) of columns with a min count
filter_by_count_number <- function (one.row, min.count=10, min.fraction=0.25) {
length(which(one.row > min.count)) / length(one.row) > min.fraction
}

get_count_filter <- function (countTable, min.count=10, min.fraction=1/6, logical=FALSE, reverse=FALSE) {
  f <- apply(countTable, 1, filter_by_count_number, min.count, min.fraction)
  if (reverse) { f <- !f }
  if (logical) { return(f) }
  else { return(which(f)) }
}

no.inf <- function(d.f, repl=NaN) {
  apply(d.f, c(1,2), function(x) if (is.infinite(x)) repl else x )
}

"%ni%" <- Negate("%in%")

# Low-level function to read from count text file
read_count_table_v3 <- function (filename, nrows) {
    read.table(textConnection(readLines(filename)[1:nrows]), stringsAsFactors=F, sep='\t', col.names=c('GeneID','x5utr','cds','x3utr','utr','not.in.db'))
}

# High-level function to make a data frame from supplied samples and fields (cds, utr etc)
read_counts <- function (samples, type=c('cds'), nrows) {
    filenames <- with(samples, paste(path, name, ext, sep=''))
    first_tb <- read_count_table_v3(filenames[1], nrows)
    if (sum(c('all','All','a','A') %in% type)) {
        read_cols = which(!colnames(first_tb) %in% 'GeneID')
    } else {
        read_cols = which(colnames(first_tb) %in% setdiff(type, 'GeneID'))
    }
    datafr <- data.frame(first_tb[c(1, read_cols)])
    for (i in 2:nrow(samples)) {
        datafr <- cbind(datafr, read_count_table_v3(filenames[i], nrows)[read_cols])
    }
    row.names(datafr) <- datafr$GeneID
    datafr <- datafr[-1]
    colnames(datafr) <- paste(rep(samples$sample,each=length(read_cols)), colnames(first_tb)[read_cols], sep='.')
    datafr
}


get_read_weights <- function(trim_data, sizes) {
  tdata <- cbind(length=1:98, count=0)
  for (i in which(trim_data$length==1):nrow(trim_data)) { tdata[trim_data$length[i],'count']=trim_data$count[i]}
  b1 <- list(nc=tdata[,'count'], ab=c(0.5,98.5), nskip=c(0))
  f1 <- ash1(b1,3)
  f1sum <- sum(f1$y[sizes])
  f1$y[sizes]/f1sum
}

strip <- function(str) {
  gsub("^\\s+|\\s+$", "", str)
}

trimean <- function(x) {
  q <- as.numeric(quantile(x, c(0.25, 0.50, 0.75)))
  (q[1] + 2*q[2] + q[3]) / 4
}

get_par_limits <- function(p=4) {
  d <- (100 + 2*p) / p
  usr <- par('usr')
  xr <- (usr[2] - usr[1]) / d
  yr <- (usr[4] - usr[3]) / d
  xlim <- c(usr[1] + xr, usr[2] - xr)
  ylim <- c(usr[3] + yr, usr[4] - yr)
  list(xlim=xlim, ylim=ylim)
}

plot_reps <- function(norm.dataset, design, label) {
  for (i in seq(1,23,2)) {
    corr = cor(norm.dataset[,i], norm.dataset[,i+1], use='complete', method='spear')
    xlab=paste(label, design[i], 'rep1', sep=', ')
    ylab=paste(label, design[i+1], 'rep2', sep=', ')
    pdf(paste(label, design[i], 'rep_scatter', 'pdf', sep='.'), useDingbats=F)
    log_scatter(norm.dataset[,i], norm.dataset[,i+1], pch=20, cex=0.8, xlab=xlab, ylab=ylab)
    abline(0,1,col='red')
    legend('topleft', legend=paste('Spearman correlation=',round(corr,digits = 2)))
    dev.off()
  }
}


if (dependency['DESeq', 'is.ok']) {
  get_cds <- function (counts, design, normfactors, ...) {
      cds <- newCountDataSet(counts, design)
      sizeFactors(cds) <- normfactors
      cds <- estimateDispersions(cds, ...)
      cds
  }
  get_dispersions <- function(cds) {
    1/fitInfo(cds)$fittedDispEsts
  }
} else {
  get_cds <- function (...) {
    stop("Original 'get_cds' function failed due to dependency problems.")
  }
  get_dispersions <- function (...) {
    stop("Original 'get_diepersions' function failed due to dependency problems.")
  }
}

# if (dependency['DESeq2', 'is.ok']) {
#   get_deseq2_matrix <- function(counts,colData,design,normfactors,...) {
#     deseq2_matrix <- DESeqDataSetFromMatrix(counts,colData=colData,design=design)
#     sizeFactors(deseq2_matrix) <- normfactors
#     deseq2_matrix <- estimateDispersions(deseq2_matrix,...)
#     deseq2_matrix
#   }
# } else {
#   get_deseq2_matrix <- function(...){
#     stop("Original 'get_cds' function failed due to dependency problems.")
#   }
# }
    

# Normalization functions
# Dependance: edgeR


if (dependency['edgeR', 'is.ok']) {
  get_norm_factors <- function(counts, group, what='normfactor', method="upperquartile", ...) {
    dge_counts <- DGEList(counts=counts, group=group)
    dge_counts <- calcNormFactors(dge_counts, method=method, ...)
    eff.libsize <- dge_counts$samples$lib.size * dge_counts$samples$norm.factors
    eff.libgeo <- exp(mean(log(eff.libsize)))
    if (what == 'normfactor') return(eff.libsize / eff.libgeo)
    if (what == 'efflibsize') return(eff.libgeo)
    if (what == 'dgenorm') return(dge_counts$samples$norm.factors)
  }
} else {
  get_norm_factors <- function (...) {
    stop("Original 'get_norm_factors' function failed due to dependency problems.")
  }
}

### Circadian analysis functions

# log-likelihood and NB harmonic regression

get_permarray <- function() {
  permarray <- matrix(0,64,12)
  permarray[1,] <- seq(1,12)
  for (i in 1:6) {
    n <- 2^(i-1)
    for (j in 1:n) {
      permline <- permarray[j,]
      old <- permline[i]
      permline[i] <- permline[i+6]
      permline[i+6] <- old
      permarray[n+j,] <- permline
    }
  }
  permarray
}

PERMARRAY <- get_permarray()

get_perm <- function (raw.data) {
  perm <- matrix(0, 64, 12)
  for (j in 1:64) {
    perm[j,] <- raw.data[PERMARRAY[j,]]
  }
  perm
}

estimate_mode <- function(x) {
  d <- density(x, na.rm=TRUE)
  d$x[which.max(d$y)]
}

get_perm_ar_periods <- function (raw.data, delta, ar.methods,
                                 exp.period, error) {
  limits = c(exp.period - error, exp.period + error)
  perm <- get_perm(as.numeric(raw.data))
  results <- NULL
  for (method in ar.methods) {
    periods <- NULL
    for (j in 1:64) {
      ar.model <- tryCatch(spec.ar(perm[j,], plot=FALSE, method=method,
                                   order= exp.period/delta), error = function(err) {
                                     print("can't fit ar model"); return(NULL) } )
      if (!is.null(ar.model)) {
        peaks <- which(diff(sign(diff(ar.model$spec)))==-2)+1
        c.periods <- 1/ar.model$freq[peaks] * delta
        periods <- c(periods, c.periods[c.periods > limits[1] & c.periods < limits[2]])
      }
    }
    if (!is.null(periods)) {
      len <- sum(!is.na(periods))
      if (len < 6) {
        est_period <- mean(periods, na.rm=TRUE)
      } else if (len < 12) {
        est_period <- mean(periods, trim=0.1, na.rm=TRUE)
      } else {
        est_period <- estimate_mode(periods)
      }
    } else {
      est_period = NaN
    }
    results <- c(results, est_period)
  }
  results
}

harmonic_regression.nb <- function (raw.data, tpoints, period, theta, link='log') {
  cosx <- cos(2*pi/period * tpoints)
  sinx <- sin(2*pi/period * tpoints)
  tryCatch(glm(raw.data ~ cosx + sinx, family=negative.binomial(theta, link=link)),
           error = function (err) {writeLines(paste('Harmonic regression failed:',err)); return(NULL)})
}

harmonic_regression <- function (raw.data, tpoints, period, weights=NULL) {
  cosx <- cos(2*pi/period * tpoints)
  sinx <- sin(2*pi/period * tpoints)
  tryCatch(lm(raw.data ~ cosx + sinx, weights=weights),
           error = function (err) {writeLines(paste('Linear harmonic regression failed:',err)); return(NULL)})
}
get_params <- function (fit, period, no.negative=T, logval=F, logbase=2)  {
  if (!is.null(fit) && !any(is.na(coef(fit)))) {
    # Calculate parameters amp, phase
    coefs <- coef(fit)
    amp <- as.numeric(sqrt(coefs[2]^2 + coefs[3]^2))
    phase <- as.numeric(atan2(coefs[3], coefs[2]) / (2*pi) * period)
    if (phase < 0) {
      phase <- phase + period
    }
    intercept <- coefs[[1]]
    minp <- intercept - amp
    if (no.negative) {
      minl <- max(mean(sort(fit$model$count)[1:2], na.rm=T), 0.1)
      minp <- max(minp, minl)
    }
    maxp <- intercept + amp
    if (isTRUE(logval)) {
      ratio <- logbase^(maxp - minp)
    } else {
      ratio <- maxp / minp
    }
    c(mean=intercept, amp=amp, min=minp, max=maxp, ratio=ratio, phase=phase)
  } else {
    c(mean=NA, amp=NA, min=NA, max=NA, ratio=NA, phase=NA)
  }
}

get_params.loglink <- function (glmfit, period)  {
  # Calculate parameters amp, phase from log-linked glm
  coefs <- coef(glmfit)
  amp <- as.numeric(sqrt(coefs[2]^2 + coefs[3]^2)) * exp(coefs[1])
  phase <- as.numeric(atan2(coefs[3], coefs[2]) / (2*pi) * period)
  if (phase < 0) {
    phase <- phase + period
  }
  cbind(amp=amp, phase=phase)
}

get_ar <- function(raw.data, tpoints, ar.methods=c('yw','mle','burg'),
                   perm = FALSE,
                   exp.period=24, theta=50, link='log', error=2) {
  if (max(raw.data) < 5) {
    results <- cbind(period=NaN, amp=NaN, phase=NaN, p.val=NaN)
  } else {
    delta <- tpoints[2] - tpoints[1]
    if (perm) {
      ar.periods <- get_perm_ar_periods(raw.data, delta=delta, ar.methods=ar.methods,
                                        exp.period=exp.period, error=error)
    } else {
      ar.periods <- NA
    }
    if (perm & sum(!is.na(ar.periods)) > 0) {
      periods <- ar.periods[!is.na(ar.periods)]
      if (!length(which(periods==exp.period))) {
        periods <- c(periods, exp.period)
      }
    } else {
      periods <- exp.period
    }
    best_fit <- NULL
    for (period in periods) {
      harmonic_fit <- harmonic_regression.nb(round(raw.data), tpoints, period, theta, link=link)
      if (!is.null(harmonic_fit) && (is.null(best_fit) || (harmonic_fit$aic < best_fit$aic))) {
        best_fit <- harmonic_fit
        best_period <- period
      }
    }
    if (!is.null(best_fit)) {
      null_fit <- update(best_fit, . ~ 1)
      pval <- pchisq(2*(logLik(best_fit) - logLik(null_fit)), 2, lower.tail=FALSE)
      results <- cbind(period=best_period, get_params.loglink(best_fit, best_period), p.val=pval)
    } else {
      results <- cbind(period=NaN, amp=NaN, phase=NaN, p.val=NaN)
    }
  }
  results
}

# JTK
get_jtk <- function (dataset) {
  res <- apply(dataset,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  #bhq <- p.adjust(unlist(res[,1]),"BH")
  #res <- cbind(bhq,res)
  colnames(res) <- c("pval","period","phase","amp")
  res
}

# Sigmoid non-linear curve fitting

SD <- function(i,m,p,time) 1/(1+exp(8/m*(time+24-i*24-(p%%24)-m/2)))
SU <- function(i,l,p,time) 1/(1+exp(-8/l*(time+24-i*24-(p%%24)+l/2)))

# helper function to return residuals for low level functions of minpack.lm
sigmoid.res <- function(pars, observed, time) observed - sigmoid.i(pars, time)

# helper function to extract datapoints from a sigmoid fit
sigmoid.i <- function (pars, time) {
  with (as.list(c(pars)), {
    SU <- function(i) {1/(1+exp(-8/l*(time+24-i*24-(p%%24)+l/2)))}
    SD <- function(i) {1/(1+exp(8/m*(time+24-i*24-(p%%24)-m/2)))}
    return(B*(1+(fc-1)*(SU(0)+SD(0)+SU(1)+SD(1)+SU(2)+SD(2)-3)))
  })
}

get_aic <- function(fit) {
  w <- fit$weights
  if (is.null(w)) w <- 1
  if (is.null(fit$residuals)) {
    r <- resid(fit)
  } else {
    r <- fit$residuals
  }
  rss <- sum( w*r^2)
  n <- 24 # hardcoded !!!
  k <- length(coef(fit))
  n * log((2 * pi)/n) + n + 2 + n * log(rss) + 2 * k
}

get_ss <- function(fit) {
  if (is.null(fit) || is.na(fit)) {
    return(NA)
  }
  w <- fit$weights
  if (is.null(w)) w <- 1
  if (is.null(fit$residuals)) {
    r <- resid(fit)
  } else {
    r <- fit$residuals
  }
  sum(w * r^2)
}
delta.aic.all <- function(fits) {
  d1 <- delta.aic(fits$sigmoid,fits$null)
  d2 <- delta.aic(fits$sinusodial, fits$null)
  d3 <- 0
  min.d <- min(d1, d2, d3, na.rm=T)
  c(d1-min.d, d2-min.d, d3-min.d)
}
delta.aic.all2 <- function(fits) {
  d1 <- delta.aic(fits$sigmoid,fits$null)
  d2 <- delta.aic(fits$sinusodial, fits$null)
  d3 <- 0
  min.d <- min(d1, d2, d3, na.rm=T)
  c(d1-min.d, d2-min.d, d3-min.d, d3-min.d)
}
delta.aic <- function(fit2, fit1) {
  if (is.null(fit1) | is.null(fit2)) {
    return(NA)
  }
  k1 <- length(coef(fit1)) + 1
  k2 <- length(coef(fit2)) + 1
  N <- 24
  N*log(get_ss(fit2)/get_ss(fit1)) + 2*(k2-k1) - 2*k1*(k1+1)/(N-k1-1) + 2*k2*(k2+1)/(N-k2-1)
}
waic <- function (aic.deltas) {
  w.sum <- sum(exp(-0.5*aic.deltas), na.rm=T)
  exp(-0.5*aic.deltas)/w.sum
}
aic.pval <- function (fit2, fit1) {
  if (is.null(fit2) | is.null(fit1)) return(NA)
  delta.aicc <- delta.aic(fit2, fit1)
  evidence <- exp(0.5*delta.aicc)
  evidence / (1 + evidence)
}

point_harmonic <- function(harmonic_fit) {
  pars <- get_params(harmonic_fit,24)
  points(seq(0,22,0.5), pars[1]+pars[2]*cos(pi/12*(seq(0,22,0.5) - pars[3])), col='green', type='l', lty=2)
}

fitModel <- function(raw.data, tpoints, period=24, theta=50, weights=NULL) {
  # 0 counts produces error in nls call - fix by setting them to 0.1 if all reps are 0!
  DataT <- data.frame(time=tpoints, cosx=cos(2*pi/period * tpoints), sinx=sin(2*pi/period * tpoints), count=as.numeric(raw.data))
  sumd <- aggregate(count ~ time, data=DataT, sum)
  zerod <- subset(sumd, count==0, time, drop=T)
  if (length(zerod)) {
    DataT[DataT$time %in% zerod,]$count <- min(1, min(sumd[sumd$time %ni% zerod,]$count/2))
  }
  ###
  DataT <- cbind(DataT, logcount = log(DataT$count))
  ###
  sinusodial <- fitSinusodial(DataT, period, theta, weights)
  sigmoid <- fitSigmoid(DataT, weights=sinusodial$sinusodial$weights)
  list(sinusodial=sinusodial$sinusodial, sigmoid=sigmoid, null=sinusodial$null, weight_mode=sinusodial$weight_mode)
}

fitSinusodial <- function(DataT, period=24, theta=50, weights=NULL) {
  if (is.null(weights)) {
    nb_fit <- tryCatch(glm(count ~ cosx + sinx, data=DataT, family=negative.binomial(theta, link='identity')),
                       error = function (err) {writeLines(paste('Negative binomial harmonic regression failed:',err, 'Using Poisson weights...')); return(NULL)})
    if (is.null(nb_fit)) {
      weights <- 1/merge(DataT, aggregate(count ~ time, DataT, mean), by=1)$count.y
      weight_mode <- 'poisson'
    } else {
      weights <- nb_fit$weights
      weight_mode <- 'negbinom'
    }
  } else {
    weight_mode <- 'manual'
  }
  lm_fit <- tryCatch(lm(count ~ cosx + sinx, data=DataT, weights=weights),
                       error = function (err) {writeLines(paste('Weighted least square harmonic regression failed:',err)); return(NULL)})
  if (is.null(lm_fit)) {
    null_fit <- lm(count ~ 1, data=DataT, weights=weights)
  } else {
    null_fit <- update(lm_fit, . ~ 1)
  }
  return(list(sinusodial=lm_fit, null=null_fit, weight_mode=weight_mode))
}

# Sigmoid fitting function without restriction up + down < 24 hrs

if (dependency['minpack.lm', 'is.ok']) {
  fitSigmoid <- function (DataT, theta=50, weights=NULL, plot.it=FALSE, debug=FALSE,
                          lower=c(fc=1,p=-12,l=3,m=3,B=5), upper=c(fc=1000,p=36,l=16,m=16,B=800000)) {
    getFit <- function (DataT, start, lower, upper, weights) {
      ### count -> logcount
      tryCatch(nlsLM(count ~ B*(1+(fc-1)*(SU(0,l,p,time)+SD(0,m,p,time)+
                             SU(1,l,p,time)+SD(1,m,p,time)+SU(2,l,p,time)+SD(2,m,p,time)-3)),
                      DataT, start=start, control=c(maxiter=100), lower=lower, upper=upper, weights=weights),
               error=function(err) { writeLines(paste('Sigmoid fit failed:', strip(err))); return(NULL)})
    }
      if (!any(is.na(DataT$count))) {
        Mu <- aggregate(count ~ time, data=DataT, mean)
        mu <- Mu$count
        # var based on negative binomial
        var <- mu + mu^2/theta
        # alternative var -direct estimation
        var2 <- tapply(DataT$count, DataT$time, var)
        # var based on poisson
        var3 <- mu
        # initial param estimations
        max_mu <- max(mu)
        B_e <- max(min(mu), 5)  # Base
        fc_e <- max(max_mu / B_e, 1.0) # Fold-change
        m_e <- 4 # downtime
        l_e <- 4  # uptime
        p_e <- Mu$time[which(mu == max_mu)] # phase
        # smooth based smart starters
        DataS <- rbind(DataT[,c(1,4)], DataT[23:24,c(1,4)])
        DataS$time[25:26] <- 24
        ss <- smooth.spline(DataS$time, DataS$count, spar=0.2, cv=NA)
        ss.min <- min(ss$y)
        ss.max <- max(ss$y)
        ss.p <- which(ss$y == ss.max)
        ss.amp <- ss.max - ss.min
        ss.exp <- (ss.max - ss$y) / ss.amp
        # search for up-time
        ord <- (ss.p  - 1:13) %% 13 + 1
        ss.up <- max((ss.p - ord[which(ss.exp[c(1:13)[ord]] > 0.7)][1]) %% 12 * 2, lower['l'])
        ss.up <- min(ss.up, upper['l'])
        # search for down-time
        ord <- (1:13 - ss.p) %% 13 + 1
        ss.down <- max((ord[which(ss.exp[c(1:13)[ord]] > 0.7)][1] - ss.p) %% 12 * 2, lower['m'])
        ss.down <- min(ss.down, upper['m'])
        ss.B <- max(ss.min, 5)
        ss.fc <- max(ss.max / ss.B, 1.0)
        ss.pp <- ss$x[ss.p]
        ss.B2 <- max(mean(sort(ss$y)[1:2]), 5)
        ss.fc2 <- max(mean(sort(ss$y, decreasing=T)[1:2]) / ss.B2 , 1.0)
        # matrix of starting values (4 start points to avoid false local minima)
        starts <- matrix(nrow=5, ncol=5, dimnames=list(NULL, c('fc', 'p', 'l', 'm', 'B')))
        starts[3,] <- c(fc=fc_e, p=p_e, l=l_e, m=m_e, B=B_e)
        starts[1,] <- c(fc=ss.fc, p=ss.pp, l=ss.up, m=ss.down, B=ss.B)
        starts[2,] <- c(fc=ss.fc2, p=ss.pp, l=ss.up, m=ss.down, B=ss.B2)
        starts[4,] <- c(fc=1.0, p=ss.pp, l=8, m=8, B=mean(mu))
        starts[5,] <- c(fc=ss.fc2, p=ss.pp, l=3, m=3, B=ss.B2)
#        alt <- c(fc=fc_e/2, p=p_e, l=l_e*2, m=m_e*2, B=mean(mu))
        # weight calculation
        if (is.null(weights)) weights <- DataT$count/DataT$count
        #weights=1/rep(var,2)/sum(var)

        if (isTRUE(debug)) {
          print(DataT)
          print(Mu)
          print(starts)
          print(weights)
#         print(var3)
        #print(var2)
        }
        fits <- lapply(1:5, function(x) getFit(DataT, starts[x,], lower, upper, weights))
        sss <- as.numeric(lapply(fits, get_ss))
        if (sum(is.na(sss)) == 5) {
          fit <- NULL
	      } else {
          fit <- fits[[which(sss == min(sss, na.rm=T))[1]]]
	      }
        if (plot.it) {
          plot(count ~ time, DataT)
          points(Mu, type='b', pch=20)
          t2 <- seq(0,22,0.25)
          #points(t, fitted(fit)[1:6], col='red', type='b', pch=20)
          points(t2, sigmoid.i(coef(fit), t2), col='red', type='l', lty=3)
          if (isTRUE(debug)) {
            print(sss)
            null.fit <- lm(count ~ 1, data=DataT, weights=weights)
            print(coef(fit))
            print(c(get_aic(fit), get_aic(null.fit)))
            #print(aic_weight)
          }
        }
        return(fit)
      } else {
        return(NULL)
      }
  }
} else {
  fitSigmoid <- function (...) {
    stop("Original 'fitSigmoid' function failed due to dependency problems.")
  }
}


get_sigmoid_pars <- function (fit) {
  if (!(is.null(fit) || is.na(fit))) {
    x <- seq(0,23.995,0.005)
    pars <- coef(fit)
    series <- sigmoid.i(pars, x)
    smax <- max(series)
    smin <- min(series)
    c(pars, min=smin, max=smax, ratio=smax/smin, phase=x[which(series == smax)])
  } else {
    c(fc=NA,l=NA,m=NA,B=NA,pval=NA,min=NA,max=NA,ratio=NA,phase=NA)
  }
}


if (dependency['gtools', 'is.ok']) {
  get_sets <- function( cyclic.set, groups=2, labels = '') {
    perms <- permutations(2,groups, c(TRUE, FALSE), repeats=TRUE)
    if (labels == '') {
      labels <- LETTERS[1:groups^2]
    }
    as.factor(apply(cyclic.set, 1, function(x) {
      if (any(is.na(x))) {
        val <- NA
      } else {
        val <- labels[which(perms[,1]==x[1] & perms[,2]==x[2])]
      }
      val
    }))
  }
} else {
  get_sets <- function(...) {
    stop("Original 'get_sets' function failed due to dependency problems.")
  }
}

my.heatmap <- function (genes = c(), set='', cyclic.sets, extra_label='', scale_mode='individual', ncolors=10) {
  hmcol = brewer.pal(ncolors,"RdBu")
  pdf.x = 7.0
  pdf.y.perline = 10 / 635
  if (length(genes) == 0) {
    genes <- names(cyclic.sets)[which(cyclic.sets == set)]
  }
  if (set == '') {
    set = 'some_genes'
  }
  genes <- genes[genes %in% rownames(norm.cds.rp)]
  genes <- genes[genes %in% rownames(norm.cds.tr)]
  pdf.y = pdf.y.perline * length(genes)
  print(paste("Printing heatmap for:", set, sep=' '))
  print(paste("Length of matrix:", length(genes), sep=' '))
  testc <- cyclic.model[genes[1],]
  if (isTRUE(testc$cyclic.rp)) {
    print("Sorting by RP")
    order_by = cyclic.model.rp[genes,'phase']
  } else if (isTRUE(testc$cyclic.tr)) {
    print("Sorting by TR")
    order_by = cyclic.model.tr[genes,'phase']
  } else {
    print("Sorting by old method of AR(rp)")
    order_by = ar.rp[genes, 'phase']
  }
  ordered_genes = genes[order(-order_by)]
  if (scale_mode == 'individual') {
    rp.matrix = scale(t(as.matrix(norm.cds.rp[ordered_genes,])))
    tr.matrix = scale(t(as.matrix(norm.cds.tr[ordered_genes,])))
  }
  if (scale_mode == 'group') {
    return("Aaa")
  }
  pdf(file=paste("set",set, extra_label, scale_mode, ncolors, "heatmap.pdf",sep='_'), useDingbats=F, width=pdf.x, height=pdf.y)
  par(mar=c(0,0,0,0))
  par(mfrow=c(1,2))
  image(tr.matrix, col=rev(hmcol), xaxt='n', yaxt='n')
  image(rp.matrix, col=rev(hmcol), xaxt='n', yaxt='n')
  #   dev.copy(pdf, file=paste("set",set, scale_mode, ncol, "heatmap.pdf",sep='_'), useDingbats=F, width=pdf.x, height=pdf.y)
  dev.off()
}


my.heatmap.l <- function (genes = c(), set='', cyclic.sets, extra_label='', scale_mode='individual', ncolors=10) {
  hmcol = brewer.pal(ncolors,"RdGy")
  pdf.x = 7.0
  pdf.y.perline = 10 / 635
  if (length(genes) == 0) {
    genes <- names(cyclic.sets)[which(cyclic.sets == set)]
  }
  if (set == '') {
    set = 'some_genes'
  }
  genes <- genes[genes %in% rownames(norm.cds.l.rp)]
  genes <- genes[genes %in% rownames(norm.cds.l.tr)]
  pdf.y = pdf.y.perline * length(genes)
  print(paste("Printing heatmap for:", set, sep=' '))
  print(paste("Length of matrix:", length(genes), sep=' '))
  testc <- l.cyclic.model[genes[1],]
  if (isTRUE(testc$cyclic.l.rp)) {
    print("Sorting by RP")
    order_by = cyclic.model.l.rp[genes,'phase']
  } else if (isTRUE(testc$cyclic.l.tr)) {
    print("Sorting by TR")
    order_by = cyclic.model.l.tr[genes,'phase']
  } else {
    print("Sorting by old method of AR(rp)")
    order_by = ar.rp[genes, 'phase']
  }
  ordered_genes = genes[order(-order_by)]
  if (scale_mode == 'individual') {
    rp.matrix = scale(t(as.matrix(norm.cds.l.rp[ordered_genes,])))
    tr.matrix = scale(t(as.matrix(norm.cds.l.tr[ordered_genes,])))
  }
  if (scale_mode == 'group') {
    return("Aaa")
  }
  pdf(file=paste("set",set, extra_label, scale_mode, ncolors,length(genes), "heatmap_liver.pdf",sep='_'), useDingbats=F, width=pdf.x, height=pdf.y)
  par(mar=c(0,0,0,0))
  par(mfrow=c(1,2))
  image(tr.matrix, col=rev(hmcol), xaxt='n', yaxt='n')
  image(rp.matrix, col=rev(hmcol), xaxt='n', yaxt='n')
  #   dev.copy(pdf, file=paste("set",set, scale_mode, ncol, "heatmap.pdf",sep='_'), useDingbats=F, width=pdf.x, height=pdf.y)
  dev.off()
}

my.heatmap.k <- function (genes = c(), set='', cyclic.sets, extra_label='', scale_mode='individual', ncolors=10) {
  hmcol = brewer.pal(ncolors,"RdGy")
  pdf.x = 7.0
  pdf.y.perline = 10 / 635
  if (length(genes) == 0) {
    genes <- names(cyclic.sets)[which(cyclic.sets == set)]
  }
  if (set == '') {
    set = 'some_genes'
  }
  genes <- genes[genes %in% rownames(norm.cds.k.rp)]
  genes <- genes[genes %in% rownames(norm.cds.k.tr)]
  pdf.y = pdf.y.perline * length(genes)
  print(paste("Printing heatmap for:", set, sep=' '))
  print(paste("Length of matrix:", length(genes), sep=' '))
  testc <- k.cyclic.model[genes[1],]
  if (isTRUE(testc$cyclic.k.rp)) {
    print("Sorting by RP")
    order_by = cyclic.model.k.rp[genes,'phase']
  } else if (isTRUE(testc$cyclic.k.tr)) {
    print("Sorting by TR")
    order_by = cyclic.model.k.tr[genes,'phase']
  } else {
    print("Sorting by old method of AR(rp)")
    order_by = ar.rp[genes, 'phase']
  }
  ordered_genes = genes[order(-order_by)]
  if (scale_mode == 'individual') {
    rp.matrix = scale(t(as.matrix(norm.cds.k.rp[ordered_genes,])))
    tr.matrix = scale(t(as.matrix(norm.cds.k.tr[ordered_genes,])))
  }
  if (scale_mode == 'group') {
    return("Aaa")
  }
  pdf(file=paste("set",set, extra_label, scale_mode, ncolors,length(genes), "heatmap_kidney.pdf",sep='_'), useDingbats=F, width=pdf.x, height=pdf.y)
  par(mar=c(0,0,0,0))
  par(mfrow=c(1,2))
  image(tr.matrix, col=rev(hmcol), xaxt='n', yaxt='n')
  image(rp.matrix, col=rev(hmcol), xaxt='n', yaxt='n')
  #   dev.copy(pdf, file=paste("set",set, scale_mode, ncol, "heatmap.pdf",sep='_'), useDingbats=F, width=pdf.x, height=pdf.y)
  dev.off()
}
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # Make the panel
  plotCols = cols   # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # #rows calculated from #cols
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    pushViewport(vplayout(curRow, curCol))
    grid.draw(plots[[i]])
    popViewport()
  }
}
proplot.liver <- function (identifier='', setfolder=FALSE, subfolder='') {
  get_fig <- function(plotdata, fitdata, max_range, ytitle, ptitle, colors, fitline='1343', sigdata=NULL, fitconf=NULL) {
    gfig <- ggplot(data=plotdata, aes(time, count, ymin=count1, ymax=count2, group=treatment, colour=treatment))
    if (! is.null(fitconf)) gfig <- gfig + geom_polygon(data=fitconf, aes(x=x, y=y), fill='#eeeeee')
    gfig <- gfig + geom_line(data=fitdata, aes(x=x,y=y), linetype=fitline, size=1.2)
    if (! is.null(sigdata)) gfig <- gfig + geom_line(data=sigdata, aes(x=x,y=y), linetype=2, size=1.2)
    gfig <- gfig +
      geom_line(size=1.5) +
      geom_pointrange(aes(shape=treatment), size=1.5) +
      scale_y_continuous(ytitle,
                         limits=max_range,
                         expand=c(0,0) ) +
      scale_x_continuous("Zeitgeber time",
                         breaks = res, expand = c(0,0)) +
      scale_color_manual(values=colors) +
      ggtitle(ptitle) +
      theme_bw() +
      theme(plot.title = element_text(size=30, lineheight=.8, face="bold")) +
      theme(axis.text.x = element_text(colour='black', size=16),
            axis.text.y = element_text(colour='black', size=16)) +
      theme(legend.position='none')
    gfig
  }
  get_empty_fig <- function(ptitle) {
    ggplot() +
      geom_text(aes(x=1, y=1, label='Not available')) +
      ggtitle(ptitle) +
      theme(plot.title = element_text(size=30, lineheight=.8, face="bold")) +
      theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  }
  get_double_axis <- function(plot1, plot2) {
    # build gtables from the two plots
    g1 <- ggplot_gtable(ggplot_build(plot1))
    g2 <- ggplot_gtable(ggplot_build(plot2))
    # extract y-axis from second and modify it
    new_y_index <- which(g2$layout$name == "axis-l")
    new_y <- g2$grobs[[new_y_index]]$children[[2]]
    new_y$widths <- rev(new_y$widths)
    new_y$grobs <- rev(new_y$grobs)
    if (is.null(new_y$grobs[[1]]$x)) new_y$grobs[[1]]$x <- rep(unit(1, "npc") - unit(0.15, "cm"), 4)
    new_y$grobs[[1]]$x <- new_y$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
    # create a new grid object from the first plot and add modified panel with new y-axis
    panel <- c(subset(g1$layout, name == "panel", t:r))
    new_g <- gtable_add_grob(g1, g1$grobs[[which(g1$layout$name == "panel")]], panel$t,
                             panel$l, panel$b, panel$l)
    new_g <- gtable_add_cols(new_g, g2$widths[g2$layout[new_y_index, ]$l], length(new_g$widths) - 1)
    new_g <- gtable_add_grob(new_g, new_y, panel$t, length(new_g$widths) - 1, panel$b)
    set_fixed_size(g=new_g)
  }
  set_fixed_size <- function(p=NULL, g=ggplotGrob(p), lwidth=unit(0.55, "in"), pwidth=unit(5, "in"), height=unit(4, "in")){
    panel_index_w<- g$layout$l[g$layout$name=="panel"]
    panel_index_h<- g$layout$t[g$layout$name=="panel"]
    g$widths[[panel_index_w]] <- pwidth
    g$heights[[panel_index_h]] <- height
    axis_l_index <- g$layout$r[g$layout$name=="axis-l"]
    #ylab_index <- g$layout$r[g$layout$name=="ylab"]
    g$widths[[axis_l_index]] <- lwidth
    #g$widths[[ylab_index]] <- lwidth
    class(g) <- c("fixed", class(g), "ggplot")
    g
  }

   if (identifier != '') {
     gene_name = annot_gene_names[identifier,1]
     if (is.na(gene_name)) {
       nidentifier <- row.names(annot_gene_names)[which(annot_gene_names == identifier)]
       if (!length(nidentifier)) stop(paste('Could not find a gene with the supplied identifier:', identifier))
       identifier <- nidentifier
       gene_name = annot_gene_names[identifier,1]
     }
   } else {
     stop('Either a gene ID or name should be supplied as identifier')
   }

  # Variables
  print(paste('Processing',identifier, gene_name, sep=':'))

  path = '.'
  if (subfolder != '') path = paste(path,subfolder,sep='/')
  if (setfolder) path=paste(path,l.cyclic.model.sets[identifier],sep='/')
  path = paste(path,'',sep='/')
  res <- seq(0,22,2)
  t_points <- res
  x_points <- seq(-2,24,0.5)
  reps <- as.data.frame(strsplit(colnames(norm.cds.l.rp), '[_.]'))[2,]
  rep1 <- which(reps == 1)
  rep2 <- which(reps == 2)
  empty_fit <- data.frame(x=NaN, y=NaN, count1=0, count2=0, treatment=c('rpf','trf'))
  empty_fit.sig <- data.frame(x=NaN, y=NaN, count1=0, count2=0, treatment=c('rpfs','trfs'))

  # RP data processing
  no_rp <- identifier %ni% rownames(rpkm.cds.l.rp)
  if (!no_rp) {
    rp <- as.numeric(rpkm.cds.l.rp[identifier,])
    rp.1 <- rp[rep1]
    rp.2 <- rp[rep2]
    rp_data <- data.frame(time=t_points, count=(rp.1 + rp.2)/2,
                            count1=rp.1, count2=rp.2, treatment='rp')
    rp_mean <- max(mean(rp_data$count), 1)
    rp_max_val <- max(rp.1, rp.2)
    rp_max_fold <- rp_max_val / rp_mean

    what_to_plot <- as.character(cyclic.model.l.rp[identifier,'decision'])
    if (what_to_plot == 'sinusoidal') {
      rp_fit <- data.frame(x=x_points, y=(model.l.rp[identifier,]$sinusoidal.mean + cos(pi/12*(x_points - model.l.rp[identifier,]$sinusoidal.phase))*model.l.rp[identifier,]$sinusoidal.amp)/l.composite.sizes[identifier, 'cds']/libsize.l.rp*1e+09,
                           count1=0, count2=0, treatment=rep('rpf',length(x_points)))
    } else if (what_to_plot == 'sigmoid') {
      rp.sig.pars <- model.l.rp[identifier,1:5]
      names(rp.sig.pars) <- c('fc', 'p', 'l', 'm', 'B')
      rp_fit <- data.frame(x=x_points, y=sigmoid.i(rp.sig.pars, x_points)/l.composite.sizes[identifier, 'cds']/libsize.l.rp*1e+09, count1=0, count2=0,
                                   treatment=rep('rpfs',length(x_points)))
    } else {
      rp_fit <- data.frame(x=x_points, y=model.l.rp[identifier,]$null.intercept/l.composite.sizes[identifier, 'cds']/libsize.l.rp*1e+09, count1=0, count2=0, treatment=rep('rpf',length(x_points)))
    }
  } else {
    rp_max_fold = 1
    rp_mean = 1
  }

  # TR data processing
  no_tr <- identifier %ni% rownames(rpkm.cds.l.tr)
  if (!no_tr) {
    tr <- as.numeric(rpkm.cds.l.tr[identifier,])
    tr.1 <- tr[rep1]
    tr.2 <- tr[rep2]
    tr_data <- data.frame(time=t_points, count=(tr.1 + tr.2)/2,
                          count1=tr.1, count2=tr.2, treatment='tr')
    tr_mean <- max(mean(tr_data$count), 1)
    tr_max_val <- max(tr.1, tr.2)
    tr_max_fold <- tr_max_val / tr_mean

    what_to_plot <- as.character(cyclic.model.l.tr[identifier,'decision'])
    if (what_to_plot =='sinusoidal') {
      tr_fit <- data.frame(x=x_points, y=(model.l.tr[identifier,]$sinusoidal.mean +
                                            cos(pi/12*(x_points - model.l.tr[identifier,]$sinusoidal.phase))*model.l.tr[identifier,]$sinusoidal.amp) /
                                            l.composite.sizes[identifier, 'cds']/libsize.l.tr*1e+09,
                           count1=0, count2=0, treatment=rep('trf',length(x_points)))
    } else if (what_to_plot == 'sigmoid') {
      tr.sig.pars <- model.l.tr[identifier,1:5]
      names(tr.sig.pars) <- c('fc', 'p', 'l', 'm', 'B')
      tr_fit <- data.frame(x=x_points, y=sigmoid.i(tr.sig.pars, x_points)/l.composite.sizes[identifier, 'cds']/libsize.l.tr*1e+09, count1=0, count2=0,
                                   treatment=rep('trf',length(x_points)))
    } else {
      tr_fit <- data.frame(x=x_points, y=model.l.tr[identifier,]$null.intercept/l.composite.sizes[identifier, 'cds']/libsize.l.tr*1e+09, count1=0, count2=0, treatment=rep('trf',length(x_points)))
    }
  } else {
    tr_max_fold = 1
    tr_mean = 1
  }

  # Setting the scales of the plot based on mRNA and premRNA range & means
  eff_fold <- ceiling(max(rp_max_fold, tr_max_fold))
  max_rp_range <- eff_fold * rp_mean
  max_tr_range <- eff_fold * tr_mean

  # Create the RP figure
  if (!no_rp) {
    if (!no_tr) {
      rp.focal_data <- rbind(rp_data, cbind(tr_data[1], tr_data[2:4]*max_rp_range/max_tr_range, tr_data[5]))
      rp.focal_fit <- rbind(rp_fit, cbind(tr_fit[1], tr_fit[2]*max_rp_range/max_tr_range, tr_fit[3:5]))
    } else {
      rp.focal_data <- rbind(rp_data, rp_data)
      rp.focal_fit <- rbind(rp_fit, rp_fit)
    }
    rp_fig.counts <- get_fig(rp.focal_data,rp.focal_fit, c(0,max_rp_range), 'RPKM', gene_name, colors=c("#0066cc","#5aadff","#cc6600","#ffad5a")) #, rp.focal_sigmoid)
  } else {
    rp_fig.counts <- get_empty_fig(paste(gene_name,"RP"))
  }
  if (!no_tr) {
    if (!no_rp) {
      tr.focal_data <- rbind(tr_data, cbind(rp_data[1], rp_data[2:4]*max_tr_range/max_rp_range, rp_data[5]))
      tr.focal_fit <- rbind(tr_fit, cbind(rp_fit[1], rp_fit[2]*max_tr_range/max_rp_range, rp_fit[3:5]))
    } else  {
      tr.focal_data <- rbind(tr_data, tr_data)
      tr.focal_fit <- rbind(tr_fit, tr_fit)
    }
    tr_fig.counts <- get_fig(tr.focal_data, tr.focal_fit, c(0,max_tr_range), 'RPKM', gene_name, colors=c("#0066cc","#5aadff","#cc6600","#ffad5a"), fitline='4313') #, tr.focal_sigmoid)

  } else {
    tr_fig.counts <- get_empty_fig(paste(gene_name,"TR"))
  }


  # Secondary Y-axis crazy hacks!
  profile_fig <- get_double_axis(tr_fig.counts, rp_fig.counts)

  # efficiency plot
  no_eff <- identifier %ni% rownames(l.eff.noinf)
  if (!no_eff) {
    eff <- as.numeric(l.eff.log.noinf[identifier,])
    eff.1 <- eff[rep1]
    eff.2 <- eff[rep2]
    eff_data <- data.frame(time=t_points, count=l.eff.log.noinf.mean[identifier,],  #(eff.1 + eff.2)/2,
                          count1=eff.1, count2=eff.2, treatment='none')
    eff_mean <- 2^mean(eff_data$count, na.rm=T)
    eff_max_val <- max(eff.1, eff.2, na.rm=T)
    eff_min_val <- min(eff.1, eff.2, na.rm=T)
    eff_max_fold <- 2^(eff_max_val - eff_min_val)
    if (l.cyclic.eff.fdr[identifier]<0.05) {
      eff_fit <- l.cyclic.eff.fits[[identifier]]
    } else {
      eff_fit <- update(l.cyclic.eff.fits[[identifier]], .~1)
    }
    if (l.cyclic.eff.decision[identifier]) {
      fitcol <- c('#000000','#444444','#eeeeee')
    } else {
      fitcol <- c('#000000','#bbbbbb','#eeeeee')
    }
    eff_predict <- as.data.frame(predict(eff_fit, data.frame(cosx=cos(pi/12 * x_points), sinx=sin(pi/12 * x_points)), interval = 'confidence'))
    eff.cyc <- l.cyclic.eff.pars[identifier,]
    if (is.na(eff.cyc$amp)) {
      eff.cyc$amp <- 0
      eff.cyc$phase <- 1
      eff.cyc$period <- 1
    }
    neff_fold <- ceiling(max(rp_max_fold, tr_max_fold, eff_max_fold))
    max_eff_range <- log2(neff_fold * eff_mean)
    min_eff_range <- log2(eff_mean / neff_fold)
    eff_fit <- data.frame(x=x_points, y=eff_predict$fit,
                         count1=0, count2=0, treatment=rep('nonef',length(x_points)))
    conf_data <- data.frame(x=c(x_points, rev(x_points)), y=c(eff_predict$upr, rev(eff_predict$lwr)), count1=0, count2=0, treatment=rep('nonez', 2*length(x_points)))
    eff_fig.counts <- set_fixed_size(p=get_fig(eff_data,eff_fit, c(min_eff_range,max_eff_range), 'log2(Ribosome occupancy)', gene_name, colors=fitcol, fitconf=conf_data))
  } else {
    eff_fig.counts <- set_fixed_size(p=get_empty_fig(paste(gene_name,"log2(Ribosome occupancy)")))
  }
  what_to_plot <- c('set', as.character(l.cyclic.model.sets[identifier]))
  pdf(paste(path, gene_name, '_liver_', paste(what_to_plot,collapse='_'), '.pdf', sep=''), width=14, height=7, useDingbats=FALSE)
  suppressMessages(multiplot(profile_fig, eff_fig.counts, cols=2))
  garbage <- dev.off()
}

proplot.kidney <- function (identifier='', setfolder=FALSE, subfolder='') {
  get_fig <- function(plotdata, fitdata, max_range, ytitle, ptitle, colors, fitline='1343', sigdata=NULL, fitconf=NULL) {
    gfig <- ggplot(data=plotdata, aes(time, count, ymin=count1, ymax=count2, group=treatment, colour=treatment))
    if (! is.null(fitconf)) gfig <- gfig + geom_polygon(data=fitconf, aes(x=x, y=y), fill='#eeeeee')
    gfig <- gfig + geom_line(data=fitdata, aes(x=x,y=y), linetype=fitline, size=1.2)
    if (! is.null(sigdata)) gfig <- gfig + geom_line(data=sigdata, aes(x=x,y=y), linetype=2, size=1.2)
    gfig <- gfig +
      geom_line(size=1.5) +
      geom_pointrange(aes(shape=treatment), size=1.5) +
      scale_y_continuous(ytitle,
                         limits=max_range,
                         expand=c(0,0) ) +
      scale_x_continuous("Zeitgeber time",
                         breaks = res, expand = c(0,0)) +
      scale_color_manual(values=colors) +
      ggtitle(ptitle) +
      theme_bw() +
      theme(plot.title = element_text(size=30, lineheight=.8, face="bold")) +
      theme(axis.text.x = element_text(colour='black', size=16),
            axis.text.y = element_text(colour='black', size=16)) +
      theme(legend.position='none')
    gfig
  }
  get_empty_fig <- function(ptitle) {
    ggplot() +
      geom_text(aes(x=1, y=1, label='Not available')) +
      ggtitle(ptitle) +
      theme(plot.title = element_text(size=30, lineheight=.8, face="bold")) +
      theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  }
  get_double_axis <- function(plot1, plot2) {
    # build gtables from the two plots
    g1 <- ggplot_gtable(ggplot_build(plot1))
    g2 <- ggplot_gtable(ggplot_build(plot2))
    # extract y-axis from second and modify it
    new_y_index <- which(g2$layout$name == "axis-l")
    new_y <- g2$grobs[[new_y_index]]$children[[2]]
    new_y$widths <- rev(new_y$widths)
    new_y$grobs <- rev(new_y$grobs)
    if (is.null(new_y$grobs[[1]]$x)) new_y$grobs[[1]]$x <- rep(unit(1, "npc") - unit(0.15, "cm"), 4)
    new_y$grobs[[1]]$x <- new_y$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
    # create a new grid object from the first plot and add modified panel with new y-axis
    panel <- c(subset(g1$layout, name == "panel", t:r))
    new_g <- gtable_add_grob(g1, g1$grobs[[which(g1$layout$name == "panel")]], panel$t,
                             panel$l, panel$b, panel$l)
    new_g <- gtable_add_cols(new_g, g2$widths[g2$layout[new_y_index, ]$l], length(new_g$widths) - 1)
    new_g <- gtable_add_grob(new_g, new_y, panel$t, length(new_g$widths) - 1, panel$b)
    set_fixed_size(g=new_g)
  }
  set_fixed_size <- function(p=NULL, g=ggplotGrob(p), lwidth=unit(0.55, "in"), pwidth=unit(5, "in"), height=unit(4, "in")){
    panel_index_w<- g$layout$l[g$layout$name=="panel"]
    panel_index_h<- g$layout$t[g$layout$name=="panel"]
    g$widths[[panel_index_w]] <- pwidth
    g$heights[[panel_index_h]] <- height
    axis_l_index <- g$layout$r[g$layout$name=="axis-l"]
    #ylab_index <- g$layout$r[g$layout$name=="ylab"]
    g$widths[[axis_l_index]] <- lwidth
    #g$widths[[ylab_index]] <- lwidth
    class(g) <- c("fixed", class(g), "ggplot")
    g
  }
  
  if (identifier != '') {
    gene_name = annot_gene_names[identifier,1]
    if (is.na(gene_name)) {
      nidentifier <- row.names(annot_gene_names)[which(annot_gene_names == identifier)]
      if (!length(nidentifier)) stop(paste('Could not find a gene with the supplied identifier:', identifier))
      identifier <- nidentifier
      gene_name = annot_gene_names[identifier,1]
    }
  } else {
    stop('Either a gene ID or name should be supplied as identifier')
  }
  
  # Variables
  print(paste('Processing',identifier, gene_name, sep=':'))
  
  path = '.'
  if (subfolder != '') path = paste(path,subfolder,sep='/')
  if (setfolder) path=paste(path,k.cyclic.model.sets[identifier],sep='/')
  path = paste(path,'',sep='/')
  res <- seq(0,22,2)
  t_points <- res
  x_points <- seq(-2,24,0.5)
  reps <- as.data.frame(strsplit(colnames(norm.cds.k.rp), '[_.]'))[2,]
  rep1 <- which(reps == 1)
  rep2 <- which(reps == 2)
  empty_fit <- data.frame(x=NaN, y=NaN, count1=0, count2=0, treatment=c('rpf','trf'))
  empty_fit.sig <- data.frame(x=NaN, y=NaN, count1=0, count2=0, treatment=c('rpfs','trfs'))
  
  # RP data processing
  no_rp <- identifier %ni% rownames(rpkm.cds.k.rp)
  if (!no_rp) {
    rp <- as.numeric(rpkm.cds.k.rp[identifier,])
    rp.1 <- rp[rep1]
    rp.2 <- rp[rep2]
    rp_data <- data.frame(time=t_points, count=(rp.1 + rp.2)/2,
                          count1=rp.1, count2=rp.2, treatment='rp')
    rp_mean <- max(mean(rp_data$count), 1)
    rp_max_val <- max(rp.1, rp.2)
    rp_max_fold <- rp_max_val / rp_mean
    
    what_to_plot <- as.character(cyclic.model.k.rp[identifier,'decision'])
    if (what_to_plot == 'sinusoidal') {
      rp_fit <- data.frame(x=x_points, y=(model.k.rp[identifier,]$sinusoidal.mean + cos(pi/12*(x_points - model.k.rp[identifier,]$sinusoidal.phase))*model.k.rp[identifier,]$sinusoidal.amp)/k.composite.sizes[identifier, 'cds']/libsize.k.rp*1e+09,
                           count1=0, count2=0, treatment=rep('rpf',length(x_points)))
    } else if (what_to_plot == 'sigmoid') {
      rp.sig.pars <- model.k.rp[identifier,1:5]
      names(rp.sig.pars) <- c('fc', 'p', 'l', 'm', 'B')
      rp_fit <- data.frame(x=x_points, y=sigmoid.i(rp.sig.pars, x_points)/k.composite.sizes[identifier, 'cds']/libsize.k.rp*1e+09, count1=0, count2=0,
                           treatment=rep('rpfs',length(x_points)))
    } else {
      rp_fit <- data.frame(x=x_points, y=model.k.rp[identifier,]$null.intercept/k.composite.sizes[identifier, 'cds']/libsize.k.rp*1e+09, count1=0, count2=0, treatment=rep('rpf',length(x_points)))
    }
  } else {
    rp_max_fold = 1
    rp_mean = 1
  }
  
  # TR data processing
  no_tr <- identifier %ni% rownames(rpkm.cds.k.tr)
  if (!no_tr) {
    tr <- as.numeric(rpkm.cds.k.tr[identifier,])
    tr.1 <- tr[rep1]
    tr.2 <- tr[rep2]
    tr_data <- data.frame(time=t_points, count=(tr.1 + tr.2)/2,
                          count1=tr.1, count2=tr.2, treatment='tr')
    tr_mean <- max(mean(tr_data$count), 1)
    tr_max_val <- max(tr.1, tr.2)
    tr_max_fold <- tr_max_val / tr_mean
    
    what_to_plot <- as.character(cyclic.model.k.tr[identifier,'decision'])
    if (what_to_plot =='sinusoidal') {
      tr_fit <- data.frame(x=x_points, y=(model.k.tr[identifier,]$sinusoidal.mean +
                                            cos(pi/12*(x_points - model.k.tr[identifier,]$sinusoidal.phase))*model.k.tr[identifier,]$sinusoidal.amp) /
                             k.composite.sizes[identifier, 'cds']/libsize.k.tr*1e+09,
                           count1=0, count2=0, treatment=rep('trf',length(x_points)))
    } else if (what_to_plot == 'sigmoid') {
      tr.sig.pars <- model.k.tr[identifier,1:5]
      names(tr.sig.pars) <- c('fc', 'p', 'l', 'm', 'B')
      tr_fit <- data.frame(x=x_points, y=sigmoid.i(tr.sig.pars, x_points)/k.composite.sizes[identifier, 'cds']/libsize.k.tr*1e+09, count1=0, count2=0,
                           treatment=rep('trf',length(x_points)))
    } else {
      tr_fit <- data.frame(x=x_points, y=model.k.tr[identifier,]$null.intercept/k.composite.sizes[identifier, 'cds']/libsize.k.tr*1e+09, count1=0, count2=0, treatment=rep('trf',length(x_points)))
    }
  } else {
    tr_max_fold = 1
    tr_mean = 1
  }
  
  # Setting the scales of the plot based on mRNA and premRNA range & means
  eff_fold <- ceiling(max(rp_max_fold, tr_max_fold))
  max_rp_range <- eff_fold * rp_mean
  max_tr_range <- eff_fold * tr_mean
  
  # Create the RP figure
  if (!no_rp) {
    if (!no_tr) {
      rp.focal_data <- rbind(rp_data, cbind(tr_data[1], tr_data[2:4]*max_rp_range/max_tr_range, tr_data[5]))
      rp.focal_fit <- rbind(rp_fit, cbind(tr_fit[1], tr_fit[2]*max_rp_range/max_tr_range, tr_fit[3:5]))
    } else {
      rp.focal_data <- rbind(rp_data, rp_data)
      rp.focal_fit <- rbind(rp_fit, rp_fit)
    }
    rp_fig.counts <- get_fig(rp.focal_data,rp.focal_fit, c(0,max_rp_range), 'RPKM', gene_name, colors=c("#0066cc","#5aadff","#cc6600","#ffad5a")) #, rp.focal_sigmoid)
  } else {
    rp_fig.counts <- get_empty_fig(paste(gene_name,"RP"))
  }
  if (!no_tr) {
    if (!no_rp) {
      tr.focal_data <- rbind(tr_data, cbind(rp_data[1], rp_data[2:4]*max_tr_range/max_rp_range, rp_data[5]))
      tr.focal_fit <- rbind(tr_fit, cbind(rp_fit[1], rp_fit[2]*max_tr_range/max_rp_range, rp_fit[3:5]))
    } else  {
      tr.focal_data <- rbind(tr_data, tr_data)
      tr.focal_fit <- rbind(tr_fit, tr_fit)
    }
    tr_fig.counts <- get_fig(tr.focal_data, tr.focal_fit, c(0,max_tr_range), 'RPKM', gene_name, colors=c("#0066cc","#5aadff","#cc6600","#ffad5a"), fitline='4313') #, tr.focal_sigmoid)
    
  } else {
    tr_fig.counts <- get_empty_fig(paste(gene_name,"TR"))
  }
  
  
  # Secondary Y-axis crazy hacks!
  profile_fig <- get_double_axis(tr_fig.counts, rp_fig.counts)
  
  # efficiency plot
  no_eff <- identifier %ni% rownames(k.eff.noinf)
  if (!no_eff) {
    eff <- as.numeric(k.eff.log.noinf[identifier,])
    eff.1 <- eff[rep1]
    eff.2 <- eff[rep2]
    eff_data <- data.frame(time=t_points, count=k.eff.log.noinf.mean[identifier,],  #(eff.1 + eff.2)/2,
                           count1=eff.1, count2=eff.2, treatment='none')
    eff_mean <- 2^mean(eff_data$count, na.rm=T)
    eff_max_val <- max(eff.1, eff.2, na.rm=T)
    eff_min_val <- min(eff.1, eff.2, na.rm=T)
    eff_max_fold <- 2^(eff_max_val - eff_min_val)
    if (k.cyclic.eff.fdr[identifier]<0.05) {
      eff_fit <- k.cyclic.eff.fits[[identifier]]
    } else {
      eff_fit <- update(k.cyclic.eff.fits[[identifier]], .~1)
    }
    if (k.cyclic.eff.decision[identifier]) {
      fitcol <- c('#000000','#444444','#eeeeee')
    } else {
      fitcol <- c('#000000','#bbbbbb','#eeeeee')
    }
    eff_predict <- as.data.frame(predict(eff_fit, data.frame(cosx=cos(pi/12 * x_points), sinx=sin(pi/12 * x_points)), interval = 'confidence'))
    eff.cyc <- k.cyclic.eff.pars[identifier,]
    if (is.na(eff.cyc$amp)) {
      eff.cyc$amp <- 0
      eff.cyc$phase <- 1
      eff.cyc$period <- 1
    }
    neff_fold <- ceiling(max(rp_max_fold, tr_max_fold, eff_max_fold))
    max_eff_range <- log2(neff_fold * eff_mean)
    min_eff_range <- log2(eff_mean / neff_fold)
    eff_fit <- data.frame(x=x_points, y=eff_predict$fit,
                          count1=0, count2=0, treatment=rep('nonef',length(x_points)))
    conf_data <- data.frame(x=c(x_points, rev(x_points)), y=c(eff_predict$upr, rev(eff_predict$lwr)), count1=0, count2=0, treatment=rep('nonez', 2*length(x_points)))
    eff_fig.counts <- set_fixed_size(p=get_fig(eff_data,eff_fit, c(min_eff_range,max_eff_range), 'log2(Ribosome occupancy)', gene_name, colors=fitcol, fitconf=conf_data))
  } else {
    eff_fig.counts <- set_fixed_size(p=get_empty_fig(paste(gene_name,"log2(Ribosome occupancy)")))
  }
  what_to_plot <- c('set', as.character(k.cyclic.model.sets[identifier]))
  pdf(paste(path, gene_name, '_kidney_', paste(what_to_plot,collapse='_'), '.pdf', sep=''), width=14, height=7, useDingbats=FALSE)
  suppressMessages(multiplot(profile_fig, eff_fig.counts, cols=2))
  garbage <- dev.off()
}

make_circ <- function (x) {
  if (is.na(x)) {
    return(NA)
  }
  if (x > -12 & x <= 12) {
    new.x <- x
  }
  if (x > 12) {
    new.x <- x - 24
  }
  if (x <= -12) {
    new.x <- x + 24
  }
  new.x
}

phase_biplot <- function(genelist, dataset='ar', xorder='tr',xproj='rp', connect=FALSE, circ=T, legend=T, pch=20, cex=0.6, title=NULL,color=F,color1='red',color2='blue',highlight=NULL) {
  tr <- get(paste(dataset,xorder,sep='.'))[genelist,]
  tr <- tr[order(tr$phase),]
  rp <- get(paste(dataset,xproj,sep='.'))[rownames(tr),] #rp is ordered like tr (so by phase)
  selected = genelist[match(rownames(tr),genelist)] #order genelist as tr
  if (isTRUE(circ)) {
    x1 = as.numeric(lapply(tr$phase, make_circ))
    x2 = as.numeric(lapply(rp$phase, make_circ))
    xaxp = c(-12,12,4)
  } else {
    x1 = tr$phase
    x2 = rp$phase
    xaxp=c(0,24,4)
  }
  print(length(x1))
  print(length(x2))
  y = 1:length(genelist)
   plot(x1, y, pch=pch, cex=cex, xlab='Phase (h)', ylab=paste('Genes, N=', length(genelist), sep=''),main=title, xaxp=xaxp, col=color1)
  if (connect) { #not used
    segments(x0=x1, y0=y, x1=x2, y1=y, col='red', lwd=0.5)
    points(x1, y, pch=pch, cex=cex, col='red') #not used (#vio)
  }
  if (!isTRUE(color)) { #if color=FALSE
    points(x2, y, pch=pch, cex=cex, col=color2)
  } else { #if color=TRUE
    nb=1
    for (i in y) {
      if (selected[i] %in% highlight[,2]){
        points(x2[i], i, pch=pch, cex=0.8, col='blue')
        text(x2[i],i,highlight[highlight$gid==selected[i],1],pos=1,cex=1)
        nb=nb+1
      } else {
        points(x2[i], i, pch=pch, cex=cex, col=color2)
      }
      }
   
#     ncols <- 6
#     mycol <- rev(brewer.pal(ncols,"RdBu"))
#     #mycol <- c('blue','green','orange','red','purple','cyan')
#     coldata <- log2(rp$ratio/tr$ratio)
#     colbox <- boxplot(coldata, plot=F)
#     cols <- rep(ncols, length(coldata))
#     for (i in (ncols-1):1) {
#       cols[which(coldata < colbox$stats[i])] = i
#     }
#     print(mycol)
#     #print(coldata)
#     #print(colbox$stats)
#     #print(cols)
#     points(x2, y, pch=pch, cex=cex, col=mycol[cols])
  }
  if (isTRUE(legend)) legend('bottom',legend=c('RNAseq','RPFseq'), col=c(color1,color2), pch=c(20,20),cex=0.8,ncol=2)
}


log_scatter <- function (x, y, steps=c(1,1), ...) {
  plot(x, y, log='xy', xaxt='n', yaxt='n', ...)
  usr <- par('usr')
  xr <- (usr[2] - usr[1]) / 27 # 27 = (100 + 2*4) / 4
  yr <- (usr[4] - usr[3]) / 27
  xlim <- c(usr[1] + xr, usr[2] - xr)
  ylim <- c(usr[3] + yr, usr[4] - yr)
  aty <- log10(axTicks(2))
  atx <- log10(axTicks(1))
  aty <- seq(min(aty), floor(ylim[2]), steps[2])
  atx <- seq(min(atx), floor(xlim[2]), steps[1])
  labelsx <- sapply(atx, function(x) as.expression(bquote(10^ .(x))))
  labelsy <- sapply(aty, function(x) as.expression(bquote(10^ .(x))))
  axis(2,at=10^aty,labels=labelsy)
  axis(1,at=10^atx,labels=labelsx)
}

log_scatter2 <- function (x, y, steps=c(1,1), log.ax='xy', ...) {
  needx <- grepl('x', log.ax)
  needy <- grepl('y', log.ax)
  if (needx) xaxt = 'n' else xaxt='s'
  if (needy) yaxt = 'n' else yaxt='s'
  plot(x, y, log=log.ax, xaxt=xaxt, yaxt=yaxt, ...)
  usr <- par('usr')
  if (needx) {
    xr <- (usr[2] - usr[1]) / 27 # 27 = (100 + 2*4) / 4
    xlim <- c(usr[1] + xr, usr[2] - xr)
    atx <- log10(axTicks(1))
    atx <- seq(min(atx), floor(xlim[2]), steps[1])
    labelsx <- sapply(atx, function(x) as.expression(bquote(10^ .(x))))
    axis(1,at=10^atx,labels=labelsx)
  }
  if (needy) {
    yr <- (usr[4] - usr[3]) / 27
    ylim <- c(usr[3] + yr, usr[4] - yr)
    aty <- log10(axTicks(2))
    aty <- seq(min(aty), floor(ylim[2]), steps[2])
    labelsy <- sapply(aty, function(x) as.expression(bquote(10^ .(x))))
    axis(2,at=10^aty,labels=labelsy)
  }
}

lmp <- function (modelobject, label='') {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  s <- summary(modelobject)
  f <- s$fstatistic
  f.p <- pf(f[1],f[2],f[3],lower.tail=F)
  t <- summary(modelobject)$coefficients
  t.p <- t[2,4]
  t.t <- t[2,3]
  b <- t[2,1]
  leg.text <- bquote(paste( .(label), ': ', italic(b),'=', .(round(b,3)), ', ', italic(t), '(', .(f[3]), ')=', .(round(abs(t.t), 2)),
                          ', ', italic(p), '=', .(signif(t.p,3)), '; ',
                          italic(R)^2,'=', .(round(s$r.squared,3)), ', ', italic(F), '(', .(f[2]), ',',
                          .(f[3]), ')=', .(round(f[1],2)), ', ', italic(p), '=', .(signif(f.p,3)) ))
  return(leg.text)
}

lm.fpval <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  s <- summary(modelobject)
  f <- s$fstatistic
  if (length(f)) return(pf(f[1],f[2],f[3],lower.tail=F))
  return(NA)
}

scatter_with_marginal_densities <- function(d.f, ordered_keys, my.cols, margin.xy='xy', multiplier = 0.8, multiplier.x=multiplier,
                                            multiplier.y=multiplier, fg.cols=NA, den.cols=NA, lwd.den=NA, main='', xlab='', ylab='', 
                                            label.it=F, xlim=c(), ylim=c(), log.xy='xy', den.solid=T, pch=21, ...) {
  dots <- list(...)
  nmdots <- names(dots)
  for (ax in c('x','y')) {
    if (grepl(ax, margin.xy)) {
      density.m <- list()
      i = 1
      for (key in ordered_keys) {
        if (key == 'All' | key == 'all') {
          print('Special key: producing *density* for *all* instead of key matching')
          density.m[[i]] <- density(get(ax, d.f), na.rm=TRUE)
        } else {
          density.m[[i]] <- density(subset(d.f, k==key, get(ax))[,1], na.rm=TRUE)
        }
        i = i + 1
      }
      assign(paste('density',ax, sep='.'), density.m)
    }
  }
  if (length(xlim) == 0) {
    xlim = c(min(d.f$x), max(d.f$x))
  }
  if (length(ylim) == 0) {
    ylim = c(min(d.f$y), max(d.f$y))
  }
  # canvas
  if (grepl('x|y', log.xy)) {
    log_scatter2(1,1, log.ax=log.xy, ylim=ylim, xlim=xlim, type='n', main=main, xlab=xlab, ylab=ylab)
    grid('grey')
  } else {
    plot(1, 1, ylim=ylim, xlim=xlim, type='n', main=main, xlab=xlab, ylab=ylab) #eqscplot for square plots
    pretty_grid(2, col='grey')
    #diagonal_grid(interval=2, col='lightgrey', lty='dotdash')
  }
  real_limits <- get_par_limits(p=2)
  xlim <- real_limits$xlim
  ylim <- real_limits$ylim
  # scatterplots
  i = 1
  if (pch < 21 | pch > 25) fg.cols = fg.cols
  for (key in ordered_keys) {
    points(subset(d.f, k==key, 1:2), , type='p', col=fg.cols[i], bg=my.cols[i], pch=pch, ...)
    i = i + 1
  }
  # marginal densities
  if (length(den.cols) == 1 && is.na(den.cols)) den.cols <- my.cols
  if (den.solid) den.cols = paste(substring(den.cols,1,7),'FF',sep='') else {
    den.cols <- if (length(den.cols[0])>7) my.cols else paste(substring(den.cols,1,7),'BF',sep='')
  }
  if (is.na(lwd.den)) {
    lwd.den <- if ("lwd" %in% nmdots) dots$lwd
    else par('lwd')
  }
  if (grepl('x', margin.xy)) {
    divider = max(sapply(density.x, function(x) max(x$y)))
    for (i in 1:length(ordered_keys)) {
      points(density.x[[i]]$x, density.x[[i]]$y/divider*multiplier.x+ylim[1], type='l', col=den.cols[i], lwd=lwd.den, ...)
    }
  }
  if (grepl('y', margin.xy)) {
    divider = max(sapply(density.y, function(x) max(x$y)))
    for (i in 1:length(ordered_keys)) {
      points(density.y[[i]]$y/divider*multiplier.y+xlim[1], density.y[[i]]$x, type='l', col=den.cols[i], lwd=lwd.den, ...)
    }
  }
  if (label.it) {
    d.f.sub <- subset(d.f, name != '')
    text(d.f.sub$x, d.f.sub$y, d.f.sub$name, pos=d.f.sub$pos, cex=0.6)
  }
}

rose24h <- function(x, color, title='', ylim=NULL) {
  ggplot(data=data.frame(x=x), aes(x=x)) +
    geom_histogram(breaks=seq(0,24), color='black',fill=color) +
    coord_polar(start=0) +
    scale_y_continuous("Count", limits=ylim) +
    scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0,24)) +
    theme_minimal() +
    theme(panel.grid.major=element_line(colour = 'darkgrey'), axis.text=element_text(size = 14), plot.title=element_text(face='bold', size=16)) +
    ggtitle(title)
}

motion_scatter_with_densities  <- function(d.f, nframes, ordered_keys, my.cols, margin.xy='xy', multiplier = 0.8,
                                           main='', xlab='', ylab='', xlim=c(), ylim=c(), log.xy='xy', den.solid=T, pch=21, ...) {
  
  tpoints <- ncol(d.f$x)
  for (frame in 1:nframes -1) {
    tiff(paste('f',sprintf("%04d", frame+1),'.tif', sep=''), width=480)
    d.f.frame <- motion_smoother(d.f, frame, nframes, tpoints)
    scatter_with_marginal_densities(d.f.frame, ordered_keys, my.cols, margin.xy, multiplier, main, xlab, ylab,
                                    xlim, ylim, log.xy, den.solid, pch, ...)
    #smooth_text(x = 4.5, y=3.5, frame, nframes, tpoints, tlabels=paste('ZT', seq(0,22,2), sep=''), sep='')
    text_pos <- smooth_clock(x=4.5, y=3, frame, nframes, tpoints, radius=0.25, radius.y=0.5, lwd.clock=2)
    smooth_text(text_pos$x, text_pos$y, frame, nframes, tpoints,  tlabels=paste('ZT', seq(0,22,2), sep=''))
    dev.off()
  }
}

# ========= Functions that depend heavily on the analysis.R data structures

get_mirids <- function (mir_assoc, mir_list = c(), uniq=T) {
  ids <- mir_assoc$GeneID[which(mir_assoc$miR %in% mir_list)]
  if (uniq) unique(ids) else ids
}

highlight_mirs <- function (rpkm.cds.tr, log.eff, expressed.list, mir.list, labels=c('All', 'miR targets'), order='last', margins='') {
  bcols = brewer.pal(n = 4, 'Paired')
  mcols = paste(bcols, '60', sep='')
  tdf.mir <- data.frame(x=rowMeans(no.inf(log10(rpkm.cds.tr[expressed.list,])), na.rm = T),
                        y=rowMeans(log.eff[expressed.list,], na.rm=T), k='All')
  tdf.mir$k <- as.character(tdf.mir$k)
  mir.ids <- get_mirids(mir_assoc, mir.list)
  tdf.mir[mir.ids[which(mir.ids %in% expressed.list)],]$k <- 'mir-targets'
  print(summary(mir.ids %in% expressed.list))
  if (order=='last') {
    plot_term <- c('All', 'mir-targets')
  } else {
    plot_term <- c('mir-targets', 'All')
  }
  if (grepl('x', margins)) xlim = c(-2,5) else xlim = c(-1,5)
  if (grepl('y', margins)) ylim = c(-6.5,4) else ylim = c(-5,4)
  scatter_with_marginal_densities(tdf.mir, plot_term, margin.xy = margins, log.xy = '',
                                  xlim=xlim, ylim=ylim, my.cols = mcols[c(1,4)] , pch=21,
                                  main='Ribosome occupancy vs transcript adundance',
                                  xlab='Transcript abundance, log10(RPKM)', ylab='log2(Ribosome occupancy)', lwd=2)
  legend('topright', labels, col=NA, pt.bg=bcols[c(1,4)], pch=21, cex=1.3, lty=0, bty='n')
}

getPCA = function(x, intgroup, ntop=500)
{
  require("lattice")
  require("genefilter")
  rv = rowVars(exprs(x))
  select = order(rv, decreasing=TRUE)[seq_len(ntop)]
  prcomp(t(exprs(x)[select,]))
}

getPCA2 = function(x, intgroup, ntop=500)
{
  require("lattice")
  require("genefilter")
  rv = rowVars(exprs(x))
  select = order(rv, decreasing=TRUE)[seq_len(ntop)]
  princomp(exprs(x)[select,])
}

add_semi_boxfunc <- function (s3d.obj) {
  s3d.obj$semibox <- function (...) {
    mem.par <- par(mar = mar, usr = usr)
    on.exit(par(mem.par))
    lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, ...)
    lines(c(x.min, x.max)+y.max*yx.f, c(z.max, z.max)+y.max*yz.f, ...)
    lines(c(x.max, x.max)+y.max*yx.f, c(z.max, z.min)+y.max*yz.f, ...)
    lines(c(x.min, x.min)+y.max*yx.f, c(z.max, z.min)+y.max*yz.f, ...)
  }
  environment(s3d.obj$semibox) <-environment(s3d.obj$box3d)
  s3d.obj
}

add_grey_plane <- function(s3d.obj) {
  s3d.obj$grey_plane <- function(...) {
    mem.par <- par(mar = mar, usr = usr)
    on.exit(par(mem.par))
    polygon(c(x.min, x.max, x.max+y.max*yx.f, x.min+y.max*yx.f),
            c(z.min, z.min, z.min+y.max*yz.f, z.min+y.max*yz.f), ...)
  }
  environment(s3d.obj$grey_plane) <-environment(s3d.obj$box3d)
  s3d.obj
}

pretty_grid <- function (distx, disty=distx, col="lightgray", lty="dotted", lwd=par('lwd')) {
  usr <- floor(par('usr'))
  atx <- c(usr[1]:usr[2])[usr[1]:usr[2] %% distx == 0]
  aty <- c(usr[3]:usr[4])[usr[3]:usr[4] %% disty == 0]
  abline(v=atx, col=col, lty=lty, lwd=lwd)
  abline(h=aty, col=col, lty=lty, lwd=lwd)
}

diagonal_grid <- function(slope=-1, interval=2, col="lightgray", lty="dotted", lwd=par('lwd')) {
  usr <- floor(par('usr'))
  atx <- c(usr[1]:usr[2])[usr[1]:usr[2] %% interval == 0]
  aty <- c(usr[3]:usr[4])[usr[3]:usr[4] %% interval == 0]
  a1 <- max(aty) - slope * min(atx)
  a2 <- min(aty) - slope * max(atx)
  a3 <- min(aty) - slope * min(atx)
  a4 <- max(aty) - slope * max(atx)
  for (a in seq(min(a1,a2,a3,a4), max(a1,a2,a3,a4), interval)) abline(a, slope, col=col, lty=lty, lwd=lwd)
}

col2hex <- function(col, maxcol=255) {
  c <- col2rgb(col)
  rgb(c[1],c[2],c[3], maxColorValue=maxcol)
}
