# Shared definitions for the profile discovery pilot (8L:16D, one-minute DAM).
profile_version <- "2026-09-09-v1"
safe_cor <- function(x, y, method = "pearson", min_n = 8L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_n || sd(x[ok]) < 1e-10 || sd(y[ok]) < 1e-10) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}
finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
finite_median <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
strict_mean <- function(x, required) if (length(x) == required && all(is.finite(x))) mean(x) else NA_real_
periodic_smooth <- function(x) (x[c(length(x), seq_len(length(x)-1))] + 2*x + x[c(2:length(x),1)]) / 4
profile_features <- function(y) {
  # Bin centres; local phases are offsets from the corresponding LD transition.
  nm <- c("Harmonic24Amplitude", "Harmonic12Amplitude", "TwoHarmonicR2",
          "MorningPeakOffset", "EveningPeakOffset", "MorningPeakAmplitude", "EveningPeakAmplitude")
  out <- setNames(rep(NA_real_, length(nm)), nm)
  if (length(y) != 48L || any(!is.finite(y))) return(out)
  h <- seq(.25, 23.75, by=.5)
  X <- cbind(1, cos(2*pi*h/24), sin(2*pi*h/24), cos(4*pi*h/24), sin(4*pi*h/24))
  b <- as.numeric(crossprod(X,y)/colSums(X^2))
  out[1:2] <- c(sqrt(sum(b[2:3]^2)),sqrt(sum(b[4:5]^2)))
  ss <- sum((y-mean(y))^2)
  if (ss > 1e-12) out[3] <- max(0,1-sum((y-X%*%b)^2)/ss)
  sm <- periodic_smooth(y)
  for (j in 1:2) {
    centre <- c(0,8)[j]
    offset <- ((h-centre+12) %% 24)-12
    ix <- which(offset >= -3 & offset < 3)
    ix <- ix[order(offset[ix])]
    v <- sm[ix]; peak <- which.max(v); amp <- max(v)-min(v)
    out[5+j] <- amp
    # Prespecified heuristic for identifiable local peaks, not a rhythm test.
    if (amp >= .10 && peak > 1 && peak < length(ix) && sum(abs(v-max(v))<1e-8)==1L)
      out[3+j] <- offset[ix[peak]]
  }
  out
}
stability_features <- function(M) {
  # Rows are complete days D2-D6, columns are 48 half-hour bins.
  nm <- c("LD24hLagCorrelation","LDInterdailyStability","LDIntradailyVariability","OddEvenProfileCorrelation")
  out <- setNames(rep(NA_real_,4),nm)
  if (!all(dim(M)==c(5L,48L)) || any(!is.finite(M))) return(out)
  D <- sweep(M,1,rowMeans(M),"-")
  out[1] <- safe_cor(as.vector(t(D[1:4,])),as.vector(t(D[2:5,])))
  x <- as.vector(t(M)); ss <- sum((x-mean(x))^2)
  if(ss>1e-12) {
    out[2] <- 5*sum((colMeans(M)-mean(x))^2)/ss
    out[3] <- mean(diff(x)^2)/mean((x-mean(x))^2)
  }
  out[4] <- safe_cor(colMeans(M[c(1,3,5),]),colMeans(M[c(2,4),]))
  out
}
batch_adjust <- function(z, response="value", min_fit=6L) {
  z <- data.table::copy(z)
  z <- z[is.finite(get(response)) & n_flies >= min_fit]
  z[, Batch:=factor(Batch)]; z[, Genotype:=factor(Genotype)]
  if(nlevels(z$Batch)<2 || nlevels(z$Genotype)<2) return(NULL)
  X <- model.matrix(~Batch+Genotype, z)
  fit <- lm.fit(X,z[[response]])
  if(fit$rank<ncol(X)) return(NULL) # disconnected designs cannot calibrate globally
  bc <- grepl("^Batch",colnames(X))
  be <- as.vector(X[,bc,drop=FALSE]%*%fit$coefficients[bc])
  effects <- unique(data.table::data.table(Batch=as.character(z$Batch),effect=be))
  effects[,effect:=effect-mean(effect)]
  z[,adjusted:=get(response)-be+mean(unique(data.table::data.table(Batch=Batch,be=be))$be)]
  # Equal-batch EMM: intercept + genotype coefficient + mean(batch coefficients).
  gm <- unique(z[,.(Genotype=as.character(Genotype))])
  gm[,value:=fit$coefficients[1]+mean(c(0,fit$coefficients[bc]))]
  gc <- grepl("^Genotype",colnames(X))
  gcoef <- setNames(fit$coefficients[gc],sub("^Genotype","",colnames(X)[gc]))
  gm[Genotype %in% names(gcoef),value:=value+unname(gcoef[Genotype])]
  list(lines=gm,batches=effects,rank=fit$rank,residual_df=nrow(X)-fit$rank)
}
replicate_pairs <- function(z, min_n=8L) {
  z <- data.table::copy(z[Genotype!="CS" & is.finite(value) & n_flies>=min_n])
  data.table::setorder(z,Genotype,Batch)
  z[, if(.N<2L) NULL else {
    ix <- combn(seq_len(.N),2L)
    .(Batch1=Batch[ix[1,]],Batch2=Batch[ix[2,]],x=value[ix[1,]],y=value[ix[2,]])
  },by=Genotype]
}
pair_rho <- function(p) {
  if (nrow(p)<4L) return(NA_real_)
  # Symmetric pairing makes the result independent of batch ordering.
  safe_cor(c(p$x,p$y),c(p$y,p$x),"spearman")
}
