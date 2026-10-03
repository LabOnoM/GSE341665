# Shared production functions, also exercised on synthetic data by CI.
construct_length_scaled_counts <- function(tpm, effective_length, raw_counts) {
  if (!identical(dimnames(tpm),dimnames(effective_length)) ||
      !identical(dimnames(tpm),dimnames(raw_counts))) stop('Abundance/count/length alignment mismatch')
  if (any(!is.finite(tpm)) || any(tpm<0) || any(!is.finite(effective_length)) ||
      any(effective_length<=0) || any(!is.finite(raw_counts)) || any(raw_counts<0))
    stop('Invalid abundance, length or count values')
  scaled <- tpm * rowMeans(effective_length)
  if (any(colSums(scaled)<=0) || any(colSums(raw_counts)<=0)) stop('Empty library')
  scaled <- sweep(scaled,2,colSums(scaled),'/')
  scaled <- sweep(scaled,2,colSums(raw_counts),'*')
  if (max(abs(colSums(scaled)/colSums(raw_counts)-1))>=1e-12) stop('Library totals changed')
  scaled
}

map_programs <- function(programs, symbols) {
  if (anyDuplicated(unlist(programs))) stop('Duplicate program gene')
  lapply(programs,function(members) {
    hits <- lapply(members,function(g) which(symbols==g))
    if (any(lengths(hits)!=1L)) stop('Program genes must map uniquely to retained genes')
    unlist(hits,use.names=FALSE)
  })
}

run_gene_set_tests <- function(v, index, design, seed, nrot) {
  if (qr(design)$rank!=ncol(design)) stop('Design is not full rank')
  if (nrow(design)-ncol(design)<3L) stop('Insufficient residual degrees of freedom')
  if (!identical(colnames(design)[ncol(design)],'genotypedcKO')) stop('Contrast sign/column mismatch')
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  set.seed(seed)
  rr <- limma::mroast(v,index=index,design=design,contrast=ncol(design),nrot=nrot,
                     set.statistic='mean',adjust.method='BH',midp=FALSE,sort='none')
  cc <- limma::camera(v,index=index,design=design,contrast=ncol(design),
                     inter.gene.cor=NA,allow.neg.cor=FALSE,use.ranks=FALSE,sort=FALSE)
  cf <- limma::camera(v,index=index,design=design,contrast=ncol(design),
                     inter.gene.cor=0.01,allow.neg.cor=FALSE,use.ranks=FALSE,sort=FALSE)
  for (result in list(rr,cc,cf))
    if (!isTRUE(all.equal(result$FDR,p.adjust(result$PValue,'BH')))) stop('BH family mismatch')
  list(ROAST=rr,CAMERA=cc,CAMERA_fixed001=cf)
}
