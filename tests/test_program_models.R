# Small synthetic test of the actual production functions; no restricted data.
script <- sub('--file=','',grep('--file=',commandArgs(FALSE),value=TRUE))
extra <- Sys.getenv('PROGRAM_R_LIBRARY','')
if(nzchar(extra)) .libPaths(c(extra,.libPaths()))
source(file.path(dirname(dirname(script)),'pipeline/program_model_functions.R'))
suppressPackageStartupMessages({library(edgeR);library(limma)})
set.seed(7);n <- 100
counts <- matrix(rnbinom(n*6,mu=100,size=10),n,6)
counts[1:10,4:6] <- counts[1:10,4:6]*4
dimnames(counts) <- list(paste0('gene',1:n),paste0('sample',1:6))
len <- counts*0+1000
tpm <- sweep(counts,2,colSums(counts),'/')*1e6
scaled <- construct_length_scaled_counts(tpm,len,counts)
stopifnot(max(abs(scaled-counts))<1e-10)
metadata <- data.frame(W_1=c(-1,1,0,-.5,.5,0),genotype=factor(rep(c('CON','dcKO'),each=3),levels=c('CON','dcKO')))
index <- map_programs(list(signal=paste0('gene',1:10),neutral=paste0('gene',11:20)),rownames(counts))
for(design in list(model.matrix(~W_1+genotype,metadata),model.matrix(~genotype,metadata))) {
  y <- calcNormFactors(DGEList(scaled));v <- voom(y,design,plot=FALSE)
  result <- run_gene_set_tests(v,index,design,11,199)
  stopifnot(result$ROAST['signal','Direction']=='Up',result$ROAST['signal','PValue']>0,
            identical(result$ROAST,run_gene_set_tests(v,index,design,11,199)$ROAST))
}
fails <- function(expression) stopifnot(inherits(try(expression,silent=TRUE),'try-error'))
fails(map_programs(list(ambiguous='gene1'),c('gene1','gene1')))
fails(map_programs(list(missing='absent'),rownames(counts)))
fails(map_programs(list(a='gene1',b='gene1'),rownames(counts)))
bad <- len;bad[1,1] <- 0;fails(construct_length_scaled_counts(tpm,bad,counts))
design <- model.matrix(~W_1+genotype,metadata);design[,2] <- design[,1]
fails(run_gene_set_tests(v,index,design,11,199))
cat('PASS: production synthetic inference, reproducible RNG, count construction and failure modes\n')
