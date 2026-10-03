args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4L || (length(args)==6L && args[5]=='--geo-root'))
script <- sub('--file=','',grep('--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'program_model_functions.R'))
data_root <- normalizePath(args[1]); out <- args[2]; program_file <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
extra_lib <- Sys.getenv('PROGRAM_R_LIBRARY','')
if(nzchar(extra_lib)) .libPaths(c(extra_lib,.libPaths()))
suppressPackageStartupMessages({library(limma);library(edgeR);library(jsonlite)})
options(digits=16)
cfg <- fromJSON(args[4])
stopifnot(cfg$retained_genes==15803L,cfg$ROAST$rotations==99999L,
          cfg$ROAST$primary_seed==20261003L,cfg$ROAST$sensitivity_seed==20261004L,
          identical(cfg$ROAST$midp,FALSE),cfg$contrast=='dcKO minus CON')
input_path <- file.path(data_root,'04_pipeline_objects/DESeq2_Processing.Rdata')
hash_file <- function(p) digest::digest(file=p,algo='sha256')
provenance <- list(original_DESeq2_refit=FALSE,program_limma_voom_fit=TRUE,
    sources=list(list(path='04_pipeline_objects/DESeq2_Processing.Rdata',sha256=hash_file(input_path),bytes=file.info(input_path)$size)),
    configuration=cfg,configuration_sha256=hash_file(args[4]),membership_json_sha256=hash_file(program_file),
    GEO_count_crosscheck=list(performed=FALSE),reconstructed_from='DEGs_dds assays counts/avgTxLength; DEGs_Out_TPM full-precision TPM/raw counts; colData RUV_x/W_1')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(capture.output(list(mroast=formals(mroast.default),camera=formals(camera.default))),file.path(out,'function_defaults.txt'))
e <- new.env()
load(file.path(data_root,'04_pipeline_objects/DESeq2_Processing.Rdata'),envir=e)
# Extract serialized numeric data without refitting or changing the DESeq2 object.
dds <- e$DEGs_dds
al <- attr(attr(attr(dds,'assays'),'data'),'listData')
counts <- al$counts; lengths <- al$avgTxLength; nf <- al$normalizationFactors
cd <- attr(attr(dds,'colData'),'listData')
samples <- colnames(counts); ids <- rownames(counts)
stopifnot(identical(samples,c('CON1','CON2','CON3','dcKO1','dcKO2','dcKO3')),
          nrow(counts)==15803L,all(is.finite(counts)),all(counts>=0),
          !anyDuplicated(samples),!anyDuplicated(ids),
          identical(dimnames(lengths),dimnames(counts)),all(lengths>0))
tab <- e$DEGs_Out_TPM[match(ids,e$DEGs_Out_TPM$ensembl_gene_id),]
stopifnot(identical(as.character(tab$ensembl_gene_id),ids))
sx <- c('CON_1','CON_2','CON_3','dcKO_1','dcKO_2','dcKO_3')
tpm <- as.matrix(tab[paste0('TPM_',sx)])
raw_counts <- as.matrix(tab[paste0('RawCounts_',sx)])
dimnames(tpm)<-dimnames(raw_counts)<-dimnames(counts)
stopifnot(all(is.finite(tpm)),all(tpm>=0),max(abs(round(raw_counts)-counts))==0)
stopifnot(max(abs(as.matrix(tab[paste0('FeatureLength_',sx)])-lengths))<1e-8)
if(length(args)==6L) {
geo_file <- file.path(args[6],'GSE_processed_gene_raw_counts.txt')
geo <- read.delim(geo_file,check.names=FALSE)
gcounts <- as.matrix(geo[match(ids,geo$ensembl_gene_id),samples])
stopifnot(!anyNA(gcounts),max(abs(gcounts-counts))==0)
provenance$GEO_count_crosscheck <- list(performed=TRUE,cells=length(counts),sha256=hash_file(geo_file))
}
genes <- data.frame(ensembl_gene_id=ids,symbol=as.character(tab$external_gene_name))
programs <- fromJSON(program_file,simplifyVector=TRUE)
index <- map_programs(programs,genes$symbol)
stopifnot(identical(as.integer(lengths(index)),c(16L,20L,14L,18L)))
membership <- do.call(rbind,lapply(names(index),function(k) data.frame(program=k,genes[index[[k]],])))
write.csv(membership,file.path(out,'program_membership.csv'),row.names=FALSE)
metadata <- data.frame(sample=samples,genotype=factor(cd$RUV_x,levels=c('CON','dcKO')),W_1=cd$W_1,row.names=samples)
stopifnot(identical(as.character(metadata$genotype),rep(c('CON','dcKO'),each=3)),all(is.finite(metadata$W_1)))
designs <- list(RUV_adjusted=model.matrix(~W_1+genotype,metadata),unadjusted=model.matrix(~genotype,metadata))
stopifnot(max(abs(unname(designs$RUV_adjusted)-unname(attr(dds,'modelMatrix'))))==0)
write.csv(metadata,file.path(out,'sample_metadata.csv'),row.names=FALSE)
# tximport lengthScaledTPM construction, restricted to the fixed 15,803-gene universe.
# Preserve abundance patterns, correct sample-specific effective length, retain count scale.
length_scaled <- construct_length_scaled_counts(tpm,lengths,raw_counts)
y <- calcNormFactors(DGEList(counts=length_scaled,genes=genes),method='TMM')
saveRDS(list(raw_counts=raw_counts,dds_counts=counts,TPM=tpm,avgTxLength=lengths,
             archived_normalizationFactors=nf,lengthScaledTPM=length_scaled,
             TMM=y$samples,metadata=metadata,designs=designs,programs=programs,index=index),
        file.path(out,'analysis_inputs.rds'))
libs <- data.frame(sample=samples,raw_count_total=colSums(raw_counts),
                   detected=rowSums(t(raw_counts)>0),TMM_factor=y$samples$norm.factors,
                   effective_library_size=y$samples$lib.size*y$samples$norm.factors)
write.csv(libs,file.path(out,'library_qc.csv'),row.names=FALSE)
pdf(file.path(out,'QC.pdf'),width=9,height=6,onefile=TRUE)
barplot(colSums(raw_counts)/1e6,names.arg=samples,col=rep(c('#377EB8','#E66101'),each=3),ylab='Estimated counts (millions)',main='Retained 15,803-gene count totals')
logcpm <- cpm(y,log=TRUE,prior.count=2)
pc <- prcomp(t(logcpm),center=TRUE,scale.=FALSE)
pv <- 100*pc$sdev^2/sum(pc$sdev^2)
plot(pc$x[,1],pc$x[,2],col=rep(c('#377EB8','#E66101'),each=3),pch=19,cex=1.5,
     xlim=extendrange(pc$x[,1],f=.08),ylim=extendrange(pc$x[,2],f=.12),
     xlab=sprintf('PC1 (%.1f%%)',pv[1]),ylab=sprintf('PC2 (%.1f%%)',pv[2]),main='Unadjusted logCPM sample PCA')
text(pc$x[,1],pc$x[,2],labels=samples,pos=3,cex=.8)
write.csv(data.frame(sample=samples,pc$x),file.path(out,'sample_PCA.csv'),row.names=FALSE)
write.csv(cor(logcpm),file.path(out,'sample_correlations.csv'))
all_results <- list(); fit_models<-list(); qc<-list()
for(i in seq_along(designs)) {
 label<-names(designs)[i]; design<-designs[[i]]; coef<-ncol(design)
 stopifnot(qr(design)$rank==ncol(design),nrow(design)-ncol(design)>=3)
 write.csv(design,file.path(out,paste0('design_',label,'.csv')))
 v<-voom(y,design,plot=FALSE,save.plot=TRUE)
 plot(v$voom.xy$x,v$voom.xy$y,pch=16,cex=.25,col='#00000055',
      xlab=v$voom.xy$xlab,ylab=v$voom.xy$ylab,main=paste('voom mean-variance trend:',label))
 lines(v$voom.line$x,v$voom.line$y,col='#D73027',lwd=2)
 fit<-eBayes(lmFit(v,design),robust=FALSE,trend=FALSE)
 residual<-v$E-fit$coefficients%*%t(design)
 boxplot(as.data.frame(residual),outline=FALSE,las=2,ylab='logCPM residual',main=paste('Model residuals:',label))
 qc[[label]]<-list(residual_df=nrow(design)-ncol(design),condition_number=kappa(design),
                  W1_genotype_correlation=cor(metadata$W_1,as.numeric(metadata$genotype)),
                  residual_RMS=sqrt(colMeans(residual^2)),weight_range=range(v$weights),
                  empirical_Bayes_df_prior=fit$df.prior,
                  weight_quantiles=quantile(v$weights,c(0,.25,.5,.75,.95,.99,1)),
                  all_weights_finite_positive=all(is.finite(v$weights) & v$weights>0),
                  residual_df_per_gene=range(fit$df.residual),
                  full_rank=qr(design)$rank==ncol(design))
 gtest<-topTable(fit,coef=coef,number=Inf,sort.by='none')
 contributions <- do.call(rbind,lapply(names(index),function(k) {
   z <- index[[k]]
   data.frame(program=k,genes[z,],logFC=fit$coefficients[z,coef],moderated_t=fit$t[z,coef],
              mean_observation_weight=rowMeans(v$weights)[z],direction=ifelse(fit$coefficients[z,coef]>0,'Up','Down'))
 }))
 write.csv(contributions,file.path(out,paste0('program_gene_contributions_',label,'.csv')),row.names=FALSE)
 write.csv(gtest,file.path(out,paste0('gene_statistics_',label,'.csv')),row.names=FALSE)
 cat('Running',label,'ROAST with 99,999 rotations\n');flush.console()
 tests <- run_gene_set_tests(v,index,design,20261003+i-1L,99999L)
 rr <- tests$ROAST;cc <- tests$CAMERA;cf <- tests$CAMERA_fixed001
 for(n in c('rr','cc','cf')) {
   nm<-switch(n,rr='ROAST',cc='CAMERA_estimated_correlation',cf='CAMERA_fixed_correlation_001')
   result<-get(n); write.csv(data.frame(Program=rownames(result),result,row.names=NULL),
                             file.path(out,paste0(nm,'_',label,'.csv')),row.names=FALSE)
 }
 result<-data.frame(Model=label,Program=names(index),NGenes=lengths(index),
                    ROAST_direction=rr$Direction,ROAST_P=rr$PValue,ROAST_FDR=rr$FDR,
                    ROAST_mixed_P=rr$PValue.Mixed,ROAST_mixed_FDR=rr$FDR.Mixed,
                    CAMERA_direction=cc$Direction,CAMERA_correlation=cc$Correlation,
                    CAMERA_P=cc$PValue,CAMERA_FDR=cc$FDR,
                    CAMERA_fixed001_P=cf$PValue,CAMERA_fixed001_FDR=cf$FDR,row.names=NULL)
 print(result);all_results[[label]]<-result
 fit_models[[label]]<-list(voom=v,fit=fit,ROAST=rr,CAMERA=cc,CAMERA_fixed001=cf)
}
dev.off()
write.csv(do.call(rbind,all_results),file.path(out,'combined_results.csv'),row.names=FALSE)
saveRDS(fit_models,file.path(out,'fitted_models.rds'))
write_json(qc,file.path(out,'model_qc.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
provenance$RNG <- RNGkind()
provenance$sample_count <- length(samples);provenance$retained_genes <- length(ids)
provenance$program_genes_unique <- length(unique(unlist(programs)))
write_json(provenance,file.path(out,'input_provenance.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
write_json(cfg,file.path(out,'analysis_config.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('COMPLETED\n')
