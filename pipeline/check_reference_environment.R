#!/usr/bin/env Rscript
# Check the recorded reference environment before recomputing annotations.
expected <- c(Seurat="5.3.1",SeuratObject="5.2.0",Matrix="1.6.5",
              ggplot2="4.0.0",patchwork="1.3.2",uwot="0.2.3")
if (as.character(getRversion()) != "4.3.2") stop("Reference reproduction requires recorded R 4.3.2")
for (p in names(expected)) {
  if (!requireNamespace(p,quietly=TRUE)) stop(paste("Required R package missing:",p))
  if (as.character(packageVersion(p)) != expected[[p]]) stop(paste("Recorded version required:",p,expected[[p]]))
}
cat("PASS: recorded R and six reference package versions\n")
here <- dirname(sub('--file=','',grep('--file=',commandArgs(FALSE),value=TRUE)))
lockfile <- file.path(here,'environment/reference-packages.tsv')
if(!file.exists(lockfile)) lockfile <- file.path(dirname(here),'environment/reference-packages.tsv')
lock <- read.delim(lockfile,stringsAsFactors=FALSE)
for(i in seq_len(nrow(lock))) {
  if(!identical(packageDescription(lock$Package[i])$Version,lock$Version[i]))
    stop(paste('Locked dependency version required:',lock$Package[i],lock$Version[i]))
}
cat('PASS:',nrow(lock),'tested R dependency versions\n')
