#!/usr/bin/env Rscript
# Restore the tested reference library; only --install authorizes installation.
args <- commandArgs(trailingOnly=TRUE)
at <- match('--library',args)
if (is.na(at) || at==length(args)) stop('Usage: Rscript environment/setup_reference.R --library /new/library [--install]')
if (as.character(getRversion())!='4.3.2') stop('Install/use R 4.3.2 before restoring this library')
lib <- path.expand(args[[at+1]])
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
lib <- normalizePath(lib)
.libPaths(c(lib,.libPaths()))
script <- sub('--file=','',grep('--file=',commandArgs(FALSE),value=TRUE))
lock <- read.delim(file.path(dirname(script),'reference-packages.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
version <- function(p) tryCatch(packageDescription(p)$Version,error=function(e) NA_character_)
remaining <- lock$Package[is.na(vapply(lock$Package,version,'')) | vapply(lock$Package,version,'') != lock$Version]
remaining <- remaining[!is.na(remaining)]
if(length(remaining) && !('--install' %in% args)) stop(paste('Missing/mismatched locked packages:',paste(remaining,collapse=', '),'; add --install to restore into the specified library'))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=300)
if(length(remaining)) {
  cat('Install plan: restore',length(remaining),'exact package versions into',lib,'\n')
  binary <- .Platform$pkgType
  available <- available.packages(type=binary)
  for(p in remaining) {
    wanted <- lock$Version[match(p,lock$Package)]
    if(p %in% rownames(available) && available[p,'Version']==wanted) {
      install.packages(p,lib=lib,type=binary,dependencies=FALSE)
    }
  }
  remaining <- remaining[vapply(remaining,function(p) is.na(version(p)) || version(p)!=lock$Version[match(p,lock$Package)],TRUE)]
  # Install archived sources after their locked hard dependencies. Base packages
  # come with R; recommended packages are restored like other CRAN packages.
  dependencies <- function(p) {
    row <- lock[lock$Package==p,]
    text <- paste(row$Depends,row$Imports,row$LinkingTo,sep=',')
    names <- trimws(gsub('\\s*\\([^)]*\\)','',strsplit(text,',')[[1]]))
    intersect(names,lock$Package)
  }
  while(length(remaining)) {
    ready <- remaining[vapply(remaining,function(p) !any(dependencies(p) %in% remaining),TRUE)]
    if(!length(ready)) stop(paste('Cannot resolve locked dependency order:',paste(remaining,collapse=', ')))
    for(p in ready) {
      wanted <- lock$Version[match(p,lock$Package)]
      current <- available.packages(type='source')
      url <- if(p %in% rownames(current) && current[p,'Version']==wanted)
        paste0('https://cran.r-project.org/src/contrib/',p,'_',wanted,'.tar.gz') else
        paste0('https://cran.r-project.org/src/contrib/Archive/',p,'/',p,'_',wanted,'.tar.gz')
      install.packages(url,lib=lib,repos=NULL,type='source')
      if(is.na(version(p)) || version(p)!=wanted) stop(paste('Failed locked installation:',p,wanted))
    }
    remaining <- setdiff(remaining,ready)
  }
}
for(i in seq_len(nrow(lock))) if(is.na(version(lock$Package[i])) || version(lock$Package[i])!=lock$Version[i]) stop(paste('Version mismatch:',lock$Package[i]))
cat('PASS:',nrow(lock),'locked R packages; library',lib,'\n')
