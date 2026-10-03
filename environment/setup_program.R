# Restore/check exact packages only into an explicitly supplied isolated library.
args <- commandArgs(trailingOnly=TRUE)
at <- match('--library',args)
if(is.na(at) || at==length(args)) stop('Usage: Rscript environment/setup_program.R --library /isolated/library [--install]')
if(as.character(getRversion())!='4.3.2') stop('Use R 4.3.2')
lib <- path.expand(args[at+1]);dir.create(lib,recursive=TRUE,showWarnings=FALSE)
lib <- normalizePath(lib);.libPaths(c(lib,.libPaths()))
script <- sub('--file=','',grep('--file=',commandArgs(FALSE),value=TRUE))
lock <- read.delim(file.path(dirname(script),'program-packages.tsv'),stringsAsFactors=FALSE)
version <- function(p) if(requireNamespace(p,quietly=TRUE)) packageDescription(p)$Version else NA_character_
missing <- lock$Package[vapply(seq_len(nrow(lock)),function(i) !identical(version(lock$Package[i]),lock$Version[i]),TRUE)]
if(length(missing) && !('--install' %in% args)) stop(paste('Missing/mismatched locked packages:',paste(missing,collapse=', ')))
options(timeout=300)
# Dependency order is explicit for the eight-package hard-dependency closure.
for(p in c('Rcpp','locfit','lattice','statmod','jsonlite','digest','limma','edgeR')) {
  if(!(p %in% missing)) next
  row <- lock[lock$Package==p,];filename <- paste0(p,'_',row$Version,'.tar.gz')
  urls <- if(row$Source=='CRAN') c(paste0('https://cran.r-project.org/src/contrib/',filename),
    paste0('https://cran.r-project.org/src/contrib/Archive/',p,'/',filename)) else
    c(paste0('https://bioconductor.org/packages/3.18/bioc/src/contrib/',filename),
      paste0('https://bioconductor.org/packages/3.18/bioc/src/contrib/Archive/',p,'/',filename))
  downloaded <- FALSE;archive <- tempfile(fileext='.tar.gz')
  for(url in urls) if(isTRUE(tryCatch({download.file(url,archive,quiet=TRUE);TRUE},error=function(e) FALSE))) {downloaded<-TRUE;break}
  if(!downloaded) stop(paste('Cannot download exact recorded package:',p,row$Version))
  install.packages(archive,lib=lib,repos=NULL,type='source');unlink(archive)
  if(!identical(version(p),row$Version)) stop(paste('Installation version mismatch:',p))
}
for(i in seq_len(nrow(lock))) if(!identical(version(lock$Package[i]),lock$Version[i])) stop(paste('Version mismatch:',lock$Package[i]))
cat('PASS:',nrow(lock),'program packages in',lib,'\n')
