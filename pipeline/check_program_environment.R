# Check only; no global or automatic package installation.
extra <- Sys.getenv('PROGRAM_R_LIBRARY','')
if(nzchar(extra)) .libPaths(c(extra,.libPaths()))
if(as.character(getRversion())!='4.3.2') stop('Program reproduction requires recorded R 4.3.2')
script <- sub('--file=','',grep('--file=',commandArgs(FALSE),value=TRUE))
lockfile <- file.path(dirname(script),'environment/program-packages.tsv')
if(!file.exists(lockfile)) lockfile <- file.path(dirname(dirname(script)),'environment/program-packages.tsv')
lock <- read.delim(lockfile,stringsAsFactors=FALSE)
for(i in seq_len(nrow(lock))) {
  p <- lock$Package[i]
  if(!requireNamespace(p,quietly=TRUE)) stop(paste('Required program package missing:',p))
  if(!identical(packageDescription(p)$Version,lock$Version[i]))
    stop(paste('Recorded program dependency version required:',p,lock$Version[i]))
}
cat('PASS: R 4.3.2 and',nrow(lock),'locked program dependencies\n')
