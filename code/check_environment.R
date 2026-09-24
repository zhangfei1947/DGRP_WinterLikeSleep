#!/usr/bin/env Rscript
packages <- c('data.table','damr','behavr','sleepr','ggplot2','ggrepel','patchwork','lme4')
x <- data.frame(package=packages,installed=vapply(packages,requireNamespace,logical(1),quietly=TRUE))
x$version <- vapply(packages,function(p) if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else NA_character_,character(1))
print(x,row.names=FALSE)
print(R.version.string)
if(!all(x$installed)) quit(status=1)
