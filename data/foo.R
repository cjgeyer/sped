
 # Run with R version 3.0.2 or earlier because that is the version
 # required in the DESCRIPTION file R package sped (where these data
 # are going

 foo <- getRversion()
 foo <- as.character(foo)
 foo
 if (utils::compareVersion(foo, "3.0.2") > 0)
     stop("must use R version 3.0.2 or earlier")

 source("dumpdata")

 ls()

 alberta <- Alberta.triplets
 colnames(alberta) <- c("ind", "pa", "ma")
 alberta <- alberta[! is.na(alberta[ , "pa"]), ]

 destdir <- "../package/sped/data"

 storage.mode(alberta) <- "integer"

 save(alberta, file = file.path(destdir, "alberta.rda"))

 nrow(alberta)
 length(sort(unique(as.vector(alberta))))
