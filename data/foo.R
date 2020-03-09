
 source("dumpdata")

 ls()

 alberta <- Alberta.triplets
 colnames(alberta) <- c("ind", "pa", "ma")
 alberta <- alberta[! is.na(alberta[ , "pa"]), ]

 destdir <- "../package/sped/data"

 save(alberta, file = file.path(destdir, "alberta.rda"))

 nrow(alberta)
 length(sort(unique(as.vector(alberta))))
