
 library(sped)

 ped <- read.table("ped.txt", header = TRUE, stringsAsFactors = FALSE)
 ped <- as.matrix(ped)
 ped

 ped[5, 2:3] <- rev(ped[5, 2:3])
 ped

 gen <- c(goofy = 1)

 # should give error
 try(descent("spot", ped, gen, check.sex = TRUE))
