
 # Run with R version 3.0.2 or earlier because that is the version
 # required in the DESCRIPTION file R package sped (where these data
 # are going

 foo <- getRversion()
 foo <- as.character(foo)
 foo
 if (utils::compareVersion(foo, "3.2.0") > 0)
     stop("must use R version 3.2.0 or earlier")

 # Figure 1 of Thompson (1986) cited in help pages and design doc

 thompson <- matrix(c("W", "R", "V",
                      "V", "T", "U",
                      "U", "I", "M",
                      "T", "J", "K",
                      "R", "N", "O",
                      "Q", "P", "O",
                      "P", "D", "L",
                      "O", "F", "G",
                      "N", "F", "H",
                      "M", "F", "H",
                      "J", "B", "C",
                      "I", "B", "C",
                      "H", "B", "C",
                      "G", "B", "A",
                      "D", "B", "A"),
                      ncol = 3, byrow = TRUE)
 colnames(thompson) <- c("ind", "pa", "ma")

 destdir <- "../package/sped/data"

 save(thompson, file = file.path(destdir, "thompson.rda"))

 nrow(thompson)
 length(sort(unique(as.vector(thompson))))


