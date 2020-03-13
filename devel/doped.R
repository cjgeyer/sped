
 library(sped, lib.loc = "../package/sped.Rcheck")

 ped <- read.table("ped.txt", header = TRUE, stringsAsFactors = FALSE)
 ped <- as.matrix(ped)
 ped

 gen <- list(goofy = 1)
 gen

 checkitout <- function(foo) {

     # case (a) empty argument list
     if (length(foo$individuals) == 0) {
         if ((foo$value != 1) || foo$type != "a")
             stop("error in case (a)")
         return(invisible("OK"))
     }
     # now we know length(foo$individuals) > 0

     # case (b) B_1 is not a founder, not repeated, and contains no genes of S
     if (foo$type == "b") {
         b1 <- foo$individuals[1]
         r <- sum(foo$individuals == b1)
         if (r != 1)
             stop("error in case (b): b1 repeated")
         if (b1 %in% names(gen))
             stop("error in case (b): b1 contains genes of S")
         b1ped <- ped[ped[,1] == b1, ]
         if (length(b1ped) == 0)
             stop("error in case (b): b1 is founder")
         if (is.matrix(b1ped))
             stop("error in case (b): b1 appears more than once in pedigree")
         pa <- b1ped[2]
         ma <- b1ped[3]
         if (length(foo$calls) != 2)
             stop("error in case (b): not 2 calls")
         indiv.rest <- foo$individuals[-1]
         if (foo$calls[[1]]$individuals != c(pa, indiv.rest))
             stop("error in case (b): wrong individuals for 1st call")
         if (foo$calls[[2]]$individuals != c(ma, indiv.rest))
             stop("error in case (b): wrong individuals for 2nd call")
         if (! all.equal(foo$value,
             (foo$calls[[1]]$value + foo$calls[[2]]$value) / 2))
             stop("error in case (b): wrong value")
         lapply(foo$calls, checkitout)
     }

     # case (c) B_1 is not a founder, is repeated, and contains no genes of S
     if (foo$type == "c") {
         b1 <- foo$individuals[1]
         r <- sum(foo$individuals == b1)
         if (r == 1)
             stop("error in case (c): b1 not repeated")
         if (b1 %in% names(gen))
             stop("error in case (c): b1 contains genes of S")
         b1ped <- ped[ped[,1] == b1, ]
         if (length(b1ped) == 0)
             stop("error in case (c): b1 is founder")
         if (is.matrix(b1ped))
             stop("error in case (c): b1 appears more than once in pedigree")
         pa <- b1ped[2]
         ma <- b1ped[3]
         if (length(foo$calls) != 2)
             stop("error in case (c): not 2 calls")
         indiv.rest <- foo$individuals[- (1:r)]
         if (foo$calls[[1]]$individuals != c(b1, indiv.rest))
             stop("error in case (c): wrong individuals for 1st call")
         if (foo$calls[[2]]$individuals != c(pa, ma, indiv.rest))
             stop("error in case (c): wrong individuals for 2nd call")
         half_to_r_minus_1 <- (1 / 2)^(r - 1)
         if (! all.equal(foo$value, half_to_r_minus_1 * foo$calls[[1]]$value +
             (1 - half_to_r_minus_1) * foo$calls[[2]]$value))
             stop("error in case (c): wrong value")
         lapply(foo$calls, checkitout)
     }

     # case (d) B_1 is a founder and contains no genes of S
     if (foo$type == "d") {
         b1 <- foo$individuals[1]
         if (b1 %in% names(gen))
             stop("error in case (d): b1 contains genes of S")
         b1ped <- ped[ped[,1] == b1, ]
         if (length(b1ped) != 0)
             stop("error in case (d): b1 is not a founder")
         if (! is.null(foo$calls))
             stop("error in case (d): not 0 calls")
         if (! all.equal(foo$value, 0))
             stop("error in case (d): wrong value")
     }

     # case (e) B_1 contains 2 genes of S
     if (foo$type == "e") {
         b1 <- foo$individuals[1]
         g1 <- gen[b1]
         if (g1 != 2)
             stop("error in case (e): b1 does not contain 2 genes of S")
         r <- sum(foo$individuals == b1)
         indiv.rest <- foo$individuals[- (1:r)]
         if (length(foo$calls) != 1)
             stop("error in case (e): not 1 call")
         if (foo$calls[[1]]$individuals != indiv.rest)
             stop("error in case (e): wrong individuals for 1st call")
         if (! all.equal(foo$value, foo$calls[[1]]$value))
             stop("error in case (e): wrong value")
         lapply(foo$calls, checkitout)
     }

     # case (f) B_1 is not founder and contains 1 gene of S
     if (foo$type == "f") {
         b1 <- foo$individuals[1]
         g1 <- gen[b1]
         if (g1 != 1)
             stop("error in case (f): b1 does not contain 1 gene of S")
         r <- sum(foo$individuals == b1)
         b1ped <- ped[ped[,1] == b1, ]
         if (length(b1ped) == 0)
             stop("error in case (f): b1 is founder")
         if (is.matrix(b1ped))
             stop("error in case (f): b1 appears more than once in pedigree")
         pa <- b1ped[2]
         ma <- b1ped[3]
         if (length(foo$calls) != 3)
             stop("error in case (f): not 3 calls")
         indiv.rest <- foo$individuals[- (1:r)]
         if (foo$calls[[1]]$individuals != indiv.rest)
             stop("error in case (g): wrong individuals for 1st call")
         if (foo$calls[[2]]$individuals != c(pa, indiv.rest))
             stop("error in case (g): wrong individuals for 2nd call")
         if (foo$calls[[3]]$individuals != c(ma, indiv.rest))
             stop("error in case (g): wrong individuals for 3rd call")
         half_to_r <- (1 / 2)^r
         if (! all.equal(foo$value, half_to_r * foo$calls[[1]]$value) +
             (1 / 2) * (1 - half_to_r) *
             (foo$calls[[2]]$value + foo$calls[[3]]$value))
             stop("error in case (g): wrong value")
         lapply(foo$calls, checkitout)
     }

     # case (g) B_1 is founder and contains 1 gene of S
     if (foo$type == "g") {
         b1 <- foo$individuals[1]
         g1 <- gen[b1]
         if (g1 != 1)
             stop("error in case (g): b1 does not contain 1 gene of S")
         r <- sum(foo$individuals == b1)
         indiv.rest <- foo$individuals[- (1:r)]
         if (length(foo$calls) != 1)
             stop("error in case (g): not 1 call")
         if (foo$calls[[1]]$individuals != indiv.rest)
             stop("error in case (g): wrong individuals for 1st call")
         half_to_r <- (1 / 2)^r
         if (! all.equal(foo$value, half_to_r * foo$calls[[1]]$value))
             stop("error in case (g): wrong value")
         lapply(foo$calls, checkitout)
     }

     stop("error: unknown case")
 }

 # case (a) empty individuals list

 descent(character(0), ped, gen, debug = TRUE)

 # case (b) B_1 is not a founder and not repeated and contains no genes in S

 foo <- descent("spot", ped, list(goofy = 1), debug = TRUE)
 foo
 checkitout(foo)

