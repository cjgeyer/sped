
library(sped, lib.loc = "../package/sped.Rcheck")

check.descent <- function(individuals, pedigree, geneset) {

    stopifnot(is.character(individuals) || is.numeric(individuals))
    if (is.numeric(individuals))
        storage.mode(individuals) <- "integer"
    if (is.integer(individuals) && any(individuals <= 0))
        stop("individuals, if integer-valued, must be positive-valued")
    stopifnot(is.matrix(pedigree))
    stopifnot(ncol(pedigree) == 3)
    stopifnot(is.character(pedigree) || is.numeric(pedigree))
    if (is.numeric(pedigree))
        storage.mode(pedigree) <- "integer"
    if (is.integer(pedigree) && any(pedigree <= 0))
        stop("pedigree, if integer-valued, must be positive-valued")
    stopifnot(typeof(individuals) == typeof(pedigree))
    stopifnot(individuals %in% pedigree)
    stopifnot(geneset %in% 0:2)
    stopifnot(length(names(geneset)) == length(geneset))

    foo <- names(geneset)
    storage.mode(foo) <- storage.mode(pedigree)
    if (! all(foo %in% pedigree))
        stop("some geneset names not in pedigree")

    from <- c(pedigree[,1], pedigree[,1])
    to <- c(pedigree[,2], pedigree[,3])
    foo <- try(pooh::tsort(from, to))
    if (inherits(foo, "try-error"))
        stop("some individual is its own ancestor")

    porder <- match(foo, pedigree[,1])
    ped.too <- pedigree[porder, ]

    # this function is to be called recursively to do the recursive
    # definition (see design doc) in R rather than C
    # note that pedigree, ped.too, geneset, etc. do not change
    # so do not need to be arguments to this function

    checker <- function(individuals) {

        # re-order individuals do offspring come before parents,
        # that is, in the foo ordering
        idx <- match(individuals, foo)
        if (any(is.na(idx)))
            stop("auxiliary function checker called with some ",
                "individual(s) not in pedigree")
        individuals <- foo[sort(idx)]

        # case (a) set of individuals is empty
        if (length(individuals) == 0) {
            return(list(value = 1.0, individuals = individuals,
                case = "a", calls = list()))
        }

        b1 <- individuals[1]
        ped.b1 <- ped.too[foo == b1, , drop = FALSE]
        if (nrow(ped.b1) != 1)
            stop("auxiliary function checker: individual(s) in pedigree",
                " multiple times")
        is.founder.b1 <- is.na(ped.b1[1, 2])
        geneset.b1 <- if (b1 %in% names(geneset)) {
            geneset[names(geneset) == b1] } else { 0 }

        # case (b) first individual is not a founder, occurs only once
        # in the vector of individuals, and contains no genes of S
        if ((! is.founder.b1) && (sum(b1 %in% individuals) == 1) &&
            (geneset.b1 == 0)) {

            pa.b1 <- ped.b1[1, 2]
            ma.b1 <- ped.b1[1, 3]

            rc1 <- checker(c(pa.b1, individuals[-1]))
            rc2 <- checker(c(ma.b1, individuals[-1]))

            value <- (rc1$value + rc2$value) / 2
            return(list(value = value, individuals = individuals,
                case = "b", calls = list(rc1, rc2)))
        }

        # case (c) first individual is not a founder, occurs repeatedly
        # in the vector of individuals, and contains no genes of S
        if ((! is.founder.b1) && (sum(b1 %in% individuals) > 1) &&
            (geneset.b1 == 0)) {

            r <- sum(b1 == individuals)
            rest <- individuals[b1 != individuals]
            pa.b1 <- ped.b1[1, 2]
            ma.b1 <- ped.b1[1, 3]

            rc1 <- checker(c(b1, rest))
            rc2 <- checker(c(pa.b1, ma.b1, rest))

            half.hat.r.minus.1 <- (1/2)^(r - 1)
            value <- half.hat.r.minus.1 * rc1$value +
                (1 - half.hat.r.minus.1) * rc2$value
            return(list(value = value, individuals = individuals,
                case = "c", calls = list(rc1, rc2)))
        }

        # case (d) first individual is a founder and contains no genes of S
        if (is.founder.b1 && (geneset.b1 == 0)) {
            return(list(value = 0, individuals = individuals,
                case = "d", calls = list()))
        }

        # case (e) first individual contains 2 genes of S
        if (geneset.b1 == 2) {
            rest <- individuals[b1 != individuals]
            rc1 <- checker(rest)
            value <- rc1$value
            return(list(value = value, individuals = individuals,
                case = "e", calls = list(rc1)))
        }

        # case (f) first individual is not a founder and contains 1 gene of S
        if ((! is.founder.b1) && (geneset.b1 == 1)) {

            r <- sum(b1 == individuals)
            rest <- individuals[b1 != individuals]
            pa.b1 <- ped.b1[1, 2]
            ma.b1 <- ped.b1[1, 3]

            rc1 <- checker(rest)
            rc2 <- checker(c(pa.b1, rest))
            rc3 <- checker(c(ma.b1, rest))

            half <- 1/2
            half.hat.r <- half^r
            value <- half.hat.r * rc1$value +
                half * (1 - half.hat.r) * (rc2$value + rc3$value)
            return(list(value = value, individuals = individuals,
                case = "f", calls = list(rc1, rc2, rc3)))
        }

        # case (g) first individual is a founder and contains 1 gene of S
        if (is.founder.b1 && (geneset.b1 == 1)) {

            r <- sum(b1 == individuals)
            rest <- individuals[b1 != individuals]

            rc1 <- checker(rest)

            half <- 1/2
            half.hat.r <- half^r
            value <- half.hat.r * rc1$value
            return(list(value = value, individuals = individuals,
                case = "g", calls = list(rc1)))
        }

    }

    checker(individuals)
}

data(alberta)

foo <- descent(c(1085, 1094, 1180, 1260, 1272), alberta, c("52"=1))
foo

bar <- check.descent(c(1085, 1094, 1180, 1260, 1272), alberta, c("52"=1))
bar$value

