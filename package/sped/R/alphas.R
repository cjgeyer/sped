
alphas <- function(individuals, pedigree) {
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

    indall <- sort(unique(pedigree))
    founders <- setdiff(indall, pedigree[ , 1])

    result <- matrix(NA_real_, length(founders), length(individuals))
    geneset <- as.integer(1)
    for (i in seq(along = founders))
        for (j in seq(along = individuals)) {
            pedj <- pedigree[pedigree[ , 1] == individuals[j], , drop = FALSE]
            if (nrow(pedj) == 0) {
                # individual j is a founder
                result[i, j] <- 0
            } else if (nrow(pedj) == 1) {
                paj <- pedj[1, 2]
                maj <- pedj[1, 3]
                names(geneset) <- founders[i]
                result[i, j] <- 2 * descent(c(paj, maj), pedigree, geneset)
            }
        }
    rownames(result) <- founders
    colnames(result) <- individuals
    result
}

inbreeding <- function(individuals, pedigree)
    colSums(alphas(individuals, pedigree))
