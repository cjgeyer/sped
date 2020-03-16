descent <- function(individuals, pedigree, geneset, check.sex=FALSE) {
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
    stopifnot(is.logical(check.sex))
    stopifnot(length(check.sex) == 1)

    foo <- names(geneset)
    storage.mode(foo) <- storage.mode(pedigree)
    if (! all(foo %in% pedigree))
        stop("some geneset names not in pedigree")

    if (check.sex) {
        foo <- intersect(pedigree[,2], pedigree[,3])
        bar <- paste(foo, collapse = ", ")
        baz <- paste("individual(s)", bar, "is/are both father and mother")
        stop(baz)
    }

    from <- c(pedigree[,1], pedigree[,1])
    to <- c(pedigree[,2], pedigree[,3])
    foo <- try(tsort(from, to))
    if (inherits(foo, "try-error"))
        stop("some individual is its own ancestor")

    porder <- match(foo, pedigree[,1])
    ped.too <- pedigree[porder, ]
    pa <- match(ped.too[,2], foo)
    ma <- match(ped.too[,3], foo)
    pa[is.na(pa)] <- 0L
    ma[is.na(ma)] <- 0L

    iargs <- match(individuals, foo)

    genes <- rep(0L, length(foo))
    genes[match(names(geneset)[geneset == 1], foo)] <- 1L
    genes[match(names(geneset)[geneset == 2], foo)] <- 2L

    .Call(C_descent, pa = pa, ma = ma, args = iargs, genes = genes)
}
