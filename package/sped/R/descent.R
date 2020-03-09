descent <- function(individuals, pedigree, geneset, check.sex=FALSE,
    debug=FALSE) {
    stopifnot(is.atomic(individuals))
    stopifnot(is.matrix(pedigree))
    stopifnot(ncol(pedigree) == 3)
    stopifnot(typeof(individuals) == typeof(pedigree))
    stopifnot(individuals %in% pedigree)
    stopifnot(geneset %in% 0:2)
    stopifnot(length(names(geneset)) == length(geneset))
    stopifnot(is.logical(check.sex))
    stopifnot(length(check.sex) == 1)
    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)

    foo <- names(geneset)
    storage.mode(foo) <- storage.mode(pedigree)
    if (! all(foo %in% pedigree))
        stop("some geneset names not in pedigree")

    if (check.sex) {
        foo <- intersect(pedigree[,2], pedigree[,3])
        bar <- paste(foo, collapse = ", ")
        baz <- paste("individuals", bar, "are both father and mother")
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

    .Call("descent", pa = pa, ma = ma, args = iargs, genes = genes,
        debug = debug, PACKAGE = "sped")
}
