
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <string.h>
#include "sped.h"

// arguments pa and ma are one-origin indexing, as always in R
//     zero value indicates absence of parent (so individual is founder)
// if individual i in zero-origin indexing is a founder, then
//     should have pa[i] == 0 && ma[i] == 0
// if individual i in zero-origin indexing is not a founder, then
//     should have pa[i] - 1 is the zero-origin index of father
//     and ma[i] - 1 is the zero-origin index of mother,
//     and in this case the requirement that individuals come before
//     their parents in the pedigree means pa[i] - 1 > i & ma[i] - 1 > i
//
// also argument args is one-origin indexing

static double my_descent(int nind, int *ipa, int *ima, int *igenes,
    int nargs, int *iargs);

SEXP descent(SEXP pa, SEXP ma, SEXP args, SEXP genes)
{
    int nind = LENGTH(pa);
    int nargs = LENGTH(args);

    if (nind < 1)
        error("number of individuals in pedigree must be at least one");

    if (! isInteger(pa))
        error("argument pa must be type integer");
    if (! isInteger(ma))
        error("argument ma must be type integer");
    if (! isInteger(args))
        error("argument args must be type integer");
    if (! isInteger(genes))
        error("argument genes must be type integer");

    if (LENGTH(ma) != nind)
        error("arguments pa and ma must have same length");
    if (LENGTH(genes) != nind)
        error("arguments pa and genes must have same length");

    int *ipa = INTEGER(pa);
    int *ima = INTEGER(ma);
    int *iargs = INTEGER(args);
    int *igenes = INTEGER(genes);

    for (int i = 0; i < nind; i++) {
        if ((ipa[i] == 0) != (ima[i] == 0))
            error("must have both parents in pedigree or none");
        if ((ipa[i] > nind) || (ima[i] > nind))
            error("some parent index points outside of individual indices");
        if ((ipa[i] != 0) && ((ipa[i] - 1 <= i) || (ima[i] - 1 <= i)))
            error("some individual comes after one of its parents in pedigree");
    }

    for (int i = 1; i < nargs; i++)
        if ((iargs[i] <= 0 || iargs[i] > nind))
            error("argument individuals not range of indices of individuals");

    for (int i = 0; i < nind; i++) {
        int genes_i = igenes[i];
        if ((genes_i != 0) && (genes_i != 1) && (genes_i != 2))
            error("genes not 0, 1, or 2");
    }

    return ScalarReal(my_descent(nind, ipa, ima, igenes, nargs, iargs));
}

static double my_descent(int nind, int *ipa, int *ima, int *igenes,
    int nargs, int *iargs) {

    const double half = 1.0 / 2.0;

    R_CheckUserInterrupt();

    // refer to Theorem 1 in design document (in inst/DesignDoc)

    // sort args
    R_isort(iargs, nargs);

    // case (a), equation (1a)
    // by convention, nargs == 0 implies result = 1.0;
    if (nargs == 0)
        return 1.0;

    int b1 = iargs[0];
    int r = 1;
    while ((r < nargs) && (iargs[r] == b1)) r = r + 1;

    // case (b), equation (1b)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 0) && (r == 1)) {
        int *myargs = (int *) R_alloc(nargs, sizeof(int));
        memcpy(myargs, iargs, nargs * sizeof(int));
        myargs[0] = ipa[b1 - 1];
        double dfoo = my_descent(nind, ipa, ima, igenes, nargs, myargs);
        memcpy(myargs, iargs, nargs * sizeof(int));
        myargs[0] = ima[b1 - 1];
        double dbar = my_descent(nind, ipa, ima, igenes, nargs, myargs);
        return half * (dfoo + dbar);
    }

    // case (c), equation (1c)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 0) && (r > 1)) {
        double half_r_minus_one = half;
        for (int i = 2; i < r; i++) half_r_minus_one *= half;
        int mynargs = nargs - r + 1;
        // allocate one extra item for second call to my_descent
        int *myargs = (int *) R_alloc(mynargs + 1, sizeof(int));
        myargs[0] = b1;
        memcpy(myargs + 1, iargs + r, (nargs - r) * sizeof(int));
        double dfoo = my_descent(nind, ipa, ima, igenes, mynargs, myargs);
        mynargs += 1;
        myargs[0] = ipa[b1 - 1];
        myargs[1] = ima[b1 - 1];
        memcpy(myargs + 2, iargs + r, (nargs - r) * sizeof(int));
        double dbar = my_descent(nind, ipa, ima, igenes, mynargs, myargs);
        return half_r_minus_one * dfoo + (1.0 - half_r_minus_one) * dbar;
    }

    // case (d), equation (1d)
    if ((ipa[b1 - 1] == 0) && (igenes[b1 - 1] == 0))
        return 0.0;

    // case (e), equation (1e)
    if (igenes[b1 - 1] == 2) {
        int mynargs = nargs - r;
        int *myargs = (int *) R_alloc(mynargs, sizeof(int));
        memcpy(myargs, iargs + r, (nargs - r) * sizeof(int));
        return my_descent(nind, ipa, ima, igenes, mynargs, myargs);
    }

    // case (f), equation (1f)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 1)) {
        double half_r = half;
        for (int i = 2; i <= r; i++) half_r *= half;
        int mynargs = nargs - r;
        // allocate one extra item for second and third call to my_descent
        int *myargs = (int *) R_alloc(mynargs + 1, sizeof(int));
        memcpy(myargs, iargs + r, (nargs - r) * sizeof(int));
        double dfoo = my_descent(nind, ipa, ima, igenes, mynargs, myargs);
        mynargs += 1;
        myargs[0] = ipa[b1 - 1];
        memcpy(myargs + 1, iargs + r, (nargs - r) * sizeof(int));
        double dbar = my_descent(nind, ipa, ima, igenes, mynargs, myargs);
        myargs[0] = ima[b1 - 1];
        memcpy(myargs + 1, iargs + r, (nargs - r) * sizeof(int));
        double dbaz = my_descent(nind, ipa, ima, igenes, mynargs, myargs);
        return half_r * dfoo + half * (1.0 - half_r) * (dbar + dbaz);
    }

    // case (g), equation (1g)
    if ((ipa[b1 - 1] == 0) && (igenes[b1 - 1] == 1)) {
        double half_r = half;
        for (int i = 2; i <= r; i++) half_r *= half;
        int mynargs = nargs - r;
        int *myargs = (int *) R_alloc(mynargs, sizeof(int));
        memcpy(myargs, iargs + r, (nargs - r) * sizeof(int));
        double dfoo = my_descent(nind, ipa, ima, igenes, mynargs, myargs);
        return half_r * dfoo;
    }

    // should never happen
    error("got to bottom of function my_descent without executing"
        " previous return statement\n");
}

