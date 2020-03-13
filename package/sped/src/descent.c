
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

static double get_value(SEXP foo);

static SEXP result_and_debug_info(double value, char *type, SEXP args,
    int ldebug, SEXP names, int ncalls, ...);

SEXP descent(SEXP pa, SEXP ma, SEXP args, SEXP genes, SEXP debug, SEXP names)
{
    const double half = 1.0 / 2.0;

    int nind = LENGTH(pa);
    int nargs = LENGTH(args);

    if (nind < 1)
        error("number of individuals in pedigree must be at least one");
    if (nargs < 0)
        error("number of individuals in arguments must be at least zero");

    if (! isLogical(debug))
        error("argument debug must be type logical");
    if (LENGTH(debug) != 1)
        error("argument debug must have length 1");
    int ldebug = LOGICAL(debug)[0];

    if (! isInteger(pa))
        error("argument pa must be type integer");
    if (! isInteger(ma))
        error("argument ma must be type integer");
    if (! isInteger(args))
        error("argument args must be type integer");
    if (! isInteger(genes))
        error("argument genes must be type integer");
    if (! (isString(names) || isInteger(names)))
        error("argument names must character or integer");

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
        if ((ipa[i] > nind) || (ima[i] > nind)) {
#ifdef BLEAT
            REprintf("i = %d\n", i);
            REprintf("nind = %d\n", nind);
            REprintf("ipa[%d] = %d\n", i, ipa[i]);
            REprintf("ima[%d] = %d\n", i, ima[i]);
#endif // BLEAT
            error("some parent index points outside of individual indices");
        }
        if ((ipa[i] != 0) && ((ipa[i] - 1 <= i) || (ima[i] - 1 <= i)))
            error("some individual comes after one of its parents in pedigree");
    }

    for (int i = 1; i < nargs; i++)
        if ((iargs[i] <= 0 || iargs[i] > nind)) {
#ifdef BLEAT
            REprintf("i = %d\n", i);
            REprintf("nind = %d\n", nind);
            REprintf("iargs[%d] = %d\n", i, iargs[i]);
#endif // BLEAT
            error("argument individuals not range of indices of individuals");
        }

    for (int i = 0; i < nind; i++) {
        int genes_i = igenes[i];
        if ((genes_i != 0) && (genes_i != 1) && (genes_i != 2))
            error("genes not 0, 1, or 2");
    }

    // now refer to Theorem 1 in design document

    // sort args
    R_isort(iargs, nargs);

    // case (a), equation (1a)
    // by convention, nargs == 0 implies result = 1.0;
    if (nargs == 0)
        return result_and_debug_info(1.0, "a", args, ldebug, names, 0);

    int b1 = iargs[0];
    int r = 1;
    while ((r < nargs) && (iargs[r] == b1)) r = r + 1;

    // case (b), equation (1b)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 0) && (r == 1)) {
        SEXP myargs, foo, bar;
        PROTECT(myargs = allocVector(INTSXP, nargs));
        memcpy(INTEGER(myargs), iargs, nargs * sizeof(int));
        INTEGER(myargs)[0] = ipa[b1 - 1];
        PROTECT(foo = descent(pa, ma, myargs, genes, debug, names));
        INTEGER(myargs)[0] = ima[b1 - 1];
        PROTECT(bar = descent(pa, ma, myargs, genes, debug, names));
        double dfoo = get_value(foo);
        double dbar = get_value(bar);
        UNPROTECT(3);
        double dvalue = half * (dfoo + dbar);
        return result_and_debug_info(dvalue, "b", args, ldebug,
            names, 2, foo, bar);
    }

    // case (c), equation (1c)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 0) && (r > 1)) {
        double half_r_minus_one = half;
        for (int i = 2; i < r; i++) half_r_minus_one *= half;
        int mynargs = nargs - r + 1;
        SEXP myargs, foo, bar;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        INTEGER(myargs)[0] = b1;
        memcpy(INTEGER(myargs) + 1, iargs + r, (nargs - r) * sizeof(int));
        PROTECT(foo = descent(pa, ma, myargs, genes, debug, names));
        mynargs += 1;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        INTEGER(myargs)[0] = ipa[b1 - 1];
        INTEGER(myargs)[1] = ima[b1 - 1];
        memcpy(INTEGER(myargs) + 2, iargs + r, (nargs - r) * sizeof(int));
        PROTECT(bar = descent(pa, ma, myargs, genes, debug, names));
        double dfoo = get_value(foo);
        double dbar = get_value(bar);
        UNPROTECT(4);
        double dvalue = half_r_minus_one * dfoo +
            (1.0 - half_r_minus_one) * dbar;
        return result_and_debug_info(dvalue, "c", args, ldebug,
            names, 2, foo, bar);
    }

    // case (d), equation (1d)
    if ((ipa[b1 - 1] == 0) && (igenes[b1 - 1] == 0))
        return result_and_debug_info(0.0, "d", args, ldebug, names, 0);

    // case (e), equation (1e)
    if ((igenes[b1 - 1] == 2) && (r > 1)) {
        int mynargs = nargs - r;
        SEXP myargs, foo;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        memcpy(INTEGER(myargs), iargs + r, (nargs - r) * sizeof(int));
        PROTECT(foo = descent(pa, ma, myargs, genes, debug, names));
        double dfoo = get_value(foo);
        UNPROTECT(2);
        return result_and_debug_info(dfoo, "e", args, ldebug, names, 1, foo);
    }

    // case (f), equation (1f)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 1)) {
        double half_r = half;
        for (int i = 2; i <= r; i++) half_r *= half;
        int mynargs = nargs - r;
        SEXP myargs, foo, bar, baz;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        memcpy(INTEGER(myargs), iargs + r, (nargs - r) * sizeof(int));
        PROTECT(foo = descent(pa, ma, myargs, genes, debug, names));
        mynargs += 1;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        INTEGER(myargs)[0] = ipa[b1 - 1];
        memcpy(INTEGER(myargs) + 1, iargs + r, (nargs - r) * sizeof(int));
        PROTECT(bar = descent(pa, ma, myargs, genes, debug, names));
        INTEGER(myargs)[0] = ima[b1 - 1];
        PROTECT(baz = descent(pa, ma, myargs, genes, debug, names));
        double dfoo = get_value(foo);
        double dbar = get_value(bar);
        double dbaz = get_value(baz);
        UNPROTECT(5);
        double dvalue = half_r * dfoo + half * (1.0 - half_r) * (dbar + dbaz);
        return result_and_debug_info(dvalue, "f", args, ldebug, names, 3,
            foo, bar, baz);
    }

    // case (g), equation (1g)
    if ((ipa[b1 - 1] == 0) && (igenes[b1 - 1] == 1)) {
        double half_r = half;
        for (int i = 2; i <= r; i++) half_r *= half;
        int mynargs = nargs - r;
        SEXP myargs, foo;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        memcpy(INTEGER(myargs), iargs + r, (nargs - r) * sizeof(int));
        PROTECT(foo = descent(pa, ma, myargs, genes, debug, names));
        double dfoo = get_value(foo);
        UNPROTECT(2);
        double dvalue = half_r * dfoo;
        return result_and_debug_info(dvalue, "g", args, ldebug, names, 1, foo);
    }

    // should never happen
    error("got to bottom of function without executing",
        " previous return statement\n");
}

static double get_value(SEXP foo)
{
    if (isVectorList(foo))
        foo = VECTOR_ELT(foo, 0);
    return REAL(foo)[0];
}

#include <stdarg.h>

static SEXP result_and_debug_info(double value, char *type, SEXP args,
    int ldebug, SEXP names, int ncalls, ...)
{
    if (! ldebug)
        return ScalarReal(value);

    // otherwise debug
    SEXP result, resultnames;
    PROTECT(result = allocVector(VECSXP, 3));
    PROTECT(resultnames = allocVector(STRSXP, 3));
    SET_STRING_ELT(resultnames, 0, mkChar("value"));
    SET_STRING_ELT(resultnames, 1, mkChar("individuals"));
    SET_STRING_ELT(resultnames, 2, mkChar("type"));
    SET_STRING_ELT(resultnames, 3, mkChar("calls"));
    namesgets(result, resultnames);
    SET_VECTOR_ELT(result, 0, ScalarReal(value));
    SET_VECTOR_ELT(result, 2, allocVector(STRSXP, 1));
    SEXP typestring = VECTOR_ELT(result, 2);
    SET_STRING_ELT(typestring, 0, mkChar(type));
    int is_integer_names = isInteger(names);
    int nargs = LENGTH(args);
    int *iargs = INTEGER(args);
    if (is_integer_names) {
        SET_VECTOR_ELT(result, 1, allocVector(INTSXP, nargs));
        SEXP mynames = VECTOR_ELT(result, 1);
        for (int i = 0; i < nargs; i++)
            INTEGER(mynames)[i] = INTEGER(names)[iargs[i] - 1];
    } else {
        SET_VECTOR_ELT(result, 1, allocVector(STRSXP, nargs));
        SEXP mynames = VECTOR_ELT(result, 1);
        for (int i = 0; i < nargs; i++)
            SET_STRING_ELT(mynames, i, STRING_ELT(names, iargs[i] - 1));
    }
    SET_VECTOR_ELT(result, 3, allocVector(VECSXP, ncalls));
    SEXP calls = VECTOR_ELT(result, 3);
    va_list argptr;
    va_start(argptr, ncalls);
    for (int i = 0; i < ncalls; i++) {
        SEXP foo = va_arg(argptr, SEXP);
        SET_VECTOR_ELT(calls, i, foo);
    }
    va_end(argptr);
    UNPROTECT(2);
    return result;
}
