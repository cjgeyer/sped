
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <string.h>

// arguments pa and ma are one-origin indexing, as always in R
// if individual i in zero-origin indexing is a founder, then
//     should have pa[i] == 0 && ma[i] == 0
// if individual i in zero-origin indexing is not a founder, then
//     should have pa[i] - 1 is the zero-origin index of father
//     and ma[i] - 1 is the zero-origin index of mother,
//     and in this case the requirement that individuals come before
//     their parents in the pedigree means pa[i] - 1 > i & ma[i] - 1 > i
//
// also argument args is one-origin indexing

SEXP descent(SEXP pa, SEXP ma, SEXP args, SEXP genes, SEXP debug)
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

    // case (a), equation (1a)
    if (nargs == 0) {
        // by convention, nargs == 0 implies result = 1.0;
        // don't bother with checks on other arguments or even on type of args
        // except we need to have checked debug
        if (! ldebug) {
            return ScalarReal(1.0);
        } else {
            SEXP result, resultnames;
            PROTECT(result = allocVector(VECSXP, 3));
            PROTECT(resultnames = allocVector(STRSXP, 3));
            SET_STRING_ELT(resultnames, 0, mkChar("value"));
            SET_STRING_ELT(resultnames, 1, mkChar("individuals"));
            SET_STRING_ELT(resultnames, 2, mkChar("type"));
            namesgets(result, resultnames);
            SET_VECTOR_ELT(result, 0, ScalarReal(1.0));
            SET_VECTOR_ELT(result, 1, allocVector(REALSXP, 0));
            SET_VECTOR_ELT(result, 2, allocVector(STRSXP, 1));
            SEXP typestring = VECTOR_ELT(result, 2);
            SET_STRING_ELT(typestring, 0, mkChar("a"));
            UNPROTECT(2);
            return result;
        }
    }

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
    int *ima = INTEGER(pa);
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

    // now refer to Theorem 1 in design document

    // sort args
    R_isort(iargs, nargs);

    int b1 = iargs[0];
    int r = 1;
    while ((r < nargs) && (iargs[r] == b1)) r = r + 1;

    // case (b), equation (1b)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 0) && (r == 1)) {
        SEXP myargs;
        PROTECT(myargs = allocVector(INTSXP, nargs));
        memcpy(INTEGER(myargs), iargs, nargs);
        INTEGER(myargs)[0] = ipa[b1 - 1];
        SEXP foo, bar;
        PROTECT(foo = descent(pa, ma, myargs, genes, debug));
        INTEGER(myargs)[0] = ima[b1 - 1];
        PROTECT(bar = descent(pa, ma, myargs, genes, debug));
        if (! ldebug) {
            double dfoo = REAL(foo)[0];
            double dbar = REAL(bar)[0];
            double dvalue = half * (dfoo + dbar);
            UNPROTECT(3);
            return ScalarReal(dvalue);
        } else {
            SEXP result, resultnames;
            PROTECT(result = allocVector(VECSXP, 4));
            PROTECT(resultnames = allocVector(STRSXP, 4));
            SET_STRING_ELT(resultnames, 0, mkChar("value"));
            SET_STRING_ELT(resultnames, 1, mkChar("individuals"));
            SET_STRING_ELT(resultnames, 2, mkChar("type"));
            SET_STRING_ELT(resultnames, 3, mkChar("calls"));
            namesgets(result, resultnames);
            double dfoo = REAL(VECTOR_ELT(foo, 0))[0];
            double dbar = REAL(VECTOR_ELT(bar, 0))[0];
            double dvalue = half * (dfoo + dbar);
            SET_VECTOR_ELT(result, 0, ScalarReal(dvalue));
            SET_VECTOR_ELT(result, 1, duplicate(args));
            SET_VECTOR_ELT(result, 2, allocVector(STRSXP, 1));
            SEXP typestring = VECTOR_ELT(result, 2);
            SET_STRING_ELT(typestring, 0, mkChar("b"));
            SET_VECTOR_ELT(result, 3, allocVector(VECSXP, 2));
            SEXP calls = VECTOR_ELT(result, 3);
            SET_VECTOR_ELT(calls, 0, foo);
            SET_VECTOR_ELT(calls, 1, bar);
            UNPROTECT(5);
            return result;
        }
    }

    // case (c), equation (1c)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 0) && (r > 1)) {
        double half_r_minus_one = half;
        for (int i = 2; i < r; i++) half_r_minus_one *= half;
        int mynargs = nargs - r + 1;
        SEXP myargs;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        INTEGER(myargs)[0] = b1;
        memcpy(INTEGER(myargs) + 1, INTEGER(args) + r, nargs - r);
        SEXP foo;
        PROTECT(foo = descent(pa, ma, myargs, genes, debug));
        mynargs += 1;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        INTEGER(myargs)[0] = ipa[b1 - 1];
        INTEGER(myargs)[1] = ima[b1 - 1];
        memcpy(INTEGER(myargs) + 2, INTEGER(args) + r, nargs - r);
        SEXP bar;
        PROTECT(bar = descent(pa, ma, myargs, genes, debug));
        if (! ldebug) {
            double dfoo = REAL(foo)[0];
            double dbar = REAL(bar)[0];
            double dvalue = half_r_minus_one * dfoo +
                (1.0 - half_r_minus_one) * dbar;
            UNPROTECT(3);
            return ScalarReal(dvalue);
        } else {
            SEXP result, resultnames;
            PROTECT(result = allocVector(VECSXP, 4));
            PROTECT(resultnames = allocVector(STRSXP, 4));
            SET_STRING_ELT(resultnames, 0, mkChar("value"));
            SET_STRING_ELT(resultnames, 1, mkChar("individuals"));
            SET_STRING_ELT(resultnames, 2, mkChar("type"));
            SET_STRING_ELT(resultnames, 3, mkChar("calls"));
            namesgets(result, resultnames);
            double dfoo = REAL(VECTOR_ELT(foo, 0))[0];
            double dbar = REAL(VECTOR_ELT(bar, 0))[0];
            double dvalue = half_r_minus_one * dfoo +
                (1.0 - half_r_minus_one) * dbar;
            SET_VECTOR_ELT(result, 0, ScalarReal(dvalue));
            SET_VECTOR_ELT(result, 1, duplicate(args));
            SET_VECTOR_ELT(result, 2, allocVector(STRSXP, 1));
            SEXP typestring = VECTOR_ELT(result, 2);
            SET_STRING_ELT(typestring, 0, mkChar("c"));
            SET_VECTOR_ELT(result, 3, allocVector(VECSXP, 2));
            SEXP calls = VECTOR_ELT(result, 3);
            SET_VECTOR_ELT(calls, 0, foo);
            SET_VECTOR_ELT(calls, 1, bar);
            UNPROTECT(5);
            return result;
        }
    }

    // case (d), equation (1d)
    if ((ipa[b1 - 1] == 0) && (igenes[b1 - 1] == 0)) {
        if (! ldebug) {
            return ScalarReal(0.0);
        } else {
            SEXP result, resultnames;
            PROTECT(result = allocVector(VECSXP, 3));
            PROTECT(resultnames = allocVector(STRSXP, 3));
            SET_STRING_ELT(resultnames, 0, mkChar("value"));
            SET_STRING_ELT(resultnames, 1, mkChar("individuals"));
            SET_STRING_ELT(resultnames, 2, mkChar("type"));
            namesgets(result, resultnames);
            SET_VECTOR_ELT(result, 0, ScalarReal(0.0));
            SET_VECTOR_ELT(result, 1, duplicate(args));
            SET_VECTOR_ELT(result, 2, allocVector(STRSXP, 1));
            SEXP typestring = VECTOR_ELT(result, 2);
            SET_STRING_ELT(typestring, 0, mkChar("d"));
            UNPROTECT(2);
            return result;
        }
    }

    // case (e), equation (1e)
    if ((igenes[b1 - 1] == 2) && (r > 1)) {
        int mynargs = nargs - r;
        SEXP myargs;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        memcpy(INTEGER(myargs), INTEGER(args) + r, nargs - r);
        SEXP foo;
        PROTECT(foo = descent(pa, ma, myargs, genes, debug));
        if (! ldebug) {
            double dvalue = REAL(foo)[0];
            UNPROTECT(2);
            return ScalarReal(dvalue);
        } else {
            SEXP result, resultnames;
            PROTECT(result = allocVector(VECSXP, 4));
            PROTECT(resultnames = allocVector(STRSXP, 4));
            SET_STRING_ELT(resultnames, 0, mkChar("value"));
            SET_STRING_ELT(resultnames, 1, mkChar("individuals"));
            SET_STRING_ELT(resultnames, 2, mkChar("type"));
            SET_STRING_ELT(resultnames, 3, mkChar("calls"));
            namesgets(result, resultnames);
            double dvalue = REAL(VECTOR_ELT(foo, 0))[0];
            SET_VECTOR_ELT(result, 0, ScalarReal(dvalue));
            SET_VECTOR_ELT(result, 1, duplicate(args));
            SET_VECTOR_ELT(result, 2, allocVector(STRSXP, 1));
            SEXP typestring = VECTOR_ELT(result, 2);
            SET_STRING_ELT(typestring, 0, mkChar("e"));
            SET_VECTOR_ELT(result, 3, allocVector(VECSXP, 1));
            SEXP calls = VECTOR_ELT(result, 3);
            SET_VECTOR_ELT(calls, 0, foo);
            UNPROTECT(4);
            return result;
        }

    }

    // case (f), equation (1f)
    if ((ipa[b1 - 1] != 0) && (igenes[b1 - 1] == 1)) {
        double half_r = half;
        for (int i = 2; i <= r; i++) half_r *= half;
        int mynargs = nargs - r;
        SEXP myargs;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        memcpy(INTEGER(myargs), INTEGER(args) + r, nargs - r);
        SEXP foo;
        PROTECT(foo = descent(pa, ma, myargs, genes, debug));
        mynargs += 1;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        INTEGER(myargs)[0] = ipa[b1 - 1];
        memcpy(INTEGER(myargs) + 1, INTEGER(args) + r, nargs - r);
        SEXP bar;
        PROTECT(bar = descent(pa, ma, myargs, genes, debug));
        INTEGER(myargs)[0] = ima[b1 - 1];
        SEXP baz;
        PROTECT(baz = descent(pa, ma, myargs, genes, debug));
        if (! ldebug) {
            double dfoo = REAL(foo)[0];
            double dbar = REAL(bar)[0];
            double dbaz = REAL(baz)[0];
            double dvalue = half_r * dfoo +
                half * (1.0 - half_r) * (dbar + dbaz);
            UNPROTECT(3);
            return ScalarReal(dvalue);
        } else {
            SEXP result, resultnames;
            PROTECT(result = allocVector(VECSXP, 4));
            PROTECT(resultnames = allocVector(STRSXP, 4));
            SET_STRING_ELT(resultnames, 0, mkChar("value"));
            SET_STRING_ELT(resultnames, 1, mkChar("individuals"));
            SET_STRING_ELT(resultnames, 2, mkChar("type"));
            SET_STRING_ELT(resultnames, 3, mkChar("calls"));
            namesgets(result, resultnames);
            double dfoo = REAL(VECTOR_ELT(foo, 0))[0];
            double dbar = REAL(VECTOR_ELT(bar, 0))[0];
            double dbaz = REAL(VECTOR_ELT(baz, 0))[0];
            double dvalue = half_r * dfoo +
                half * (1.0 - half_r) * (dbar + dbaz);
            SET_VECTOR_ELT(result, 0, ScalarReal(dvalue));
            SET_VECTOR_ELT(result, 1, duplicate(args));
            SET_VECTOR_ELT(result, 2, allocVector(STRSXP, 1));
            SEXP typestring = VECTOR_ELT(result, 2);
            SET_STRING_ELT(typestring, 0, mkChar("f"));
            SET_VECTOR_ELT(result, 3, allocVector(VECSXP, 3));
            SEXP calls = VECTOR_ELT(result, 3);
            SET_VECTOR_ELT(calls, 0, foo);
            SET_VECTOR_ELT(calls, 1, bar);
            SET_VECTOR_ELT(calls, 2, baz);
            UNPROTECT(5);
            return result;
        }
    }

    // case (g), equation (1g)
    if ((ipa[b1 - 1] == 0) && (igenes[b1 - 1] == 1)) {
        double half_r = half;
        for (int i = 2; i <= r; i++) half_r *= half;
        int mynargs = nargs - r;
        SEXP myargs;
        PROTECT(myargs = allocVector(INTSXP, mynargs));
        memcpy(INTEGER(myargs), INTEGER(args) + r, nargs - r);
        SEXP foo;
        PROTECT(foo = descent(pa, ma, myargs, genes, debug));
        if (! ldebug) {
            double dfoo = REAL(foo)[0];
            double dvalue = half_r * dfoo;
            UNPROTECT(2);
            return ScalarReal(dvalue);
        } else {
            SEXP result, resultnames;
            PROTECT(result = allocVector(VECSXP, 4));
            PROTECT(resultnames = allocVector(STRSXP, 4));
            SET_STRING_ELT(resultnames, 0, mkChar("value"));
            SET_STRING_ELT(resultnames, 1, mkChar("individuals"));
            SET_STRING_ELT(resultnames, 2, mkChar("type"));
            SET_STRING_ELT(resultnames, 3, mkChar("calls"));
            namesgets(result, resultnames);
            double dfoo = REAL(VECTOR_ELT(foo, 0))[0];
            double dvalue = half_r * dfoo;
            SET_VECTOR_ELT(result, 0, ScalarReal(dvalue));
            SET_VECTOR_ELT(result, 1, duplicate(args));
            SET_VECTOR_ELT(result, 2, allocVector(STRSXP, 1));
            SEXP typestring = VECTOR_ELT(result, 2);
            SET_STRING_ELT(typestring, 0, mkChar("g"));
            SET_VECTOR_ELT(result, 3, allocVector(VECSXP, 1));
            SEXP calls = VECTOR_ELT(result, 3);
            SET_VECTOR_ELT(calls, 0, foo);
            UNPROTECT(4);
            return result;
        }
    }

    // should never happen
    error("got to bottom of function without executing",
        " previous return statement\n");
}

