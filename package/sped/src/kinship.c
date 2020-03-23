
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

SEXP kinship(SEXP pa, SEXP ma)
{
    int nind = LENGTH(pa);

    if (nind < 1)
        error("number of individuals in pedigree must be at least one");

    if (! isInteger(pa))
        error("argument pa must be type integer");
    if (! isInteger(ma))
        error("argument ma must be type integer");

    if (LENGTH(ma) != nind)
        error("arguments pa and ma must have same length");

    int *ipa = INTEGER(pa);
    int *ima = INTEGER(ma);

    for (int i = 0; i < nind; i++) {
        if ((ipa[i] == 0) != (ima[i] == 0))
            error("must have both parents in pedigree or none");
        if ((ipa[i] > nind) || (ima[i] > nind))
            error("some parent index points outside of individual indices");
        if ((ipa[i] != 0) && ((ipa[i] - 1 <= i) || (ima[i] - 1 <= i)))
            error("some individual comes after one of its parents in pedigree");
    }

    SEXP result;
    PROTECT(result = allocMatrix(REALSXP, nind, nind));
    for (int i = 0; i < nind * nind; i++)
        REAL(result)[i] = R_NaReal;

#define RESULT(i, j) REAL(result)[(i) + nind * (j)]

    for (int j = nind - 1; j >= 0; j--) {
        if (ipa[j] == 0) {
            RESULT(j, j) = 0.5;
        } else {
            int idx_j_pa = ipa[j] - 1;
            int idx_j_ma = ima[j] - 1;
            RESULT(j, j) = (1.0 + RESULT(idx_j_pa, idx_j_ma)) / 2.0;
        }

        for (int i = j - 1; i >= 0; i--) {
            if (ipa[i] == 0) {
                RESULT(i, j) = 0.0;
                RESULT(j, i) = 0.0;
            } else {
                int idx_i_pa = ipa[i] - 1;
                int idx_i_ma = ima[i] - 1;
                double foo = (RESULT(idx_i_pa, j) + RESULT(idx_i_ma, j)) / 2.0;
                RESULT(i, j) = foo;
                RESULT(j, i) = foo;
            }
        }

        R_CheckUserInterrupt();
    }

    UNPROTECT(1);
    return result;
}

