
#include <R.h>
#include <R_ext/Utils.h>
#include <string.h>

// arguments pa and ma are one-origin indexing, as per R
// if individual i in zero-origin indexing is a founder, then
//     should have pa[i] == 0 && ma[i] == 0
// if individual i in zero-origin indexing is not a founder, then
//     should have pa[i] - 1 is the zero-origin index of father
//     and ma[i] - 1 is the zero-origin index of mother,
//     and in this case the requirement that individuals come before
//     their parents in the pedigree means pa[i] - 1 > i & ma[i] - 1 > i

// also argument args is one-origin indexing

void descent(int *nind_in, int *pa, int *ma, int *nargs_in, int *args,
    int *genes, double *result)
{
    int nind = nind_in[0];
    int nargs = nargs_in[0];
    const double half = 1.0 / 2.0;

    if (nind < 1)
        error("number of individuals in pedigree must be at least one");
    if (nargs < 0)
        error("number of individuals in arguments must be at least zero");

    for (int i = 0; i < nind; i++) {
        if ((pa[i] == 0) != (ma[i] == 0))
            error("must have both parents in pedigree or none");
        if ((pa[i] > nind) || (ma[i] > nind)) {
#ifdef BLATHER
            Rprintf("nind = %d\n", nind);
            Rprintf("i = %d (zero-origin indexing)\n", i);
            Rprintf("pa[i] = %d (one-origin indexing)\n", pa[i]);
            Rprintf("ma[i] = %d (one-origin indexing)\n", ma[i]);
#endif // BLATHER
            error("some parent index points outside of individual indices");
        }
        if ((pa[i] != 0) && ((pa[i] - 1 <= i) || (ma[i] - 1 <= i))) {
#ifdef BLATHER
            Rprintf("nind = %d\n", nind);
            Rprintf("i = %d (zero-origin indexing)\n", i);
            Rprintf("pa[i] = %d (one-origin indexing)\n", pa[i]);
            Rprintf("ma[i] = %d (one-origin indexing)\n", ma[i]);
#endif // BLATHER
            error("some individual comes after one of its parents in pedigree");
        }
    }

    for (int i = 1; i < nargs; i++)
        if ((args[i] <= 0 || args[i] > nind)) {
#ifdef OUCH
            Rprintf("nargs = %d\n", nargs);
            Rprintf("nind = %d\n", nind);
            for (int i = 0; i < nargs; i++)
                Rprintf("args[%d] = %d\n", i, args[i]);
#endif // OUCH
            error("argument individuals not range of indices of individuals");
        }

    for (int i = 0; i < nind; i++) {
        int genes_i = genes[i];
        if ((genes_i != 0) && (genes_i != 1) && (genes_i != 2))
            error("genes not 0, 1, or 2");
    }

    // end of checks

    // now refer to Theorem 1 in design document

    // by convention, nargs == 0 implies result = 1.0;
    if (nargs == 0) {
        result[0] = 1.0;
        return;
    }

    // sort args
    isort(args, nargs);

    int b1 = args[0];
    int r = 1;
    while ((r < nargs) && (args[r] == b1)) r = r + 1;

    // case (a), equation (1a)
    if ((pa[b1 - 1] != 0) && (genes[b1 - 1] == 0) && (r == 1)) {
        double foo, bar;
        int *myargs = (int *) R_alloc(nargs, sizeof(int));
        memcpy(myargs, args, nargs);
        myargs[0] = pa[b1 - 1];
        descent(&nind, pa, ma, &nargs, myargs, genes, &foo);
        myargs[0] = ma[b1 - 1];
        descent(&nind, pa, ma, &nargs, args, genes, &bar);
        result[0] = half * (foo + bar);
#ifdef WOOF
        Rprintf("nargs = %d, ", nargs);
        if (nargs == 0) {
            Rprintf("descent called with no argument individuals");
        } else {
            Rprintf("descent called with argument individuals %d", args[0]);
            for (int i = 1; i < nargs; i++)
                Rprintf(", %d", args[i]);
        }
        Rprintf(" (case a), result = %10.5g\n", result[0]);
#endif // WOOF
        return;
    }

    // case (b), equation (1b)
    if ((pa[b1 - 1] != 0) && (genes[b1 - 1] == 0) && (r > 1)) {
        double half_r_minus_one = half;
        for (int i = 2; i < r; i++) half_r_minus_one *= half;
        double foo, bar;
        int mynargs = nargs - r + 1;
        // add one extra at end for next use
        int *myargs = (int *) R_alloc(mynargs + 1, sizeof(int));
        myargs[0] = b1;
        memcpy(myargs + 1, args + r, nargs - r);
        descent(&nind, pa, ma, &mynargs, myargs, genes, &foo);
        mynargs += 1;
        myargs[0] = pa[b1 - 1];
        myargs[1] = ma[b1 - 1];
        memcpy(myargs + 2, args + r, nargs - r);
        descent(&nind, pa, ma, &mynargs, myargs, genes, &bar);
        result[0] = half_r_minus_one * foo + (1.0 - half_r_minus_one) * bar;
#ifdef WOOF
        if (nargs == 0) {
            Rprintf("descent called with no argument individuals");
        } else {
            Rprintf("descent called with argument individuals %d", args[0]);
            for (int i = 1; i < nargs; i++)
                Rprintf(", %d", args[i]);
        }
        Rprintf(" (case b), result = %10.5g\n", result[0]);
#endif // WOOF
        return;
    }

    // case (c), equation (1c)
    if ((pa[b1 - 1] == 0) && (genes[b1 - 1] == 0)) {
        result[0] = 0.0;
#ifdef WOOF
        if (nargs == 0) {
            Rprintf("descent called with no argument individuals");
        } else {
            Rprintf("descent called with argument individuals %d", args[0]);
            for (int i = 1; i < nargs; i++)
                Rprintf(", %d", args[i]);
        }
        Rprintf(" (case c), result = %10.5g\n", result[0]);
#endif // WOOF
        return;
    }

    // case (d), equation (1d)
    if ((genes[b1 - 1] == 2) && (r > 1)) {
        double foo;
        int mynargs = nargs - r;
        // again, by convention, nargs == 0 implies result = 1.0;
        if (mynargs == 0) {
            result[0] = 1.0;
#ifdef WOOF
            if (nargs == 0) {
                Rprintf("descent called with no argument individuals");
            } else {
                Rprintf("descent called with argument individuals %d", args[0]);
                for (int i = 1; i < nargs; i++)
                    Rprintf(", %d", args[i]);
            }
            Rprintf(" (case d), result = %10.5g\n", result[0]);
#endif // WOOF
            return;
        }
        int *myargs = (int *) R_alloc(mynargs, sizeof(int));
        memcpy(myargs, args + r, nargs - r);
        descent(nind_in, pa, ma, &mynargs, myargs, genes, &foo);
        result[0] = foo; 
#ifdef WOOF
        if (nargs == 0) {
            Rprintf("descent called with no argument individuals");
        } else {
            Rprintf("descent called with argument individuals %d", args[0]);
            for (int i = 1; i < nargs; i++)
                Rprintf(", %d", args[i]);
        }
        Rprintf(" (case d), result = %10.5g\n", result[0]);
#endif // WOOF
        return;
    }

    // case (e), equation (1e)
    if ((pa[b1 - 1] != 0) && (genes[b1 - 1] == 1)) {
        double half_r = half;
        for (int i = 2; i <= r; i++) half_r *= half;
        double foo, bar, baz;
        int mynargs = nargs - r;
        // add one extra at end for next use
        int *myargs = (int *) R_alloc(mynargs + 1, sizeof(int));
        memcpy(myargs, args + r, nargs - r);
        // again, by convention, nargs == 0 implies result = 1.0;
        if (mynargs == 0)
            foo = 1.0;
        else
            descent(nind_in, pa, ma, &mynargs, myargs, genes, &foo);
        mynargs += 1;
        myargs[0] = pa[b1 - 1];
        memcpy(myargs + 1, args + r, nargs - r);
        descent(nind_in, pa, ma, &mynargs, myargs, genes, &bar);
        myargs[0] = ma[b1 - 1];
        descent(nind_in, pa, ma, &mynargs, myargs, genes, &baz);
        result[0] = half_r * foo + half * (1.0 - half_r) * (bar + baz);
#ifdef WOOF
        if (nargs == 0) {
            Rprintf("descent called with no argument individuals");
        } else {
            Rprintf("descent called with argument individuals %d", args[0]);
            for (int i = 1; i < nargs; i++)
                Rprintf(", %d", args[i]);
        }
        Rprintf(" (case e), result = %10.5g\n", result[0]);
#endif // WOOF
        return;
    }

    // case (f), equation (1f)
    if ((pa[b1 - 1] == 0) && (genes[b1 - 1] == 1)) {
        double half_r = half;
        for (int i = 2; i <= r; i++) half_r *= half;
        double foo;
        int mynargs = nargs - r;
        // again, by convention, nargs == 0 implies result = 1.0;
        if (mynargs == 0) {
            foo = 1.0;
        } else {
            int *myargs = (int *) R_alloc(mynargs, sizeof(int));
            memcpy(myargs, args + r, nargs - r);
            descent(nind_in, pa, ma, &mynargs, myargs, genes, &foo);
        }
        result[0] = half_r * foo;
#ifdef WOOF
        if (nargs == 0) {
            Rprintf("descent called with no argument individuals");
        } else {
            Rprintf("descent called with argument individuals %d", args[0]);
            for (int i = 1; i < nargs; i++)
                Rprintf(", %d", args[i]);
        }
        Rprintf(" (case f), result = %10.5g\n", result[0]);
#endif // WOOF
        return;
    }

    // should never happen
    error("got to bottom of function without executing",
        " previous return statement\n");
}

