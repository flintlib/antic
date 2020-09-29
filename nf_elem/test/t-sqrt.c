/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2020 William Hart

******************************************************************************/

#include <stdio.h>
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("sqrt....");
    fflush(stdout);

    flint_randinit(state);

    /* Regression test */
    {
        nf_t nf;
        nf_elem_t a, b, c, d;
        int is_square;
        fmpq_poly_t f;

        /* f = x^3 - 1953*x^2 - x + 1 */
        fmpq_poly_init(f);
        fmpq_poly_set_coeff_si(f, 0, 1);
        fmpq_poly_set_coeff_si(f, 1, -1);
        fmpq_poly_set_coeff_si(f, 2, -1953);
        fmpq_poly_set_coeff_si(f, 3, 1);

        nf_init(nf, f);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        nf_elem_init(d, nf);

        /* a = -26977/6*x^2+40549/6 */
        fmpq_poly_set_coeff_si(NF_ELEM(a), 0, 40549);
        fmpq_poly_set_coeff_si(NF_ELEM(a), 2, -26977);
        fmpz_set_ui(NF_ELEM_DENREF(a), 6);

        nf_elem_mul(b, a, a, nf);

        is_square = nf_elem_sqrt(c, b, nf);

        nf_elem_mul(d, c, c, nf);

        result = is_square && nf_elem_equal(d, b, nf);
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("d = "); nf_elem_print_pretty(d, nf, "x"); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        nf_elem_clear(d, nf);

        nf_clear(nf);

        fmpq_poly_clear(f);
    }

    /* Regression test */
    {
        nf_t nf;
        nf_elem_t a, b, c, d;
        int is_square;
        fmpq_poly_t f;

        /* f = x^3 + 44*x^2 - 3*x + 18 */
        fmpq_poly_init(f);
        fmpq_poly_set_coeff_si(f, 0, 18);
        fmpq_poly_set_coeff_si(f, 1, -3);
        fmpq_poly_set_coeff_si(f, 2, 44);
        fmpq_poly_set_coeff_si(f, 3, 1);

        nf_init(nf, f);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        nf_elem_init(d, nf);

        /* a = -1/80916*x^2 + 1/80916*x + 2/80916 */
        fmpq_poly_set_coeff_si(NF_ELEM(a), 0, 2);
        fmpq_poly_set_coeff_si(NF_ELEM(a), 1, 1);
        fmpq_poly_set_coeff_si(NF_ELEM(a), 2, -1);
        fmpz_set_ui(NF_ELEM_DENREF(a), 80916);

        nf_elem_mul(b, a, a, nf);

        is_square = nf_elem_sqrt(c, b, nf);

        nf_elem_mul(d, c, c, nf);

        result = is_square && nf_elem_equal(d, b, nf);
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("d = "); nf_elem_print_pretty(d, nf, "x"); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        nf_elem_clear(d, nf);

        nf_clear(nf);

        fmpq_poly_clear(f);
    }

    /* Regression test */
    {
        nf_t nf;
        nf_elem_t a, b, c, d;
        int is_square;
        fmpq_poly_t f;

        /* f = x^3 + 55*x^2 - 10*x - 1 */
        fmpq_poly_init(f);
        fmpq_poly_set_coeff_si(f, 0, -1);
        fmpq_poly_set_coeff_si(f, 1, -10);
        fmpq_poly_set_coeff_si(f, 2, 55);
        fmpq_poly_set_coeff_si(f, 3, 1);

        nf_init(nf, f);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        nf_elem_init(d, nf);

        /* a = 862/16*x^2 - 149/16*x - 13/16 */
        fmpq_poly_set_coeff_si(NF_ELEM(a), 0, -13);
        fmpq_poly_set_coeff_si(NF_ELEM(a), 1, -149);
        fmpq_poly_set_coeff_si(NF_ELEM(a), 2, 862);
        fmpz_set_ui(NF_ELEM_DENREF(a), 16);

        nf_elem_mul(b, a, a, nf);

        is_square = nf_elem_sqrt(c, b, nf);

        nf_elem_mul(d, c, c, nf);

        result = is_square && nf_elem_equal(d, b, nf);
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("d = "); nf_elem_print_pretty(d, nf, "x"); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        nf_elem_clear(d, nf);

        nf_clear(nf);

        fmpq_poly_clear(f);
    }

    /* Regression test */
    {
        nf_t nf;
        nf_elem_t a, b, c, d;
        int is_square;
        fmpq_poly_t f;

        /* f = x^8 - 8 */
        fmpq_poly_init(f);
        fmpq_poly_set_coeff_si(f, 0, -8);
        fmpq_poly_set_coeff_si(f, 8, 1);

        nf_init(nf, f);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        nf_elem_init(d, nf);

        /* a = 4050*x^6 */
        fmpq_poly_set_coeff_si(NF_ELEM(a), 6, 4050);
        fmpz_set_ui(NF_ELEM_DENREF(a), 1);

        nf_elem_mul(b, a, a, nf);

        is_square = nf_elem_sqrt(c, b, nf);

        nf_elem_mul(d, c, c, nf);

        result = is_square && nf_elem_equal(d, b, nf);
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("d = "); nf_elem_print_pretty(d, nf, "x"); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        nf_elem_clear(d, nf);

        nf_clear(nf);

        fmpq_poly_clear(f);
    }

    /* test sqrt(a^2) monic defining poly */
    for (i = 0; i < 100 * antic_test_multiplier(); )
    {
        nf_t nf;
        nf_elem_t a, b, c, d;
        int is_square, num_facs;
        slong flen, fbits, abits;
        fmpz_poly_factor_t fac;
        fmpz_poly_t pol; /* do not clear */

        flen = n_randint(state, 30) + 2;
        fbits = n_randint(state, 30) + 1;
        abits = n_randint(state, 30) + 1;
        
        nf_init_randtest(nf, state, flen, fbits);

        fmpz_poly_factor_init(fac);

        pol->coeffs = nf->pol->coeffs;
        pol->length = nf->pol->length;
        pol->alloc = nf->pol->alloc;

        fmpz_poly_factor(fac, pol);

        num_facs = fac->num*fac->exp[0];

        fmpz_poly_factor_clear(fac);

        if (nf->flag & NF_MONIC && num_facs == 1)
        {
           i++;

           nf_elem_init(a, nf);
           nf_elem_init(b, nf);
           nf_elem_init(c, nf);
           nf_elem_init(d, nf);
           
           nf_elem_randtest_bounded(a, state, abits, nf);

           nf_elem_mul(b, a, a, nf);

           is_square = nf_elem_sqrt(c, b, nf);
           
           nf_elem_mul(d, c, c, nf);

           result = is_square && nf_elem_equal(d, b, nf);
           if (!result)
           {
              printf("FAIL:\n");
              printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
              printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
              printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
              printf("d = "); nf_elem_print_pretty(d, nf, "x"); printf("\n");
              abort();
           }

           nf_elem_clear(a, nf);
           nf_elem_clear(b, nf);
           nf_elem_clear(c, nf);
           nf_elem_clear(d, nf);
        }

        nf_clear(nf);
    }

     /* test sqrt(a^2) non-monic defining poly */
    for (i = 0; i < 100 * antic_test_multiplier(); )
    {
        nf_t nf;
        nf_elem_t a, b, c, d;
        int is_square, num_facs;
        slong flen, fbits, abits;
        fmpz_poly_factor_t fac;
        fmpz_poly_t pol; /* do not clear */

        flen = n_randint(state, 10) + 2;
        fbits = n_randint(state, 10) + 1;
        abits = n_randint(state, 10) + 1;
        
        nf_init_randtest(nf, state, flen, fbits);

        fmpz_poly_factor_init(fac);

        pol->coeffs = nf->pol->coeffs;
        pol->length = nf->pol->length;
        pol->alloc = nf->pol->alloc;

        fmpz_poly_factor(fac, pol);

        num_facs = fac->num*fac->exp[0];

        fmpz_poly_factor_clear(fac);

        if (num_facs == 1)
        {
           i++;

           nf_elem_init(a, nf);
           nf_elem_init(b, nf);
           nf_elem_init(c, nf);
           nf_elem_init(d, nf);
           
           nf_elem_randtest_bounded(a, state, abits, nf);

           nf_elem_mul(b, a, a, nf);

           is_square = nf_elem_sqrt(c, b, nf);
           
           nf_elem_mul(d, c, c, nf);

           result = is_square && nf_elem_equal(d, b, nf);
           if (!result)
           {
              printf("FAIL:\n");
              printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
              printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
              printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
              printf("d = "); nf_elem_print_pretty(d, nf, "x"); printf("\n");
              abort();
           }

           nf_elem_clear(a, nf);
           nf_elem_clear(b, nf);
           nf_elem_clear(c, nf);
           nf_elem_clear(d, nf);
        }

        nf_clear(nf);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
