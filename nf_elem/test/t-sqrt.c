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

    /* test sqrt(a^2) */
    for (i = 0; i < 100 * antic_test_multiplier(); )
    {
        nf_t nf;
        nf_elem_t a, b, c, d;
        int is_square, num_facs;
        fmpz_poly_factor_t fac;
        fmpz_poly_t pol; /* do not clear */

        nf_init_randtest(nf, state, 30, 30);

        fmpz_poly_factor_init(fac);

        pol->coeffs = nf->pol->coeffs;
        pol->length = nf->pol->length;
        pol->alloc = nf->pol->alloc;

        fmpz_poly_factor(fac, pol);

        num_facs = fac->num*fac->exp[0];

        fmpz_poly_factor_clear(fac);

        if (nf->pol->length > 3 && (nf->flag & NF_MONIC) && num_facs == 1)
        {
           i++;

           nf_elem_init(a, nf);
           nf_elem_init(b, nf);
           nf_elem_init(c, nf);
           nf_elem_init(d, nf);

           nf_elem_randtest(a, state, 30, nf);

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
        }

        nf_clear(nf);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
