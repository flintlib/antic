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

    flint_printf("is_square....");
    fflush(stdout);

    flint_randinit(state);

     /* test a^2 is a square */
    for (i = 0; i < 100 * antic_test_multiplier(); )
    {
        nf_t nf;
        nf_elem_t a, b;
        int num_facs;
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

           nf_elem_randtest_bounded(a, state, abits, nf);

           nf_elem_mul(b, a, a, nf);

           result = nf_elem_is_square(b, nf);
           
           if (!result)
           {
              printf("FAIL:\n");
              printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
              printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
              abort();
           }

           nf_elem_clear(a, nf);
           nf_elem_clear(b, nf);
        }

        nf_clear(nf);
    }
    
    /* test non-squares */
    for (i = 0; i < 100 * antic_test_multiplier(); )
    {
        nf_t nf;
        nf_elem_t a, b, z;
        int num_facs;
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
           nf_elem_init(z, nf);

           do {
               nf_elem_randtest_bounded(z, state, abits, nf);
           } while (nf_elem_is_square(z, nf));

           do {
               nf_elem_randtest_bounded(a, state, abits, nf);
           } while (nf_elem_is_zero(a, nf));

           nf_elem_mul(b, a, a, nf);
           nf_elem_mul(b, b, z, nf);

           result = !nf_elem_is_square(b, nf);
           
           if (!result)
           {
              printf("FAIL:\n");
              printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
              printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
              printf("z = "); nf_elem_print_pretty(z, nf, "x"); printf("\n");
              abort();
           }

           nf_elem_clear(a, nf);
           nf_elem_clear(b, nf);
           nf_elem_clear(z, nf);
        }

        nf_clear(nf);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
