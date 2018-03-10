/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 William Hart

******************************************************************************/

#include <stdio.h>
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, j, result;
    flint_rand_t state;

    flint_printf("pow....");
    fflush(stdout);

    flint_randinit(state);

    /* test pow(a, e) = e*e*...*e */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, p1, p2;
        slong exp;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 40, 20);
        } while (fmpq_poly_degree(pol) < 1);
        
        nf_init(nf, pol);
        
        nf_elem_init(a, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);

        nf_elem_randtest(a, state, 20, nf);
        
        exp = n_randint(state, 10);

        nf_elem_pow(p1, a, exp, nf);
        nf_elem_one(p2, nf);

        for (j = 0; j < exp; j++)
           nf_elem_mul(p2, p2, a, nf);
        
        result = (nf_elem_equal(p1, p2, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
           printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
           flint_printf("exp = %w\n", exp);
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(pol);
    }
    
    /* test aliasing a and res */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, p1, p2;
        slong exp;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 40, 20);
        } while (fmpq_poly_degree(pol) < 1);
        
        nf_init(nf, pol);
        
        nf_elem_init(a, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);

        nf_elem_randtest(a, state, 20, nf);
        
        exp = n_randint(state, 10);

        nf_elem_pow(p1, a, exp, nf);
        nf_elem_set(p2, a, nf);
        nf_elem_pow(p2, p2, exp, nf);
        
        result = (nf_elem_equal(p1, p2, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
           printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
           flint_printf("exp = %w\n", exp);
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(pol);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
