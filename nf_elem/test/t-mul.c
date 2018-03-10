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
    int i, result;
    flint_rand_t state;

    flint_printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    /* test a*(b + c) = a*b + a*c */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c, s, p, p1, p2;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 40, 200);
        } while (fmpq_poly_degree(pol) < 1);
        
        nf_init(nf, pol);
        
        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        nf_elem_init(s, nf);
        nf_elem_init(p, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);

        nf_elem_randtest(a, state, 200, nf);
        nf_elem_randtest(b, state, 200, nf);
        nf_elem_randtest(c, state, 200, nf);
        
        nf_elem_add(s, b, c, nf);
        nf_elem_mul(p1, a, b, nf);
        nf_elem_mul(p2, a, c, nf);

        nf_elem_mul(p, a, s, nf);
        nf_elem_add(s, p1, p2, nf);

        result = (nf_elem_equal(p, s, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("s = "); nf_elem_print_pretty(s, nf, "x"); printf("\n");
           printf("p = "); nf_elem_print_pretty(p, nf, "x"); printf("\n");
           printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
           printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        nf_elem_clear(s, nf);
        nf_elem_clear(p, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(pol);
    }
    
    /* test aliasing a and b */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 40, 200);
        } while (fmpq_poly_degree(pol) < 1);
        
        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(b, state, 200, nf);
        nf_elem_randtest(c, state, 200, nf);
        
        nf_elem_mul(a, b, c, nf);
        nf_elem_mul(b, b, c, nf);
        
        result = (nf_elem_equal(a, b, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(pol);
    }

    /* test aliasing a and c */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 40, 200);
        } while (fmpq_poly_degree(pol) < 1);
        
        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(b, state, 200, nf);
        nf_elem_randtest(c, state, 200, nf);
        
        nf_elem_mul(a, b, c, nf);
        nf_elem_mul(c, b, c, nf);
        
        result = (nf_elem_equal(a, c, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(pol);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
