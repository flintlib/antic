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

    Copyright (C) 2014 William Hart

******************************************************************************/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "nf.h"
#include "nf_elem.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("div....");
    fflush(stdout);

    flint_randinit(state);

    /* test a*^-1 = 1 */
    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        fmpq_poly_t g, pol;
        nf_t nf;
        nf_elem_t a, b, c;

        fmpq_poly_init(g);
        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 25, 100);
        } while (fmpq_poly_degree(pol) < 1);
   
        nf_init(nf, pol);
        
        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        
        do {
           nf_elem_randtest_not_zero(a, state, 100, nf);
           fmpq_poly_gcd(g, NF_ELEM(a), pol);
        } while (!fmpq_poly_is_one(g));
        nf_elem_randtest(b, state, 100, nf);
           
        nf_elem_div(c, b, a, nf);
        nf_elem_mul(c, c, a, nf);
        
        result = (nf_elem_equal(b, c, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print(a, nf); printf("\n");
           printf("b = "); nf_elem_print(b, nf); printf("\n");
           printf("c = "); nf_elem_print(c, nf); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(g);
        fmpq_poly_clear(pol);
    }
    
    /* test aliasing a and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;

        fmpq_poly_init(pol);
        fmpq_poly_randtest_not_zero(pol, state, 25, 100);
        
        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        
        nf_elem_randtest(b, state, 100, nf);
        nf_elem_randtest(c, state, 100, nf);
        
        nf_elem_div(a, b, c, nf);
        nf_elem_div(b, b, c, nf);
        
        result = (nf_elem_equal(a, b, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print(a, nf); printf("\n");
           printf("b = "); nf_elem_print(b, nf); printf("\n");
           printf("c = "); nf_elem_print(c, nf); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(pol);
    }

    /* test aliasing a and c */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;

        fmpq_poly_init(pol);
        fmpq_poly_randtest_not_zero(pol, state, 25, 100);
        
        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        
        nf_elem_randtest(b, state, 100, nf);
        nf_elem_randtest(c, state, 100, nf);
        
        nf_elem_div(a, b, c, nf);
        nf_elem_div(c, b, c, nf);
        
        result = (nf_elem_equal(a, c, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print(a, nf); printf("\n");
           printf("d = "); nf_elem_print(b, nf); printf("\n");
           printf("c = "); nf_elem_print(c, nf); printf("\n");
           abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
         
        nf_clear(nf);

        fmpq_poly_clear(pol);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
