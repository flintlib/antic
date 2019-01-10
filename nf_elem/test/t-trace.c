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
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("trace....");
    fflush(stdout);

    flint_randinit(state);

    /* test trace(a + b) = trace(a) + trace(b) */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;
        fmpq_t atrace, btrace, ctrace, ctrace2;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 25, 200);
        } while (fmpq_poly_degree(pol) < 1);
        
        nf_init(nf, pol);
        
        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        
        fmpq_init(atrace);
        fmpq_init(btrace);
        fmpq_init(ctrace);
        fmpq_init(ctrace2);

        nf_elem_randtest(a, state, 200, nf);
        nf_elem_randtest(b, state, 200, nf);
        
        nf_elem_add(c, a, b, nf);
        nf_elem_trace(atrace, a, nf);
        nf_elem_trace(btrace, b, nf);
        nf_elem_trace(ctrace, c, nf);
        fmpq_add(ctrace2, atrace, btrace);

        result = (fmpq_equal(ctrace, ctrace2));
        if (!result)
        {
           printf("FAIL:\n");
           printf("nf->pol = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("trace(a) = "); fmpq_print(atrace); printf("\n");
           printf("trace(b) = "); fmpq_print(btrace); printf("\n");
           printf("trace(a + b) = "); fmpq_print(ctrace); printf("\n");
           printf("trace(a) + trace(b) = "); fmpq_print(ctrace2); printf("\n");
           abort();
        }

        fmpq_clear(atrace);
        fmpq_clear(btrace);
        fmpq_clear(ctrace);
        fmpq_clear(ctrace2);

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
