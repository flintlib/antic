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

    flint_printf("norm....");
    fflush(stdout);

    flint_randinit(state);

    /* test norm(a*b) = norm(a)*norm(b) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;
        fmpq_t anorm, bnorm, cnorm, cnorm2;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 25, 200);
        } while (fmpq_poly_degree(pol) < 1);
        
        nf_init(nf, pol);
        
        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        
        fmpq_init(anorm);
        fmpq_init(bnorm);
        fmpq_init(cnorm);
        fmpq_init(cnorm2);

        nf_elem_randtest(a, state, 200, nf);
        nf_elem_randtest(b, state, 200, nf);
        
        nf_elem_mul(c, a, b, nf);
        nf_elem_norm(anorm, a, nf);
        nf_elem_norm(bnorm, b, nf);
        nf_elem_norm(cnorm, c, nf);
        fmpq_mul(cnorm2, anorm, bnorm);

        result = (fmpq_equal(cnorm, cnorm2));
        if (!result)
        {
           printf("FAIL:\n");
           printf("nf->pol = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
           printf("a = "); nf_elem_print(a, nf); printf("\n");
           printf("b = "); nf_elem_print(b, nf); printf("\n");
           printf("c = "); nf_elem_print(c, nf); printf("\n");
           printf("norm(a) = "); fmpq_print(anorm); printf("\n");
           printf("norm(b) = "); fmpq_print(bnorm); printf("\n");
           printf("norm(a*b) = "); fmpq_print(cnorm); printf("\n");
           printf("norm(a)*norm(b) = "); fmpq_print(cnorm2); printf("\n");
           abort();
        }

        fmpq_clear(anorm);
        fmpq_clear(bnorm);
        fmpq_clear(cnorm);
        fmpq_clear(cnorm2);

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
