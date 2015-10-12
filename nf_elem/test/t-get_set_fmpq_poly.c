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

    flint_printf("get/set fmpq_poly....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        fmpq_poly_t f;
        fmpq_poly_t g;
        nf_t nf;
        nf_elem_t a;

        fmpq_poly_init(pol);
        fmpq_poly_init(f);
        fmpq_poly_init(g);

        do {
           fmpq_poly_randtest_not_zero(pol, state, 40, 200);
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);

        fmpq_poly_randtest(f, state, fmpq_poly_degree(pol) - 1, 200);

        nf_elem_init(a, nf);
        nf_elem_set_fmpq_poly(a, f, nf);
        nf_elem_get_fmpq_poly(g, a, nf);
        
        result = fmpq_poly_equal(f, g);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "a");
            printf("\n");
            flint_printf("f = "); fmpq_poly_print_pretty(f, "x");
            printf("\n");
            flint_printf("g = "); fmpq_poly_print_pretty(f, "x");
            abort();
        }

        nf_elem_clear(a, nf);
        
        nf_clear(nf);

        fmpq_poly_clear(pol);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
