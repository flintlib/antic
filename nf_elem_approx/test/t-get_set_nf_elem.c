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
#include "nf_elem_approx.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("get/set_nf_elem....");
    fflush(stdout);

    flint_randinit(state);

    /* exact->approx->exact */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b;
        nf_elem_approx_t c;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 4, 20);
        } while (fmpq_poly_degree(pol) != 3
           || !_fmpz_poly_is_squarefree(fmpq_poly_numref(pol), 4));

        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_approx_init(c, nf);

        nf_elem_randtest(b, state, 20, nf);
        
        nf_elem_approx_set_nf_elem(c, b, nf);
        c->exact = 0; /* force it to convert back to exact element */
        nf_elem_approx_get_nf_elem(a, c, nf);

        result = (nf_elem_equal(a, b, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print(a, nf); printf("\n");
           printf("b = "); nf_elem_print(b, nf); printf("\n");
           abort();
        }

        nf_elem_approx_clear(c, nf);
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
