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

    Copyright (C) 2018 Tommy Hofmann

******************************************************************************/

#include <stdio.h>
#include "flint/fmpq_mat.h"
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    flint_printf("rep_mat_fmpz_mat_den....");
    fflush(stdout);

    flint_randinit(state);

    /* test mul_gen(b) = a * b, where a is the generator */
    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, p1, p2, t;
        slong d;
        slong j, k;
        fmpz_mat_t R;
        fmpz_t den;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 20, 100);
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);

        d = fmpq_poly_degree(pol);

        fmpz_mat_init(R, d, d);

        fmpz_init(den);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);
        nf_elem_init(t, nf);

        nf_elem_randtest(b, state, 100, nf);

        nf_elem_rep_mat_fmpz_mat_den(R, den, b, nf);

        /* fmpz_mat_print_pretty(R); */

        for (j = 0; j < d; j++)
        {
            nf_elem_gen(a, nf);
            nf_elem_pow(a, a, j, nf);
            nf_elem_mul(p1, b, a, nf);

            nf_elem_zero(p2, nf);

            for (k = 0; k < d; k++)
            {
                nf_elem_gen(t, nf);
                nf_elem_pow(t, t, k, nf);
                nf_elem_scalar_mul_fmpz(t, t, fmpz_mat_entry(R, j, k), nf);
                nf_elem_add(p2, p2, t, nf);
            }

            nf_elem_scalar_div_fmpz(p2, p2, den, nf);

            if (!nf_elem_equal(p1, p2, nf))
            {
                printf("FAIL:\n");
                printf("R = "); fmpz_mat_print_pretty(R); printf("\n");
                printf("d = "); fmpz_print(den); printf("\n");
                printf("K = "); nf_print(nf); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
                printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
                abort();
            }
        }


        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);
        nf_elem_clear(t, nf);
        fmpz_mat_clear(R);
        fmpz_clear(den);

        nf_clear(nf);

        fmpq_poly_clear(pol);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
