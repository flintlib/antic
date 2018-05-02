/*=============================================================================

    This file is part of ANTIC.

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
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("get_nmod_poly....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        slong j;
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a;
        nmod_poly_t reduced_elem;
        ulong mod;
        fmpz_t coeff;

        mod = n_randtest_not_zero(state);

        fmpz_init(coeff);
        nmod_poly_init(reduced_elem, mod);

        fmpq_poly_init(pol);
        do {
            fmpq_poly_randtest_not_zero(pol, state, 40, 200);
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);

        nf_elem_init(a, nf);

        nf_elem_randtest(a, state, 200, nf);

        nf_elem_get_nmod_poly_den(reduced_elem, a, nf, 0);

        for (j = 0; j < fmpq_poly_degree(pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            result = (nmod_poly_get_coeff_ui(reduced_elem, j) == fmpz_fdiv_ui(coeff, mod));
            if (!result)
            {
                printf("FAIL: Reducing without denominator\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); flint_printf("%u\n", mod);
                printf("a mod n = "); nmod_poly_print_pretty(reduced_elem, "x"); printf("\n");
                abort();
            }
        }

        nf_elem_clear(a, nf);
        nmod_poly_clear(reduced_elem);
        fmpz_clear(coeff);

        nf_clear(nf);

        fmpq_poly_clear(pol);
    }

    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        slong j;
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a;
        nmod_poly_t reduced_elem;
        fmpz_t coeff, den;
        ulong mod, d_mod, d_modinv;

        do {
            mod = n_randtest_not_zero(state);
        } while (mod == 1);

        fmpz_init(coeff);
        fmpz_init(den);

        nmod_poly_init(reduced_elem, mod);

        fmpq_poly_init(pol);
        do {
            fmpq_poly_randtest_not_zero(pol, state, 40, 200);
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);

        nf_elem_init(a, nf);

        do {
            nf_elem_randtest(a, state, 200, nf);
            nf_elem_get_den(den, a, nf);
            d_mod = fmpz_fdiv_ui(den, mod);
        } while (n_gcd(d_mod, mod) != 1);

        nf_elem_get_nmod_poly(reduced_elem, a, nf);

        for (j = 0; j < fmpq_poly_degree(pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            d_modinv = n_invmod(d_mod, mod);
            result = (nmod_poly_get_coeff_ui(reduced_elem, j) == nmod_mul(fmpz_fdiv_ui(coeff, mod), d_modinv, reduced_elem->mod));
            if (!result)
            {
                printf("FAIL: Reducing element with denominator\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); flint_printf("%u\n", mod);
                printf("a mod n = "); nmod_poly_print_pretty(reduced_elem, "x"); printf("\n");
                abort();
            }
        }

        fmpz_clear(den);
        fmpz_clear(coeff);
        nf_elem_clear(a, nf);
        nmod_poly_clear(reduced_elem);
        nf_clear(nf);
        fmpq_poly_clear(pol);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
