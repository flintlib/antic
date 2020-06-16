/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
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

    flint_printf("mod_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        slong j;
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b;
        fmpz_t coeff, mod, reduced_coeff;

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, 2 * FLINT_BITS);
        fmpz_add_ui(mod, mod, 2);

        fmpz_init(coeff);
        fmpz_init(reduced_coeff);

        fmpq_poly_init(pol);
        do {
            fmpq_poly_randtest_not_zero(pol, state, 40, 200);
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);

        nf_elem_randtest(a, state, 200, nf);

        nf_elem_mod_fmpz_den(b, a, mod, nf, 0);

        for (j = 0; j < fmpq_poly_degree(pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            fmpz_mod(coeff, coeff, mod);
            nf_elem_get_coeff_fmpz(reduced_coeff, b, j, nf);
            result = fmpz_equal(reduced_coeff, coeff);
            if (!result)
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                abort();
            }
        }

        nf_elem_smod_fmpz_den(b, a, mod, nf, 0);

        for (j = 0; j < fmpq_poly_degree(pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            fmpz_smod(coeff, coeff, mod);
            nf_elem_get_coeff_fmpz(reduced_coeff, b, j, nf);
            result = fmpz_equal(reduced_coeff, coeff);
            if (!result)
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                abort();
            }
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        fmpz_clear(coeff);
        fmpz_clear(reduced_coeff);
        fmpz_clear(mod);
        nf_clear(nf);
        fmpq_poly_clear(pol);
    }

    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        slong j;
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;
        fmpz_t coeff, mod, den;

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, 2 * FLINT_BITS);
        fmpz_add_ui(mod, mod, 2);

        fmpz_init(coeff);
        fmpz_init(den);

        fmpq_poly_init(pol);
        do {
            fmpq_poly_randtest_not_zero(pol, state, 4, 2);
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(a, state, 2, nf);

        nf_elem_mod_fmpz(b, a, mod, nf);

        nf_elem_sub(c, b, a, nf);

        for (j = 0; j < fmpq_poly_degree(pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, c, j, nf);
            fmpz_mod(coeff, coeff, mod);
            result = fmpz_is_zero(coeff);
            if (!result || !nf_elem_den_is_one(c, nf))
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                abort();
            }
        }

        nf_elem_smod_fmpz(b, a, mod, nf);

        nf_elem_sub(c, b, a, nf);

        for (j = 0; j < fmpq_poly_degree(pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, c, j, nf);
            fmpz_smod(coeff, coeff, mod);
            result = fmpz_is_zero(coeff);
            if (!result || !nf_elem_den_is_one(c, nf))
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                abort();
            }
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        fmpz_clear(coeff);
        fmpz_clear(mod);
        fmpz_clear(den);
        nf_clear(nf);
        fmpq_poly_clear(pol);
    }

    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;
        fmpz_t mod, den, gcd;

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, 2 * FLINT_BITS);
        fmpz_add_ui(mod, mod, 2);

        fmpz_init(den);
        fmpz_init(gcd);

        fmpq_poly_init(pol);
        do {
            fmpq_poly_randtest_not_zero(pol, state, 4, 2);
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(a, state, 2, nf);

        nf_elem_coprime_den(b, a, mod, nf);

        nf_elem_sub(c, b, a, nf);

        nf_elem_get_den(den, c, nf);
        fmpz_gcd(gcd, den, mod);

        if (!fmpz_is_one(gcd))
        {
            printf("FAIL: Coprime denominators\n");
            printf("f = "); fmpq_poly_print_pretty(pol, "x"); printf("\n");
            printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
            printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
            printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
            printf("n = "); fmpz_print(mod); printf("\n");
            abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        fmpz_clear(mod);
        fmpz_clear(den);
        fmpz_clear(gcd);
        nf_clear(nf);
        fmpq_poly_clear(pol);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
