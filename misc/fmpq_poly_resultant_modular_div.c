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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2015 Claus Fieker

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_resultant_div(fmpz_t rnum, fmpz_t rden, 
                          const fmpz *poly1, const fmpz_t den1, slong len1, 
                          const fmpz *poly2, const fmpz_t den2, slong len2,
                          const fmpz_t divisor, slong num_primes)
{
    fmpz_t div, l;

    if (len2 == 1)
    {
        if (len1 == 1)
        {
            fmpz_one(rnum);
            fmpz_one(rden);
        }
        else if (len1 == 2)
        {
            fmpz_set(rnum, poly2);
            fmpz_set(rden, den2);
        }
        else
        {
            fmpz_pow_ui(rnum, poly2, len1 - 1);
            if (fmpz_is_one(den2))
            {
                fmpz_one(rden);
            }
            else
            {
                fmpz_pow_ui(rden, den2, len1 - 1);
            }
        }
        fmpz_divexact(rnum, rnum, divisor);
    }
    else  /* len1 >= len2 >= 2 */
    {
        fmpz_t c1, c2, t;
        fmpz *prim1, *prim2, *g;
        slong lenG = len2;

        fmpz_init(c1);
        fmpz_init(c2);

        _fmpz_vec_content(c1, poly1, len1);
        _fmpz_vec_content(c2, poly2, len2);

        prim1 = _fmpz_vec_init(len1);
        prim2 = _fmpz_vec_init(len2);

        _fmpz_vec_scalar_divexact_fmpz(prim1, poly1, len1, c1);
        _fmpz_vec_scalar_divexact_fmpz(prim2, poly2, len2, c2);

        fmpz_init(l);
        if (!fmpz_is_one(c1))
        {
           fmpz_pow_ui(l, c1, len2 - 1);
           fmpz_init(div);
           fmpz_gcd(div, l, divisor);
           fmpz_divexact(div, divisor, div);
        } else {
           fmpz_init_set(div, divisor);
        }

        if (!fmpz_is_one(c2))
        {
           fmpz_pow_ui(l, c2, len1 - 1);
           fmpz_gcd(l, l, div);
           fmpz_divexact(div, div, l);
        }


        _fmpz_poly_resultant_modular_div(rnum, prim1, len1, prim2, len2, div, num_primes);

        fmpz_init(t);
        if (!fmpz_is_one(c1))
        {
            fmpz_pow_ui(t, c1, len2 - 1);
            fmpz_mul(rnum, rnum, t);
        }
        if (!fmpz_is_one(c2))
        {
            fmpz_pow_ui(t, c2, len1 - 1);
            fmpz_mul(rnum, rnum, t);
        }

        if (fmpz_is_one(den1))
        {
            if (fmpz_is_one(den2))
                fmpz_one(rden);
            else
                fmpz_pow_ui(rden, den2, len1 - 1);
        }
        else
        {
            if (fmpz_is_one(den2))
                fmpz_pow_ui(rden, den1, len2 - 1);
            else
            {
                fmpz_pow_ui(rden, den1, len2 - 1);
                fmpz_pow_ui(t,    den2, len1 - 1);
                fmpz_mul(rden, rden, t);
            }
        }
        _fmpq_canonicalise(rnum, rden);
        fmpz_clear(t);

        fmpz_clear(c1);
        fmpz_clear(c2);
        fmpz_clear(div);
        fmpz_clear(l);
        _fmpz_vec_clear(prim1, len1);
        _fmpz_vec_clear(prim2, len2);
    }
}

void fmpq_poly_resultant_div(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g, const fmpz_t divisor, slong num_primes)
{
    const slong len1 = f->length;
    const slong len2 = g->length;

    if (len1 == 0 || len2 == 0)
    {
        fmpq_zero(r);
    }
    else
    {
        if (len1 >= len2)
        {
            _fmpq_poly_resultant_div(fmpq_numref(r), fmpq_denref(r), 
                                 f->coeffs, f->den, len1, 
                                 g->coeffs, g->den, len2,
                                 divisor, num_primes);
        }
        else
        {
            _fmpq_poly_resultant_div(fmpq_numref(r), fmpq_denref(r), 
                                 g->coeffs, g->den, len2, 
                                 f->coeffs, f->den, len1,
                                 divisor, num_primes);

            if (((len1 | len2) & WORD(1)) == WORD(0))
                fmpq_neg(r, r);
        }
    }
}

