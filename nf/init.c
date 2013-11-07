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

    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart

******************************************************************************/

#include "nf.h"

void nf_init(nf_t nf, fmpq_poly_t pol)
{
    fmpq_poly_init(nf->pol);
    fmpq_poly_set(nf->pol, pol);

    if (fmpz_is_one(fmpq_poly_denref(pol)) /* denominator is one and numerator is monic */
     && fmpz_is_one(fmpq_poly_numref(pol) + pol->length - 1))
    {
       fmpz_poly_init2(nf->pinv.zz, pol->length);
       _fmpz_poly_preinvert(nf->pinv.zz->coeffs, fmpq_poly_numref(pol), pol->length);
       nf->pinv.zz->length = pol->length;

      nf->flag = NF_MONIC;
    } else
    {
       fmpz_preinvn_init(nf->pinv.qq, fmpq_poly_numref(pol) + pol->length - 1);
       nf->flag = NF_GENERIC;
    }

    if (pol->length <= NF_POWERS_CUTOFF && pol->length > 1) /* compute powers of generator mod pol */
    {
       slong i;
       fmpq_poly_t pow, p;
       
       fmpq_poly_init2(pow, pol->length);
       fmpq_poly_one(pow);
       fmpq_poly_init2(p, pol->length - 1);
       
       nf->powers = flint_malloc(sizeof(fmpq_poly_t)*(2*pol->length - 3));

       for (i = 0; i < 2*pol->length - 3; i++)
       {
          fmpq_poly_init(nf->powers + i);

          if (pow->length == pol->length) /* reduce power mod pol */
          {
             fmpz_set(fmpq_poly_denref(p), fmpq_poly_numref(pol) + pol->length - 1);
             _fmpz_vec_scalar_mul_fmpz(fmpq_poly_numref(p), fmpq_poly_numref(pol), 
                   pol->length - 1, fmpq_poly_numref(pow) + pow->length - 1);
             _fmpq_poly_set_length(p, pol->length - 1);
             _fmpq_poly_normalise(p);
             fmpq_poly_canonicalise(p);
             fmpq_poly_sub(pow, pow, p);
             _fmpq_poly_set_length(pow, pol->length - 1);
             _fmpq_poly_normalise(pow);
             fmpq_poly_canonicalise(pow);
          }

          fmpq_poly_set(nf->powers + i, pow);
          fmpq_poly_shift_left(pow, pow, 1);
       }

       fmpq_poly_clear(p);
       fmpq_poly_clear(pow);
    }
}

