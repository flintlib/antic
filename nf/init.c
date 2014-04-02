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
       nf->flag = NF_MONIC;
    else
    {
       fmpz_preinvn_init(nf->pinv.qq, fmpq_poly_numref(pol) + pol->length - 1);
       nf->flag = NF_GENERIC;
    }

    if (pol->length == 3) /* quadratic case */
       nf->flag |= NF_QUADRATIC;

    if (pol->length <= NF_POWERS_CUTOFF && pol->length > 1) /* compute powers of generator mod pol */
    {
       if (nf->flag & NF_MONIC)
       {
          nf->powers.zz->powers = _fmpz_poly_powers_precompute(fmpq_poly_numref(pol), 
                                       pol->length);
          nf->powers.zz->len = pol->length;
       }
       else
       {
          nf->powers.qq->powers = _fmpq_poly_powers_precompute(fmpq_poly_numref(pol), 
                                       fmpq_poly_denref(pol), pol->length);
          nf->powers.qq->len = pol->length;
       }
   }
}

