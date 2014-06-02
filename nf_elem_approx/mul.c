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

#include "nf_elem_approx.h"

void nf_elem_approx_mul(nf_elem_approx_t a, nf_elem_approx_t b, 
                                            nf_elem_approx_t c, const nf_t nf)
{
   slong i, deg = fmpq_poly_degree(nf->pol);
   fmpz_t pow;
   fmpz * lead;

   fmpz_mul(NF_ELEM_DENREF(a), 
      NF_ELEM_DENREF(b), NF_ELEM_DENREF(c));
   
   for (i = 0; i < deg; i++)
      acb_mul(a->conj + i, b->conj + i, c->conj + i, nf->Vprec);

   lead = fmpq_poly_numref(nf->pol) + deg;

   if (!fmpz_is_one(lead)) /* non-monic defining poly */
   {
      fmpz_init(pow);
      fmpz_pow_ui(pow, lead, deg - 1);

      fmpz_mul(NF_ELEM_DENREF(a), NF_ELEM_DENREF(a), pow);

      for (i = 0; i < deg; i++)
         acb_mul_fmpz(a->conj + i, a->conj + i, pow, nf->Vprec);

      fmpz_clear(pow);
   } 
}