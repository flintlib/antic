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

void nf_elem_approx_add(nf_elem_approx_t a, nf_elem_approx_t b, 
                                            nf_elem_approx_t c, const nf_t nf)
{
   fmpz_t g, b2, c2;
   slong i, deg = fmpq_poly_degree(nf->pol);
   acb_ptr t;

   fmpz_init(g);
   fmpz_init_set_ui(b2, 1);
   fmpz_init_set_ui(c2, 1);

   if (a == b || a == c)
      t = _acb_vec_init(deg);
   else
      t = a->conj;

   if (!fmpz_is_one(NF_ELEM_DENREF(b)) && !fmpz_is_one(NF_ELEM_DENREF(c)))
      fmpz_gcd(g, NF_ELEM_DENREF(b), NF_ELEM_DENREF(c));
   else
      fmpz_set_ui(g, 1);

   if (!fmpz_equal(NF_ELEM_DENREF(b), NF_ELEM_DENREF(c)))
   {
      if (!fmpz_is_one(g))
      {
         fmpz_divexact(b2, NF_ELEM_DENREF(b), g);
         fmpz_divexact(c2, NF_ELEM_DENREF(c), g);
      } else
      {
         fmpz_set(b2, NF_ELEM_DENREF(b));
         fmpz_set(c2, NF_ELEM_DENREF(c));
      }
   }

   if (fmpz_is_one(c2))
   {
      if (fmpz_is_one(b2))
      {
         for (i = 0; i < deg; i++)
            acb_add(t + i, b->conj + i, c->conj + i, nf->Vprec);
      } else
      {
         for (i = 0; i < deg; i++)
         {
            acb_set(t + i, b->conj + i);
            acb_addmul_fmpz(t + i, c->conj + i, b2, nf->Vprec);
         }
      }
   } else
   {
      if (fmpz_is_one(b2))
      {
         for (i = 0; i < deg; i++)
         {
            acb_set(t + i, c->conj + i);
            acb_addmul_fmpz(t + i, b->conj + i, c2, nf->Vprec);
         }
      } else
      {
         for (i = 0; i < deg; i++)
         {
            acb_mul_fmpz(t + i, b->conj + i, c2, nf->Vprec);
            acb_addmul_fmpz(t + i, c->conj + i, b2, nf->Vprec);
         }
      }
   }

   fmpz_mul(NF_ELEM_DENREF(a), 
      NF_ELEM_DENREF(b), c2);

   if (a == b || a == c)
   {
      _acb_vec_clear(a->conj, deg);
      a->conj = t;
   }

   fmpz_clear(b2);
   fmpz_clear(c2);
   fmpz_clear(g);
}