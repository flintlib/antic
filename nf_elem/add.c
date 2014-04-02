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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 William Hart

******************************************************************************/

#include "nf_elem.h"

void nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can)
{
   fmpz_t d;

   const fmpz * a1 = QNF_ELEM(b)->a;
   const fmpz * b1 = QNF_ELEM(b)->b;
   const fmpz * den1 = QNF_ELEM(b)->den;
   
   const fmpz * a2 = QNF_ELEM(c)->a;
   const fmpz * b2 = QNF_ELEM(c)->b;
   const fmpz * den2 = QNF_ELEM(c)->den;
   
   fmpz * a3 = QNF_ELEM(a)->a;
   fmpz * b3 = QNF_ELEM(a)->b;
   fmpz * den3 = QNF_ELEM(a)->den;

   fmpz_init(d);
   fmpz_one(d);

   if (fmpz_equal(den1, den2))
   {
      fmpz_add(a3, a1, a2);
      fmpz_add(b3, b1, b2);
      fmpz_set(den3, den1);

      if (can && !fmpz_is_one(den3))
      {
         fmpz_gcd(d, a3, b3);
         if (!fmpz_is_one(d))
         {
            fmpz_gcd(d, d, den3);

            if (!fmpz_is_one(d))
            {
               fmpz_divexact(a3, a3, d);
               fmpz_divexact(b3, b3, d);
               fmpz_divexact(den3, den3, d);
            }
         }

         if (fmpz_sgn(den3) < 0)
         {
            nf_elem_neg(a, a, nf);
            fmpz_neg(den3, den3);
         }
      }

      fmpz_clear(d);

      return;
   }

   if (!fmpz_is_one(den1) && !fmpz_is_one(den2))
      fmpz_gcd(d, den1, den2);

   if (fmpz_is_one(d))
   {
      fmpz_mul(a3, a1, den2);
      fmpz_mul(b3, b1, den2);
      fmpz_addmul(a3, a2, den1);
      fmpz_addmul(b3, b2, den1);
      fmpz_mul(den3, den1, den2);
   } else
   {
      fmpz_t den11;
      fmpz_t den22;
      
      fmpz_init(den11);
      fmpz_init(den22);
      
      fmpz_divexact(den11, den1, d);
      fmpz_divexact(den22, den2, d);
        
      fmpz_mul(a3, a1, den22);
      fmpz_mul(b3, b1, den22);
      fmpz_addmul(a3, a2, den11);
      fmpz_addmul(b3, b2, den11);
      
      if (fmpz_is_zero(a3) && fmpz_is_zero(b3))
         fmpz_one(den3);
      else
      {
         if (can)
         {
            fmpz_t e;
            
            fmpz_init(e);
              
            fmpz_gcd(e, a3, b3);
            if (!fmpz_is_one(e))
               fmpz_gcd(e, e, d);
            
            if (fmpz_is_one(e))
               fmpz_mul(den3, den1, den22);
            else
            {
                fmpz_divexact(a3, a3, e);
                fmpz_divexact(b3, b3, e);
                fmpz_divexact(den11, den1, e);
                fmpz_mul(den3, den11, den22);
            }
            
            fmpz_clear(e);
         } else
            fmpz_mul(den3, den1, den22);
      }

      fmpz_clear(den11);
      fmpz_clear(den22);
   }

   fmpz_clear(d);
}
