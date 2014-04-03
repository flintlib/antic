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

void _nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can)
{
   fmpz_t d;

   const fmpz * const bnum = QNF_ELEM_NUMREF(b);
   const fmpz * const bden = QNF_ELEM_DENREF(b);
   
   const fmpz * const cnum = QNF_ELEM_NUMREF(c);
   const fmpz * const cden = QNF_ELEM_DENREF(c);
   
   fmpz * const anum = QNF_ELEM_NUMREF(a);
   fmpz * const aden = QNF_ELEM_DENREF(a);

   fmpz_init(d);
   fmpz_one(d);

   if (fmpz_equal(bden, cden))
   {
      fmpz_add(anum, bnum, cnum);
      fmpz_add(anum + 1, bnum + 1, cnum + 1);
      fmpz_set(aden, bden);

      if (can && !fmpz_is_one(aden))
      {
         fmpz_gcd(d, anum, anum + 1);
         if (!fmpz_is_one(d))
         {
            fmpz_gcd(d, d, aden);

            if (!fmpz_is_one(d))
            {
               fmpz_divexact(anum, anum, d);
               fmpz_divexact(anum + 1, anum + 1, d);
               fmpz_divexact(aden, aden, d);
            }
         }
      }

      fmpz_clear(d);

      return;
   }

   if (!fmpz_is_one(bden) && !fmpz_is_one(cden))
      fmpz_gcd(d, bden, cden);

   if (fmpz_is_one(d))
   {
      fmpz_mul(anum, bnum, cden);
      fmpz_mul(anum + 1, bnum + 1, cden);
      fmpz_addmul(anum, cnum, bden);
      fmpz_addmul(anum + 1, cnum + 1, bden);
      fmpz_mul(aden, bden, cden);
   } else
   {
      fmpz_t bden1;
      fmpz_t cden1;
      
      fmpz_init(bden1);
      fmpz_init(cden1);
      
      fmpz_divexact(bden1, bden, d);
      fmpz_divexact(cden1, cden, d);
        
      fmpz_mul(anum, bnum, cden1);
      fmpz_mul(anum + 1, bnum + 1, cden1);
      fmpz_addmul(anum, cnum, bden1);
      fmpz_addmul(anum + 1, cnum + 1, bden1);
      
      if (fmpz_is_zero(anum) && fmpz_is_zero(anum + 1))
         fmpz_one(aden);
      else
      {
         if (can)
         {
            fmpz_t e;
            
            fmpz_init(e);
              
            fmpz_gcd(e, anum, anum + 1);
            if (!fmpz_is_one(e))
               fmpz_gcd(e, e, d);
            
            if (fmpz_is_one(e))
               fmpz_mul(aden, bden, cden1);
            else
            {
                fmpz_divexact(anum, anum, e);
                fmpz_divexact(anum + 1, anum + 1, e);
                fmpz_divexact(bden1, bden, e);
                fmpz_mul(aden, bden1, cden1);
            }
            
            fmpz_clear(e);
         } else
            fmpz_mul(aden, bden, cden1);
      }

      fmpz_clear(bden1);
      fmpz_clear(cden1);
   }

   fmpz_clear(d);
}

void nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
   if (a == c)
   {
      nf_elem_t t;

      nf_elem_init(t, nf);

      _nf_elem_add_qf(t, b, c, nf, 1);
      nf_elem_swap(t, a, nf);

      nf_elem_clear(t, nf);
   } else
      _nf_elem_add_qf(a, b, c, nf, 1);
}
