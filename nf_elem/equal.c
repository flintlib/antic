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

    Copyright (C) 2013 William Hart

******************************************************************************/

#include "nf_elem.h"

int _nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_QUADRATIC)
   {
      slong d, bits1, bits2;
      int res = 1;

      const fmpz * a1 = QNF_ELEM(a)->a;
      const fmpz * b1 = QNF_ELEM(a)->b;
   
      const fmpz * a2 = QNF_ELEM(b)->a;
      const fmpz * b2 = QNF_ELEM(b)->b;

      const fmpz * den1 = QNF_ELEM(a)->den;
      const fmpz * den2 = QNF_ELEM(b)->den;

      fmpz_t t1, t2;

      if (fmpz_equal(den1, den2))
         return fmpz_equal(a1, a2) && fmpz_equal(b1, b2);
      
      d = fmpz_bits(den1) - fmpz_bits(den2) + 1;
      
      bits1 = fmpz_bits(b1);
      bits2 = fmpz_bits(b2);
      if (!(bits1 == 0 && bits2 == 0) && (ulong) (bits1 - bits2 + d) > 2)
         return 0;

      bits1 = fmpz_bits(a1);
      bits2 = fmpz_bits(a2);
      if (!(bits1 == 0 && bits2 == 0) && (ulong) (bits1 - bits2 + d) > 2)
         return 0;

      fmpz_init(t1);
      fmpz_init(t2);

      fmpz_mul(t1, a1, den2);
      fmpz_mul(t2, a2, den1);

      if (!fmpz_equal(t1, t2))
      {
         res = 0;
         goto cleanup;
      }

      fmpz_mul(t1, b1, den2);
      fmpz_mul(t2, b2, den1);

      if (!fmpz_equal(t1, t2))
      {
         res = 0;
         goto cleanup;
      }

cleanup:

      fmpz_clear(t1);
      fmpz_clear(t2);

      return res;
   } else
   {  
      const slong len1 = NF_ELEM(a)->length;
      const slong len2 = NF_ELEM(b)->length;
        
      if (len1 != len2)
         return 0;

      if (fmpz_equal(fmpq_poly_denref(NF_ELEM(a)), fmpq_poly_denref(NF_ELEM(b))))
         return _fmpz_vec_equal(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(b), len1);
      else
      {
          slong i;
          slong d = fmpz_bits(fmpq_poly_denref(NF_ELEM(b)))
                  - fmpz_bits(fmpq_poly_denref(NF_ELEM(a))) + 1;
          fmpz * p1 = NF_ELEM_NUMREF(a);
          fmpz * p2 = NF_ELEM_NUMREF(b);
          fmpz_t gcd, den1, den2;
          fmpz * t1, * t2;
          int equal;
  
          for (i = 0; i < len1; i++)
          {
             slong b1 = fmpz_bits(p1 + i);
             slong b2 = fmpz_bits(p2 + i);
             if (!(b1 == 0 && b2 == 0) && (ulong) (b1 - b2 + d) > 2)
                return 0;
          }
 
          fmpz_init(gcd);
          fmpz_init(den1);
          fmpz_init(den2);

          /* TODO: possibly only compute GCD if it will save time */
          fmpz_gcd(gcd, fmpq_poly_denref(NF_ELEM(a)), fmpq_poly_denref(NF_ELEM(b)));
          fmpz_divexact(den1, fmpq_poly_denref(NF_ELEM(a)), gcd);
          fmpz_divexact(den2, fmpq_poly_denref(NF_ELEM(b)), gcd);
  
          t1 = _fmpz_vec_init(len1);
          t2 = _fmpz_vec_init(len1);

          _fmpz_vec_scalar_mul_fmpz(t1, p1, len1, den2);
          _fmpz_vec_scalar_mul_fmpz(t2, p2, len2, den1);

          equal = _fmpz_vec_equal(t1, t2, len1);

          fmpz_clear(gcd);
          fmpz_clear(den1);
          fmpz_clear(den2);

          _fmpz_vec_clear(t1, len1);
          _fmpz_vec_clear(t2, len1);
 
          return equal;
      }
   }
}

int nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_QUADRATIC)
   {
      if (!fmpz_equal(QNF_ELEM(a)->den, QNF_ELEM(b)->den))
         return 0;

      if (!fmpz_equal(QNF_ELEM(a)->a, QNF_ELEM(b)->a))
         return 0;

      if (!fmpz_equal(QNF_ELEM(a)->b, QNF_ELEM(b)->b))
         return 0;

      return 1;
   } else
   {
      const slong len1 = NF_ELEM(a)->length;
      const slong len2 = NF_ELEM(b)->length;
        
      if (len1 != len2)
         return 0;

      if (fmpz_equal(fmpq_poly_denref(NF_ELEM(a)), fmpq_poly_denref(NF_ELEM(b))))
         return _fmpz_vec_equal(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(b), len1);
      else
         return 0;
   }
}