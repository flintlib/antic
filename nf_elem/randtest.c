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

#include "nf_elem.h"
#include "ulong_extras.h"

void nf_elem_randtest(nf_elem_t a, flint_rand_t state, 
                                               mp_bitcnt_t bits, const nf_t nf)
{
    if (nf->flag & NF_QUADRATIC)
    {
        fmpz_randtest(QNF_ELEM(a)->a, state, bits);
        fmpz_randtest(QNF_ELEM(a)->b, state, bits);

        if (n_randint(state, 2))
        {
           fmpz_t d;
           
           fmpz_randtest_not_zero(QNF_ELEM(a)->den, state, bits);
           fmpz_abs(QNF_ELEM(a)->den, QNF_ELEM(a)->den);

           fmpz_init(d);
           fmpz_gcd(d, QNF_ELEM(a)->a, QNF_ELEM(a)->b);
           if (!fmpz_is_one(d))
           {
              fmpz_gcd(d, d, QNF_ELEM(a)->den);

              if (!fmpz_is_one(d))
              {
                 fmpz_divexact(QNF_ELEM(a)->a, QNF_ELEM(a)->a, d);
                 fmpz_divexact(QNF_ELEM(a)->b, QNF_ELEM(a)->b, d);
                 fmpz_divexact(QNF_ELEM(a)->den, QNF_ELEM(a)->den, d);
              }
           }
        } else
           fmpz_one(QNF_ELEM(a)->den);
    }
    else
    {
        fmpq_poly_randtest(NF_ELEM(a), state, nf->pol->length - 1, bits);
    }
}

void nf_elem_randtest_not_zero(nf_elem_t a, flint_rand_t state, 
                                               mp_bitcnt_t bits, const nf_t nf)
{
   if (nf->flag & NF_QUADRATIC)
   {
       do {
          nf_elem_randtest(a, state, bits, nf);
       } while (fmpz_is_zero(QNF_ELEM(a)->a) && fmpz_is_zero(QNF_ELEM(a)->b));
   } else
   {
      do {
          nf_elem_randtest(a, state, bits, nf);
       } while (fmpq_poly_is_zero(NF_ELEM(a)));
   }
}