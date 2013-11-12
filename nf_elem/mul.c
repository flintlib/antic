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

void _nf_elem_mul(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
    const slong len1 = NF_ELEM(b)->length;
    const slong len2 = NF_ELEM(c)->length;
    const slong len = nf->pol->length;
    slong plen = len1 + len2 - 1;
      
    if (len1 == 0 || len2 == 0)
    {
       nf_elem_zero(a, nf);

       return;
    }

    if (len1 >= len2)
    {
       _fmpz_poly_mul(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(b), len1,
          NF_ELEM_NUMREF(c), len2);
    }
    else
    {
        _fmpz_poly_mul(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(c), len2,
           NF_ELEM_NUMREF(b), len1);
    }

    fmpz_mul(fmpq_poly_denref(NF_ELEM(a)), fmpq_poly_denref(NF_ELEM(b)),
       fmpq_poly_denref(NF_ELEM(c)));

    _fmpq_poly_set_length(NF_ELEM(a), plen);

    if (nf->flag & NF_MONIC)
    {
        if (plen > len - 1)
        {
           if (len <= NF_POWERS_CUTOFF)
           {
              _fmpz_poly_rem_powers_precomp(NF_ELEM_NUMREF(a), plen,
                 fmpq_poly_numref(nf->pol), len, nf->powers.zz->powers);

              _fmpq_poly_set_length(NF_ELEM(a), len - 1);
              _fmpq_poly_normalise(NF_ELEM(a));
              
           } else
           {
              fmpz * q = _fmpz_vec_init(plen - len + 1);
              fmpz * r = _fmpz_vec_init(plen);
              _fmpz_vec_set(r, NF_ELEM_NUMREF(a), plen);

              _fmpz_poly_divrem(q, NF_ELEM_NUMREF(a), r, plen, 
                  fmpq_poly_numref(nf->pol), len);

              _fmpz_vec_clear(r, plen);
              _fmpz_vec_clear(q, plen - len + 1);
           
              NF_ELEM(a)->length = len - 1;
           }
        }
    }
    else
    {
        fmpq_poly_t t;
        
        if (NF_ELEM(a)->length >= len)
        {
           if (nf->pol->length <= NF_POWERS_CUTOFF)
           {
              _fmpq_poly_rem_powers_precomp(NF_ELEM_NUMREF(a), 
                 fmpq_poly_denref(NF_ELEM(a)), plen,
                 fmpq_poly_numref(nf->pol), fmpq_poly_denref(nf->pol), 
                 len, nf->powers.qq->powers);

              _fmpq_poly_set_length(NF_ELEM(a), len - 1);
              _fmpq_poly_normalise(NF_ELEM(a));
           } else
           {
              fmpq_poly_init2(t, NF_ELEM(a)->length);
        
              _fmpq_poly_rem(t->coeffs, t->den,
                 NF_ELEM(a)->coeffs, NF_ELEM(a)->den, NF_ELEM(a)->length, 
                 nf->pol->coeffs, nf->pol->den, len, nf->pinv.qq); 
           
              _fmpq_poly_set_length(t, len - 1);
              _fmpq_poly_normalise(t);
        
              fmpq_poly_swap(t, NF_ELEM(a));
              fmpq_poly_clear(t);
           }
        }
    }
}

void nf_elem_mul(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
   nf_elem_t t;
   
   if (a == b || a == c)
   {
      nf_elem_init(t, nf);

      _nf_elem_mul(t, b, c, nf);
      fmpq_poly_swap(NF_ELEM(t), NF_ELEM(a));

      nf_elem_clear(t, nf);
   }
   else
      _nf_elem_mul(a, b, c, nf);

   fmpq_poly_canonicalise(NF_ELEM(a));
}