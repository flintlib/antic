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

void nf_elem_mul(nf_elem_t a, nf_elem_t b, nf_elem_t c, nf_t nf)
{
    const slong len1 = NF_ELEM(b)->length;
    const slong len2 = NF_ELEM(c)->length;
        
    if (len1 == 0 || len2 == 0)
    {
       nf_elem_zero(a, nf);

       return;
    }

    if (nf->flag & NF_MONIC)
    {
        slong plen = len1 + len2 - 1;
        const slong len = nf->pol->length;
        
        fmpq_poly_mul(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c));
       
        if (plen > len - 1)
        {
           fmpz * q = _fmpz_vec_init(plen - len + 1);
           fmpz * r = _fmpz_vec_init(plen);
           _fmpz_vec_set(r, NF_ELEM_NUMREF(a), plen);

           _fmpz_poly_divrem(q, NF_ELEM_NUMREF(a), r, plen, 
                  fmpq_poly_numref(nf->pol), len);

           _fmpz_vec_clear(r, plen);
           _fmpz_vec_clear(q, plen - len + 1);
           
           NF_ELEM(a)->length = len - 1;

           fmpq_poly_canonicalise(NF_ELEM(a));
        }
    }
    else
    {
        fmpq_poly_t t;
        
        fmpq_poly_mul(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c));
        
        if (NF_ELEM(a)->length > nf->pol->length - 1)
        {
           if (nf->pol->length <= NF_POWERS_CUTOFF)
           {
              fmpq_poly_t prod;
              fmpz_t den;
              slong i;

              fmpz_init(den);
              fmpq_poly_init2(prod, nf->pol->length - 1);
             
              fmpz_set(den, fmpq_poly_denref(NF_ELEM(a)));
              
              for (i = nf->pol->length - 1; i < NF_ELEM(a)->length; i++)
              {
                 _fmpz_vec_scalar_mul_fmpz(fmpq_poly_numref(prod), 
                    fmpq_poly_numref(nf->powers + i),
                    nf->powers[i].length,
                    NF_ELEM_NUMREF(a) + i);
                 fmpz_mul(fmpq_poly_denref(prod), fmpq_poly_denref(nf->powers + i), den);
                 _fmpq_poly_add_can(NF_ELEM_NUMREF(a), fmpq_poly_denref(NF_ELEM(a)),
                    NF_ELEM_NUMREF(a), fmpq_poly_denref(NF_ELEM(a)), nf->pol->length - 1,
                    fmpq_poly_numref(prod), fmpq_poly_denref(prod), nf->powers[i].length, 0);
              }

              _fmpq_poly_set_length(NF_ELEM(a), nf->pol->length - 1);
              _fmpq_poly_normalise(NF_ELEM(a));
              fmpq_poly_canonicalise_weak(NF_ELEM(a));

              fmpz_clear(den);
              fmpq_poly_clear(prod);
           } else
           {
              fmpq_poly_init2(t, NF_ELEM(a)->length);
        
              _fmpq_poly_rem(t->coeffs, t->den,
                 NF_ELEM(a)->coeffs, NF_ELEM(a)->den, NF_ELEM(a)->length, 
                 nf->pol->coeffs, nf->pol->den, nf->pol->length, nf->pinv.qq); 
           
              _fmpq_poly_set_length(t, nf->pol->length - 1);
              _fmpq_poly_normalise(t);
        
              fmpq_poly_swap(t, NF_ELEM(a));
              fmpq_poly_clear(t);
           }
        }
    }
}
