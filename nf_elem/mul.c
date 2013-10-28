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

           _fmpz_poly_divrem_preinv(q, NF_ELEM_NUMREF(a), plen, 
                  fmpq_poly_numref(nf->pol), nf->pinv.zz->coeffs, len);

           _fmpz_vec_clear(q, plen - len + 1);
           
           NF_ELEM(a)->length = len - 1;

           fmpq_poly_canonicalise(NF_ELEM(a));
        }
    }
    else
    {
        fmpq_poly_mul(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c));
        fmpq_poly_rem(NF_ELEM(a), NF_ELEM(a), nf->pol);
    }
}
