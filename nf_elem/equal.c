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

int nf_elem_equal(nf_elem_t a, nf_elem_t b, nf_t nf)
{
    const slong len1 = NF_ELEM(a)->length;
    const slong len2 = NF_ELEM(b)->length;
        
    if (len1 != len2)
       return 0;

    if (fmpz_equal(fmpq_poly_denref(NF_ELEM(a)), fmpq_poly_denref(NF_ELEM(b))))
    {
       if (nf->flag & NF_MONIC)
          return _fmpz_vec_equal(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(b), len1);
       else
          return fmpq_poly_equal(NF_ELEM(a), NF_ELEM(b));
    } else
    {
        slong i;
        slong d = fmpz_bits(fmpq_poly_denref(NF_ELEM(b)))
                - fmpz_bits(fmpq_poly_denref(NF_ELEM(a))) + 1;
        fmpz * p1 = NF_ELEM_NUMREF(a);
        fmpz * p2 = NF_ELEM_NUMREF(b);

        for (i = 0; i < len1; i++)
        {
           slong b1 = fmpz_bits(p1 + i);
           slong b2 = fmpz_bits(p2 + i);
           if (!(b1 == 0 && b2 == 0) && (ulong) (b1 - b2 + d) > 2)
              return 0;
        }

        fmpq_poly_canonicalise(NF_ELEM(a));
        fmpq_poly_canonicalise(NF_ELEM(b));
        
        return fmpq_poly_equal(NF_ELEM(a), NF_ELEM(b));
    }
}
