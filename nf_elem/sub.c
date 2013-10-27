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

void nf_elem_sub(nf_elem_t a, nf_elem_t b, nf_elem_t c, nf_t nf)
{
    if (nf->flag & NF_MONIC)
    {
        const slong len1 = NF_ELEM(b)->length;
        const slong len2 = NF_ELEM(c)->length;
        slong len;

        _fmpz_poly_sub(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(b), len1,
                                          NF_ELEM_NUMREF(c), len2);

        len = FLINT_MAX(len1, len2);

        /* normalise */
        while (len && fmpz_is_zero(NF_ELEM_NUMREF(a) + len - 1)) len--;
        NF_ELEM(a)->length = len;
    }
    else
    {
        fmpq_poly_sub(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c));
    }
}
