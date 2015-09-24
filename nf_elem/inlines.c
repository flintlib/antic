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

    Copyright (C) 2015 William Hart

******************************************************************************/

#define NF_ELEM_INLINES_C

#include "nf_elem.h"

void nf_elem_get_den(fmpz * d, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
     fmpz_set(d, LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {  
     fmpz_set(d, QNF_ELEM_DENREF(b));
   } else
   {
     fmpz_set(d, NF_ELEM_DENREF(b));
   }
}

NF_ELEM_INLINE
void nf_elem_set_den(nf_elem_t b, fmpz * d, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
     fmpz_set(LNF_ELEM_DENREF(b), d);
   } else if (nf->flag & NF_QUADRATIC)
   {  
     fmpz_set(QNF_ELEM_DENREF(b), d);
   } else
   {
     fmpz_set(NF_ELEM_DENREF(b), d);
   }
}

