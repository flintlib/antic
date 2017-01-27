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

#include "nf_elem.h"

void nf_elem_si_sub(nf_elem_t a, slong c, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz * den = LNF_ELEM_DENREF(a);
	  fmpz * num = LNF_ELEM_NUMREF(a);
	  
      nf_elem_set(a, b, nf);
	  
	  fmpz_neg(num, num);
	  if (c >= 0)
	     fmpz_addmul_ui(num, den, c);
	  else
	     fmpz_submul_ui(num, den, -c);
	  _fmpq_canonicalise(num, den);
   }
   else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * den = QNF_ELEM_DENREF(a);
	  fmpz * num = QNF_ELEM_NUMREF(a); 
	  slong len = 2, i;
	  
	  nf_elem_set(a, b, nf);
	  
	  while (len != 0 && fmpz_is_zero(num + len - 1))
	     len--;
	  
      if (c >= 0)
	     fmpz_submul_ui(num, den, c);
	  else
	     fmpz_addmul_ui(num, den, -c);

          for(i=0; i<= len; i++) 
             fmpz_neg(num+i, num+i);
	  _fmpq_poly_canonicalise(num, den, len);
   } else
   {
      fmpq_poly_si_sub(NF_ELEM(a), c, NF_ELEM(b));
   }
}
