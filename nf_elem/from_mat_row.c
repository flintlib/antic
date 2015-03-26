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
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 William Hart
    Copyright (C) 2015 Claus Fieker

******************************************************************************/

#include "nf_elem.h"
#include "fmpq_poly.h"

void nf_elem_from_mat_row(nf_elem_t b, const fmpz_mat_t M, const fmpz * d, const int i, const nf_t nf)
{
  if (nf->flag & NF_LINEAR)
  {
    fmpz_set(LNF_ELEM_DENREF(b), d);
    fmpz_set(LNF_ELEM_NUMREF(b), fmpz_mat_entry(M, i, 0));
  } else if (nf->flag & NF_QUADRATIC)
  {
    fmpz * const bnum = QNF_ELEM_NUMREF(b);
    fmpz_set(QNF_ELEM_DENREF(b), d);
    fmpz_set(bnum, fmpz_mat_entry(M, i, 0));
    fmpz_set(bnum + 1, fmpz_mat_entry(M, i, 1));
  } else
  {
    int j;
    fmpz_set(NF_ELEM_DENREF(b), d);
    for (j=nf->pol->length-2; j>=0; j--)
      if (!fmpz_is_zero(fmpz_mat_entry(M, i, j)))
        break;
    _fmpq_poly_set_length(NF_ELEM(b), j+1);
    for (; j>=0; j--)
      fmpq_poly_set_coeff_fmpz(NF_ELEM(b),  j, fmpz_mat_entry(M, i, j));
  }
}
