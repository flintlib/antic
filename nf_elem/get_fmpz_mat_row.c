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

void nf_elem_get_fmpz_mat_row(fmpz_mat_t M, const int i, const nf_elem_t b, const nf_t nf)
{
  if (nf->flag & NF_LINEAR)
  {
    fmpz_set(fmpz_mat_entry(M, i, 0), LNF_ELEM_NUMREF(b));
  } else if (nf->flag & NF_QUADRATIC)
  {
    const fmpz * const bnum = QNF_ELEM_NUMREF(b);
    fmpz_set(fmpz_mat_entry(M, i, 0), bnum);
    fmpz_set(fmpz_mat_entry(M, i, 1), bnum + 1);
  } else
  {
    int j;
    for (j=0; j< NF_ELEM(b)->length; j++) 
    {
      fmpz_set(fmpz_mat_entry(M, i, j), NF_ELEM_NUMREF(b) + j);
    }
    for (; j< nf->pol->length-1; j++)
    {
      fmpz_zero(fmpz_mat_entry(M, i, j));
    }
  }
}
