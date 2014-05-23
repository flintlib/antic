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

    Copyright (C) 2014 William Hart

******************************************************************************/

#include "nf_elem_approx.h"

void __to_approx(acb_ptr b, const fmpz * a, const nf_t nf)
{
    acb_mat_t A, B;
    slong i, n, prec;

    n = acb_mat_nrows(nf->V);
    prec = nf->Vprec;

    acb_mat_init(A, n, 1);
    acb_mat_init(B, n, 1);

    for (i = 0; i < n; i++)
        acb_set_round_fmpz(acb_mat_entry(A, i, 0), a + i, prec);

    acb_mat_mul(B, nf->V, A, prec);

    for (i = 0; i < n; i++)
        acb_set(b + i, acb_mat_entry(B, i, 0));

    acb_mat_clear(A);
    acb_mat_clear(B);
}

void nf_elem_approx_to_approx(nf_elem_approx_t a, const nf_t nf)
{
   __to_approx(a->conj, NF_ELEM_NUMREF(a), nf);
}