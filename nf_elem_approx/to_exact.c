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

void __to_exact(fmpz * a, acb_srcptr b, const nf_t nf)
{
    acb_mat_t A, B;
    slong i, n, prec;

    n = fmpq_poly_degree(nf->pol);
    prec = nf->Vprec;

    acb_mat_init(A, n, 1);
    acb_mat_init(B, n, 1);

    for (i = 0; i < n; i++)
        acb_set(acb_mat_entry(B, i, 0), b + i);

    acb_mat_mul(A, nf->Vinv, B, prec);

    for (i = 0; i < n; i++)
    {
        if (!arb_contains_zero(acb_imagref(acb_mat_entry(A, i, 0))) ||
            !arb_get_unique_fmpz(a + i, acb_realref(acb_mat_entry(A, i, 0))))
        {
            printf("a most grave fault happened! the precision had insufficient greatness to facilitate an unambiguous retrieval of exact information\n");
            fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
            abort();
        }
    }

    acb_mat_clear(A);
    acb_mat_clear(B);
}

void nf_elem_approx_to_exact(nf_elem_approx_t a, const nf_t nf)
{
   __to_exact(NF_ELEM_NUMREF(a), a->conj, nf);
   _fmpq_poly_set_length(NF_ELEM(a), fmpq_poly_degree(nf->pol));
   _fmpq_poly_normalise(NF_ELEM(a));
   a->exact = 1;
}