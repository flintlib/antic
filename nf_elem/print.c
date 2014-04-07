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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "nf_elem.h"

void nf_elem_print(const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        flint_printf("(");
        fmpz_print(LNF_ELEM_NUMREF(a));
        flint_printf(" ");
        fmpz_print(LNF_ELEM_DENREF(a));
        flint_printf(")");        
    } else if (nf->flag & NF_QUADRATIC)
    {
        const fmpz * const anum = QNF_ELEM_NUMREF(a);
        const fmpz * const aden = QNF_ELEM_DENREF(a);
        
        flint_printf("(");
        fmpz_print(anum + 1);
        flint_printf(" ");
        fmpz_print(anum);
        flint_printf(" ");
        fmpz_print(aden);
        flint_printf(")");
    } else if (nf->flag & NF_MONIC)
    {
       _fmpz_poly_fprint_pretty(stdout, NF_ELEM_NUMREF(a), NF_ELEM(a)->length, "x");
    } else
    {
        fmpq_poly_print_pretty(NF_ELEM(a), "x");
    }
}

