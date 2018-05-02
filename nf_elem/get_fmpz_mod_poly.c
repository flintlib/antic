/*=============================================================================

    This file is part of ANTIC.

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

    Copyright (C) 2018 Tommy Hofmann

******************************************************************************/

#include "nf_elem.h"

void _nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf)
{
    if (nf_elem_is_zero(a, nf))
    {
        fmpz_mod_poly_zero(pol);
        return;
    }
    if (nf->flag & NF_LINEAR)
    {
        {
            fmpz_mod_poly_fit_length(pol, 1);
            fmpz_mod(pol->coeffs + 0, LNF_ELEM_NUMREF(a), &(pol->p));
            _fmpz_mod_poly_set_length(pol, 1);
            _fmpz_mod_poly_normalise(pol);

        }
    } else if (nf->flag & NF_QUADRATIC)
    {
        fmpz_mod_poly_fit_length(pol, 3);
        fmpz_mod(pol->coeffs + 0, QNF_ELEM_NUMREF(a), &(pol->p));
        fmpz_mod(pol->coeffs + 1, QNF_ELEM_NUMREF(a) + 1, &(pol->p));
        fmpz_mod(pol->coeffs + 2, QNF_ELEM_NUMREF(a) + 2, &(pol->p));
        _fmpz_mod_poly_set_length(pol, 3);
        _fmpz_mod_poly_normalise(pol);
    } else
    {
        slong len = NF_ELEM(a)->length;
        slong i;
        fmpz_mod_poly_fit_length(pol, len);
        for (i = 0; i < len; i++)
            fmpz_mod(pol->coeffs + i, NF_ELEM_NUMREF(a) + i, &(pol->p));
        _fmpz_mod_poly_set_length(pol, len);
        _fmpz_mod_poly_normalise(pol);
    }
}

void nf_elem_get_fmpz_mod_poly_den(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf, int den)
{
    _nf_elem_get_fmpz_mod_poly(pol, a, nf);
    if (den)
    {
        if (nf->flag & NF_LINEAR)
            fmpz_mod_poly_scalar_div_fmpz(pol, pol, LNF_ELEM_DENREF(a));
        else if (nf->flag & NF_QUADRATIC)
            fmpz_mod_poly_scalar_div_fmpz(pol, pol, QNF_ELEM_DENREF(a));
        else
            fmpz_mod_poly_scalar_div_fmpz(pol, pol, NF_ELEM_DENREF(a));
    }
}

void nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf)
{
    nf_elem_get_fmpz_mod_poly_den(pol, a, nf, 1);
}
