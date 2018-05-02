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

void
_nf_elem_mod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
{
    if (nf_elem_is_zero(a, nf))
    {
        nf_elem_zero(res, nf);
        return;
    }
    if (nf->flag & NF_LINEAR)
    {
        fmpz_mod(LNF_ELEM_NUMREF(res), LNF_ELEM_NUMREF(a), mod);
        fmpz_one(LNF_ELEM_DENREF(res));
    }
    else if (nf->flag & NF_QUADRATIC)
    { 
        _fmpz_vec_scalar_mod_fmpz(QNF_ELEM_NUMREF(res), QNF_ELEM_NUMREF(a), 3, mod);
        fmpz_one(QNF_ELEM_DENREF(res));
    }
    else
    {
        fmpq_poly_fit_length(NF_ELEM(res), fmpq_poly_length(NF_ELEM(a)));
        _fmpq_poly_set_length(NF_ELEM(res), fmpq_poly_length(NF_ELEM(a)));
        _fmpz_vec_scalar_mod_fmpz(NF_ELEM(res)->coeffs, NF_ELEM(a)->coeffs, fmpq_poly_length(NF_ELEM(a)), mod);
        fmpz_one(NF_ELEM_DENREF(res));
    }
    nf_elem_canonicalise(res, nf);
}

void
nf_elem_mod_fmpz_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int den)
{
    if (!den || nf_elem_den_is_one(a, nf))
    {
        _nf_elem_mod_fmpz(res, a, mod, nf);
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);

        nf_elem_get_den(t, a, nf);
        fmpz_mul(t, t, mod);

        _nf_elem_mod_fmpz(res, a, t, nf);

        fmpz_clear(t);
    }
}

void nf_elem_mod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
{
    nf_elem_mod_fmpz_den(res, a, mod, nf, 1);
}
