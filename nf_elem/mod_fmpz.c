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

        if (nf->flag & NF_LINEAR)
        {
            nf_elem_scalar_div_fmpz(res, res, LNF_ELEM_DENREF(a), nf);
        }
        else if (nf->flag & NF_QUADRATIC)
        {
            nf_elem_scalar_div_fmpz(res, res, QNF_ELEM_DENREF(a), nf);
        }
        else 
        {
            nf_elem_scalar_div_fmpz(res, res, NF_ELEM_DENREF(a), nf);
        }

        fmpz_clear(t);
    }
}

void nf_elem_mod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
{
    nf_elem_mod_fmpz_den(res, a, mod, nf, 1);
}
    
static void
_fmpz_ppio(fmpz_t ppi, fmpz_t ppo, const fmpz_t a, const fmpz_t b)
{
    fmpz_t c, n, g;

    fmpz_init(c);
    fmpz_init(n);
    fmpz_init(g);

    fmpz_gcd(c, a, b);
    fmpz_divexact(n, a, c);
    fmpz_gcd(g, c, n);

    while (!fmpz_is_one(g))
    {
        fmpz_mul(c, c, g);
        fmpz_divexact(n, n, g);
        fmpz_gcd(g, c, n);
    }
    fmpz_set(ppi, c);
    fmpz_set(ppo, n);

    fmpz_clear(c);
    fmpz_clear(n);
    fmpz_clear(g);
}

void
nf_elem_coprime_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
{
    if (nf_elem_is_zero(a, nf))
    {
        nf_elem_zero(res, nf);
        return;
    }

    if (nf_elem_den_is_one(a, nf))
    {
        nf_elem_mod_fmpz_den(res, a, mod, nf, 0);
        return ;
    }

    if (nf->flag & NF_LINEAR)
    {
        fmpz_t c, nc;
        fmpz_init(c);
        fmpz_init(nc);

        _fmpz_ppio(c, nc, LNF_ELEM_DENREF(a), mod);
        fmpz_mul(LNF_ELEM_DENREF(res), mod, c);
        fmpz_invmod(nc, nc, LNF_ELEM_DENREF(res));
        fmpz_mul(LNF_ELEM_NUMREF(res), LNF_ELEM_NUMREF(a), nc);
        fmpz_mod(LNF_ELEM_NUMREF(res), LNF_ELEM_NUMREF(res), LNF_ELEM_DENREF(res));
        fmpz_set(LNF_ELEM_DENREF(res), c);

        fmpz_clear(c);
        fmpz_clear(nc);
    }
    else if (nf->flag & NF_QUADRATIC)
    { 
        fmpz_t c, nc;
        fmpz_init(c);
        fmpz_init(nc);

        _fmpz_ppio(c, nc, QNF_ELEM_DENREF(a), mod);
        fmpz_mul(QNF_ELEM_DENREF(res), mod, c);
        fmpz_invmod(nc, nc, QNF_ELEM_DENREF(res));
        _fmpz_vec_scalar_mul_fmpz(QNF_ELEM_NUMREF(res), QNF_ELEM_NUMREF(a), 3, nc);
        _fmpz_vec_scalar_mod_fmpz(QNF_ELEM_NUMREF(res), QNF_ELEM_NUMREF(res), 3, QNF_ELEM_DENREF(res));
        fmpz_set(QNF_ELEM_DENREF(res), c);
        
        fmpz_clear(c);
        fmpz_clear(nc);
    }
    else
    {
        fmpz_t c, nc;
        fmpz_init(c);
        fmpz_init(nc);

        fmpq_poly_fit_length(NF_ELEM(res), fmpq_poly_length(NF_ELEM(a)));
        _fmpz_ppio(c, nc, NF_ELEM_DENREF(a), mod);
        fmpz_mul(NF_ELEM_DENREF(res), mod, c);
        fmpz_invmod(nc, nc, NF_ELEM_DENREF(res));
        _fmpz_vec_scalar_mul_fmpz(NF_ELEM(res)->coeffs, NF_ELEM(a)->coeffs, fmpq_poly_length(NF_ELEM(a)), nc);
        _fmpz_vec_scalar_mod_fmpz(NF_ELEM(res)->coeffs, NF_ELEM(res)->coeffs, fmpq_poly_length(NF_ELEM(a)), NF_ELEM_DENREF(res));
        fmpz_set(NF_ELEM_DENREF(res), c);
        _fmpq_poly_set_length(NF_ELEM(res), fmpq_poly_length(NF_ELEM(a)));

        fmpz_clear(c);
        fmpz_clear(nc);
    }
    nf_elem_canonicalise(res, nf);
}
