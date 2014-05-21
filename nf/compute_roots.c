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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "nf.h"
#include "acb_poly.h"

#define DEBUG 0

static __inline__ int
check_accuracy(acb_srcptr vec, slong len, slong prec)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (mag_cmp_2exp_si(arb_radref(acb_realref(vec + i)), -prec) >= 0
         || mag_cmp_2exp_si(arb_radref(acb_imagref(vec + i)), -prec) >= 0)
            return 0;
    }

    return 1;
}

void
_nf_compute_roots(acb_ptr roots, const fmpz_poly_t poly,
    slong initial_prec, slong target_prec)
{
    slong prec, deg, isolated, maxiter;
    acb_poly_t cpoly;

    deg = poly->length - 1;

    acb_poly_init(cpoly);
    roots = _acb_vec_init(deg);

    for (prec = initial_prec; ; prec *= 2)
    {
        acb_poly_set_fmpz_poly(cpoly, poly, prec);
        maxiter = FLINT_MIN(FLINT_MAX(deg, 32), prec);

        isolated = acb_poly_find_roots(roots, cpoly,
            prec == initial_prec ? NULL : roots, maxiter, prec);

        if (isolated == deg && check_accuracy(roots, deg, target_prec))
            break;
    }

#if DEBUG
    {
        slong i;

        printf("roots:\n");
        for (i = 0; i < deg; i++)
        {
            acb_printd(roots + i, 30);
            printf("\n");
        }
    }
#endif

    acb_poly_clear(cpoly);
    _acb_vec_clear(roots, deg);
}

void
nf_compute_roots(nf_t nf, slong prec)
{
    fmpz_poly_t t;
    fmpz_poly_init(t);

    fmpq_poly_get_numerator(t, nf->pol);

    if (!fmpz_poly_is_squarefree(t))
    {
        printf("nf_compute_roots: polynomial must be squarefree!\n");
    }
    else
    {
        if (nf->roots == NULL)
            nf->roots = _acb_vec_init(fmpz_poly_degree(t));

        _nf_compute_roots(nf->roots, t, 32, prec);
        nf->roots_prec = prec;
    }

    fmpz_poly_clear(t);
}

