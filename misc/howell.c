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

    Copyright (C) 2015 Tommy Hofmann 

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "fmpz_mat.h"

/* This is copied from nmod_mat/lu_classic.c */
static __inline__ int
nmod_mat_pivot(nmod_mat_t A, slong start_row, slong col)
{
    slong j;
    mp_ptr u;

    if (nmod_mat_entry(A, start_row, col) != 0)
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (nmod_mat_entry(A, j, col) != 0)
        {
            u = A->rows[j];
            A->rows[j] = A->rows[start_row];
            A->rows[start_row] = u;

            return -1;
        }
    }
    return 0;
}

static __inline__ int
fmpz_mat_pivot(fmpz_mat_t A, slong start_row, slong col)
{
    slong j;
    fmpz * u;

    if (fmpz_is_zero(fmpz_mat_entry(A, start_row, col)) == 0)
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (nmod_mat_entry(A, j, col) != 0)
        {
            u = A->rows[j];
            A->rows[j] = A->rows[start_row];
            A->rows[start_row] = u;

            return -1;
        }
    }
    return 0;
}


void
_fmpz_mat_swap_rows(fmpz_mat_t A, slong i, slong j)
{
    fmpz * u;
    u = A->rows[i];
    A->rows[i] = A->rows[j];
    A->rows[j] = u;
}

void
_nmod_mat_swap_rows(nmod_mat_t A, slong i, slong j)
{
    mp_ptr u;
    u = A->rows[i];
    A->rows[i] = A->rows[j];
    A->rows[j] = u;
}

void
n_ppio(mp_ptr ppi, mp_ptr ppo, mp_limb_t a, mp_limb_t b)
{
    mp_limb_t c, n, m, g;

    c = n_gcd(a, b);
    n = a/c;
    m = c;
    g = n_gcd(c, n);
    while( g != 1 )
    {
        c = c * g;
        n = n/g;
        g = n_gcd(c, n);
    }
    *ppi = c;
    *ppo = n;
}

void
fmpz_ppio(fmpz_t ppi, fmpz_t ppo, fmpz_t a, fmpz_t b)
{
    fmpz_t c, n, g;

    fmpz_init(c);
    fmpz_init(n);
    fmpz_init(g);

    fmpz_gcd(c, a, b);
    fmpz_divexact(n, a, c);
    fmpz_gcd(g, c, n);

    while ( fmpz_is_one(g) != 1 )
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

mp_limb_t
n_stab(mp_limb_t a, mp_limb_t b, nmod_t N)
{
    mp_limb_t g, s, t;
    g = n_gcd(a, b);
    b = n_gcd(g, N.n);
    n_ppio(&s, &t, N.n/b, a/b);
    return t;
}

void
fmpz_stab(fmpz_t t, fmpz_t a, fmpz_t b, fmpz_t N)
{
    fmpz_t g, gg, s, aa, bb;

    fmpz_init(g);
    fmpz_init(gg);
    fmpz_init(s);
    fmpz_init(aa);
    fmpz_init(bb);

    fmpz_gcd(g, a, b);
    fmpz_gcd(gg, g, N);

    fmpz_divexact(bb, N, b);
    fmpz_divexact(aa, a, gg);

    fmpz_ppio(s, t, bb, aa);

    fmpz_clear(g);
    fmpz_clear(gg);
    fmpz_clear(s);
    fmpz_clear(aa);
    fmpz_clear(bb);
}

mp_limb_t
n_unit(mp_limb_t a, nmod_t N)
{
    mp_limb_t g, s, t, l, d;

    g = n_xgcd(&s, &t, a, N.n);

    if (g == 1)
    {
        return s;
    }
    else
    {
        l = N.n/g;
        d = n_stab(s, l, N);
        return nmod_add(s, nmod_mul(d, l, N), N);
    }
}

void
fmpz_unit(fmpz_t u, fmpz_t a, fmpz_t N)
{
    fmpz_t g, s, t, l, d, k;

    fmpz_init(g);
    fmpz_init(s);
    fmpz_init(t);
    fmpz_init(l);
    fmpz_init(d);
    fmpz_init(k);

    fmpz_xgcd(g, s, t, a, N);

    if (fmpz_is_one(g) == 1)
    {
        fmpz_set(u, s);
    }
    else
    {
        fmpz_divexact(l, N, g);
        fmpz_stab(d, s, l, N);
        fmpz_mul(u, d, l);
        fmpz_add(u, u, s);
        fmpz_mod(u, u, N);
    }
}

mp_limb_t
_n_unit(mp_limb_t a, mp_limb_t N)
{
    nmod_t mod;

    nmod_init(&mod, N);

    return n_unit(a, mod);
}

/* Number of rows has to be larger then the number of columns */
slong
_nmod_mat_howell(nmod_mat_t A)
{
    mp_limb_t s, t, u, v, q, t1, t2, g;
    slong m, n, row, col, i, k, l;
    mp_limb_t **r;
    nmod_t mod;
    mp_ptr extra_row;


    n = A->r;
    m = A->c;
    r = A->rows;
    mod = A->mod;

    extra_row = _nmod_vec_init(m);

    row = col = 0;

    while (row < n && col < m)
    {
        if (nmod_mat_pivot(A, row, col) == 0)
        {
            col++;
            continue;
        }
        for (i = row + 1; i < n; i++)
        {
            g = n_xgcd(&s, &t, nmod_mat_entry(A, row, col), nmod_mat_entry(A, i, col));
            /* now g = a*x - b*y 
             a,b < x < mod.n */
            t = nmod_neg(t, mod);
            u = (nmod_mat_entry(A, i, col))/g;
            u = nmod_neg(u, mod);
            v = (nmod_mat_entry(A, row, col))/g;
            /* now g = a*x + b*y and 0 = sv - tu = 1 modulo mod.n */
            
            for (k = col; k < m; k++)
            {
                t1 = nmod_add(nmod_mul(s, nmod_mat_entry(A, row, k), mod), nmod_mul(t, nmod_mat_entry(A, i, k), mod), mod);
                t2 = nmod_add(nmod_mul(u, nmod_mat_entry(A, row, k), mod), nmod_mul(v, nmod_mat_entry(A, i, k), mod), mod);
                nmod_mat_entry(A, row, k) = t1;
                nmod_mat_entry(A, i, k) = t2;
            }
        }
        row++;
        col++;
    }

    
    for (col = 0; col < m; col++)
    {
        if (nmod_mat_entry(A, col, col) != 0)
        {
            u = n_unit(nmod_mat_entry(A, col, col), mod);
            for (k = col; k < m; k++)
            {
                nmod_mat_entry(A, col, k) = nmod_mul(u, nmod_mat_entry(A, col, k), mod);
            }
            for (row = 0; row < col ; row++)
            {
                
                
                q = nmod_mat_entry(A, row, col)/nmod_mat_entry(A, col, col);

                for (l = row; l< m; l++)
                {
                    s = nmod_sub(nmod_mat_entry(A, row, l), nmod_mul(q, nmod_mat_entry(A, col, l), mod), mod);
                    nmod_mat_entry(A, row, l) = s;
                }
            }
            g = n_gcd(mod.n, nmod_mat_entry(A, col, col));
            if (g == 1)
            {
                continue;

            }
            g = mod.n/g;
            _nmod_vec_scalar_mul_nmod(extra_row, r[col], m, g, mod);
        }
        else
        {
          _nmod_vec_set(extra_row, r[col], m);
        }

        for (row = col + 1; row < m; row++)
        {
            g = n_xgcd(&s, &t, nmod_mat_entry(A, row, row), extra_row[row]);
            if (g == 0)
            {
                continue;
            }
            t = nmod_neg(t, mod);
            u = extra_row[row]/g;
            u = nmod_neg(u, mod);
            v = (nmod_mat_entry(A, row, row))/g;
            /* now g = a*x + b*y and 0 = sv - tu = 1 modulo mod.n */
            
            for (k = row; k < m; k++)
            {
                t1 = nmod_add(nmod_mul(s, nmod_mat_entry(A, row, k), mod), nmod_mul(t, extra_row[k], mod), mod);
                t2 = nmod_add(nmod_mul(u, nmod_mat_entry(A, row, k), mod), nmod_mul(v, extra_row[k], mod), mod);
                nmod_mat_entry(A, row, k) = t1;
                extra_row[k] = t2;
            }
        }
    }
    _nmod_vec_clear(extra_row);

    return 0;
}

slong
_fmpz_mat_howell(fmpz_mat_t A, fmpz_t mod)
{
    fmpz_t s, t, q, u, v, t1, t2, g;
    slong m, n, row, col, i, k, l;
    fmpz  ** r;
    fmpz * extra_row;

    n = A->r;
    m = A->c;
    r = A->rows;

    extra_row = _fmpz_vec_init(m);



    row = col = 0;

    for(row = 0; row < n; row++)
    {
        for(col = 0; col < m; col++)
        {
            fmpz_mod(fmpz_mat_entry(A, row, col), fmpz_mat_entry(A, row, col), mod);
        }
    }

    row = col = 0;

    while (row < n && col < m)
    {
        if (fmpz_mat_pivot(A, row, col) == 0)
        {
            col++;
            continue;
        }
        for (i = row + 1; i < n; i++)
        {
            fmpz_xgcd(g, s, t, fmpz_mat_entry(A, row, col), fmpz_mat_entry(A, i, col));
            fmpz_divexact(u, fmpz_mat_entry(A, i, col), g);
            fmpz_neg(u, u);
            fmpz_divexact(v, fmpz_mat_entry(A, row, col), g);
            
            for (k = col; k < m; k++)
            {
                fmpz_mul(t1, s, fmpz_mat_entry(A, row, k));
                fmpz_addmul(t1, t, fmpz_mat_entry(A, i, k));
                fmpz_mod(t1, t1, mod);
                fmpz_mul(t2, u, fmpz_mat_entry(A, row, k));
                fmpz_addmul(t2, v, fmpz_mat_entry(A, i, k));
                fmpz_mod(t2, t2, mod);
                fmpz_set(fmpz_mat_entry(A, row, k), t1);
                fmpz_set(fmpz_mat_entry(A, i, k), t2);
            }
        }
        row++;
        col++;
    }
    
    for (col = 0; col < m; col++)
    {
        if (fmpz_is_zero(fmpz_mat_entry(A, col, col)) != 1)
        {
            fmpz_unit(u, fmpz_mat_entry(A, col, col), mod);
            for (k = col ; k < m; k++)
            {
                fmpz_mul(fmpz_mat_entry(A, col, k), u, fmpz_mat_entry(A, col, k));
                fmpz_mod(fmpz_mat_entry(A, col, k), fmpz_mat_entry(A, col, k), mod);
            }
            for (row = 0; row < col ; row++)
            {
                
                
                fmpz_divexact(q, fmpz_mat_entry(A, row, col), fmpz_mat_entry(A, col, col));

                for (l = row; l< m; l++)
                {
                    fmpz_submul(fmpz_mat_entry(A, row, l), q, fmpz_mat_entry(A, col, l));
                    fmpz_mod(fmpz_mat_entry(A, row, l), fmpz_mat_entry(A, row, l), mod);
                }
            }
            fmpz_gcd(g, mod, fmpz_mat_entry(A, col, col));

            if (fmpz_is_one(g) == 1)
            {
                continue;

            }
            fmpz_divexact(g, mod, g);
            _fmpz_vec_scalar_mul_fmpz(extra_row, r[col], m, g);
            _fmpz_vec_scalar_mod_fmpz(extra_row, extra_row, m, mod);
        }
        else
        {
          _fmpz_vec_set(extra_row, r[col], m);
        }

        for (row = col + 1; row < m; row++)
        {
            fmpz_xgcd(g, s, t, fmpz_mat_entry(A, row, row), extra_row + row);
            if (fmpz_is_zero(g))
            {
                continue;
            }
            fmpz_divexact(u, extra_row + row, g);
            fmpz_neg(u, u);
            fmpz_divexact(v, fmpz_mat_entry(A, row, row), g);
            
            for (k = row; k < m; k++)
            {
                fmpz_mul(t1, s, fmpz_mat_entry(A, row, k));
                fmpz_addmul(t1, t, extra_row + k);
                fmpz_mod(t1, t1, mod);
                fmpz_mul(t2, u, fmpz_mat_entry(A, row, k));
                fmpz_addmul(t2, v, extra_row + k);
                fmpz_mod(t2, t2, mod);
                fmpz_set(fmpz_mat_entry(A, row, k), t1);
                fmpz_set(extra_row + k, t2);
            }
        }
    }
    _fmpz_vec_clear(extra_row, m);

    return 0;
}
