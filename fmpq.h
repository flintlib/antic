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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#ifndef FMPQ_H
#define FMPQ_H

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"

typedef struct
{
    fmpz num;
    fmpz den;
}
fmpq;

typedef fmpq fmpq_t[1];


static __inline__ void fmpq_init(fmpq_t x)
{
    x->num = 0L;
    x->den = 1L;
}

static __inline__ void fmpq_clear(fmpq_t x)
{
    fmpz_clear(&x->num);
    fmpz_clear(&x->den);
}

static __inline__ void fmpq_xzero(fmpq_t res)
{
    fmpz_zero(&res->num);
    fmpz_set_ui(&res->den, 1UL);
}

void _fmpq_canonicalise(fmpz_t num, fmpz_t den);

void fmpq_canonicalise(fmpq_t res);


void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, long p, ulong q);

void fmpq_set_si(fmpq_t res, long p, ulong q);

void _fmpq_print(fmpz_t num, fmpz_t den);

void fmpq_print(const fmpq_t x);




void _fmpq_add(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2);


void _fmpq_sub(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2);


void _fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);


void _fmpq_addmul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);


void _fmpq_submul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

void fmpq_submul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);



int _fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod);

int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod);

int _fmpq_reconstruct_fmpz(fmpz_t num, fmpz_t den, const fmpz_t a, const fmpz_t m);

int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m);

#endif