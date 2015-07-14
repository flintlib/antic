/*=============================================================================

    This file is part of ANTIC

    ANTIC is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ANTIC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Claus Fieker

******************************************************************************/

#ifndef ANTIC_MISC_H
#define ANTIC_MISC_H

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <fmpz_poly.h>
#include <fmpq_poly.h>

#ifdef __cplusplus
 extern "C" {
#endif

FLINT_DLL void _fmpz_poly_resultant_modular_div(fmpz_t res, const fmpz * poly1, slong len1, 
       const fmpz * poly2, slong len2, const fmpz_t divisor, slong num_primes);

FLINT_DLL void fmpz_poly_resultant_modular_div(fmpz_t res, const fmpz_poly_t poly1,
                                                      const fmpz_poly_t poly2,
                                      const fmpz_t divisor, slong num_primes);


FLINT_DLL void _fmpq_poly_resultant_div(fmpz_t rnum, fmpz_t rden, 
                          const fmpz *poly1, const fmpz_t den1, slong len1, 
                          const fmpz *poly2, const fmpz_t den2, slong len2,
                          const fmpz_t divisor, slong num_primes);

FLINT_DLL void fmpq_poly_resultant_div(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g, const fmpz_t divisor, slong num_primes);

#ifdef __cplusplus
}
#endif
#endif

