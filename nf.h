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

    Copyright (C) 2013 William Hart

******************************************************************************/

#ifndef NF_H
#define NF_H

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct {
   fmpq_poly_t pol;  /* defining polynomial */
   union {
      fmpz_poly_t zz; /* precomputed inverse for reduction mod pol ZZ case */
      fmpz_preinvn_t qq; /* precomputed inverse for leading coeff of num(pol), QQ case */
   } pinv;
   fmpq_poly_struct * powers; /* powers of the generator mod pol */
   ulong flag;       /* 1 = pol monic over ZZ, 2 = quadratic field */
} nf_struct;

typedef nf_struct nf_t[1];

#define NF_POWERS_CUTOFF 12 /* maximum length of pol where we precompute powers */

#define NF_GENERIC 0
#define NF_MONIC 1
#define NF_QUADRATIC 2

/******************************************************************************

    Initialisation

******************************************************************************/

void nf_init(nf_t nf, fmpq_poly_t pol);

void nf_clear(nf_t nf);

void nf_print(const nf_t nf);

#ifdef __cplusplus
}
#endif

#endif
