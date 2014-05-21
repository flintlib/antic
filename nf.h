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
#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct {
   fmpq_poly_t pol;  /* defining polynomial */
   union {
      /* insert any precomputed inverse for zz case here */
      fmpz_preinvn_t qq; /* precomputed inverse for leading coeff of num(pol), QQ case */
   } pinv;
   union { /* powers of the generator mod pol */
      fmpq_poly_powers_precomp_t qq; 
      fmpz_poly_powers_precomp_t zz;
   } powers;
   ulong flag;       /* 1 = pol monic over ZZ, 2 = linear, 4 = quadratic field */
   acb_ptr roots;  /* numerical approximations of polynomial roots */
   slong roots_prec; /* precision of roots */
   acb_mat_t V; /* V*[a0, a1, a2]~ = conjugates of alpha if alpha = a0 + a1*x + ... */
   acb_mat_t Vinv; /* inverse of V */
   slong Vprec; /* precision of V */
} nf_struct;

typedef nf_struct nf_t[1];

#define NF_POWERS_CUTOFF 30 /* maximum length of pol where we precompute powers */

#define NF_GENERIC 0
#define NF_MONIC 1
#define NF_LINEAR 2
#define NF_QUADRATIC 4

/******************************************************************************

    Initialisation

******************************************************************************/

void nf_init(nf_t nf, fmpq_poly_t pol);

void nf_clear(nf_t nf);

void nf_print(const nf_t nf);

void _nf_compute_roots(acb_ptr roots, const fmpz_poly_t poly,
    slong initial_prec, slong target_prec);

void nf_compute_roots(nf_t nf, slong prec);

#ifdef __cplusplus
}
#endif

#endif
