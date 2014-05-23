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

    Copyright (C) 2014 William Hart

******************************************************************************/

#ifndef NF_ELEM_APPROX_H
#define NF_ELEM_APPROX_H

#include <mpir.h>
#include "flint.h"
#include "fmpq_poly.h"
#include "nf.h"
#include "nf_elem.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* 
   element in a number field (usual polynomial representation
   nf_elem_t plus conjugates as vector of complex numbers) 
*/
typedef struct 
{
   fmpq_poly_t elem; /* exact repn. (same as nf_elem_t general case) */
   acb_ptr conj; /* appr. conjugates of elem (den of elem is part of repn.) */
   int exact; /* whether exact representation is known for this element */
} nf_elem_approx_struct;

typedef nf_elem_approx_struct nf_elem_approx_t[1];

/*
   NF_ELEM,
   NF_ELEM_NUMREF and
   NF_ELEM_DENREF also work on an nf_elem_approx_t
*/

/******************************************************************************

    Initialisation

******************************************************************************/

void nf_elem_approx_init(nf_elem_approx_t a, const nf_t nf);

void nf_elem_approx_clear(nf_elem_approx_t a, const nf_t nf);

/******************************************************************************

    Conversion

******************************************************************************/

void nf_elem_approx_to_approx(nf_elem_approx_t a, const nf_t nf);

void nf_elem_approx_to_exact(nf_elem_approx_t a, const nf_t nf);

void nf_elem_approx_set_nf_elem(nf_elem_approx_t a, 
                                             const nf_elem_t b, const nf_t nf);

void nf_elem_approx_get_nf_elem(nf_elem_t a,
                                            nf_elem_approx_t b, const nf_t nf);


/******************************************************************************

    Comparison

******************************************************************************/


/******************************************************************************

    I/O

******************************************************************************/


/******************************************************************************

    Arithmetic

******************************************************************************/


void nf_elem_approx_add(nf_elem_approx_t a, nf_elem_approx_t b, 
                                            nf_elem_approx_t c, const nf_t nf);

void nf_elem_approx_sub(nf_elem_approx_t a, nf_elem_approx_t b, 
                                            nf_elem_approx_t c, const nf_t nf);

void nf_elem_approx_mul(nf_elem_approx_t a, nf_elem_approx_t b, 
                                            nf_elem_approx_t c, const nf_t nf);

void nf_elem_approx_inv(nf_elem_approx_t a, nf_elem_approx_t b, const nf_t nf);

void nf_elem_approx_div(nf_elem_approx_t a, nf_elem_approx_t b,
                                            nf_elem_approx_t c, const nf_t nf);

#ifdef __cplusplus
}
#endif

#endif
