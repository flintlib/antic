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

#ifndef NF_ELEM_H
#define NF_ELEM_H

#include <mpir.h>
#include "flint.h"
#include "fmpq_poly.h"
#include "nf.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct /* element of a linear number field */
{
   fmpz_t num;
   fmpz_t den;
} lnf_elem_struct;

typedef lnf_elem_struct lnf_elem_t[1];

typedef struct /* element of a quadratic number field */
{
   fmpz num[3];
   fmpz_t den;
} qnf_elem_struct;

typedef qnf_elem_struct qnf_elem_t[1];

typedef union /* element in a number field (specified by an nf_t) */
{
   fmpq_poly_t elem; /* general case */
   lnf_elem_t lelem; /* linear number field */
   qnf_elem_t qelem; /* quadratic number field */
} nf_elem_struct;

typedef nf_elem_struct nf_elem_t[1];

#define NF_ELEM_NUMREF(xxx) fmpq_poly_numref((xxx)->elem)
#define NF_ELEM_DENREF(xxx) fmpq_poly_denref((xxx)->elem)
#define LNF_ELEM_NUMREF(xxx) ((xxx)->lelem->num)
#define LNF_ELEM_DENREF(xxx) ((xxx)->lelem->den)
#define QNF_ELEM_NUMREF(xxx) ((xxx)->qelem->num)
#define QNF_ELEM_DENREF(xxx) ((xxx)->qelem->den)
#define NF_ELEM(xxx) (xxx)->elem
#define LNF_ELEM(xxx) (xxx)->lelem
#define QNF_ELEM(xxx) (xxx)->qelem

/******************************************************************************

    Initialisation

******************************************************************************/

void nf_elem_init(nf_elem_t a, const nf_t nf);

void nf_elem_clear(nf_elem_t a, const nf_t nf);

void nf_elem_randtest(nf_elem_t a, flint_rand_t state, 
                                              mp_bitcnt_t bits, const nf_t nf);

void nf_elem_randtest_not_zero(nf_elem_t a, flint_rand_t state, 
                                              mp_bitcnt_t bits, const nf_t nf);

static __inline__
void nf_elem_canonicalise(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR) {
      _fmpq_canonicalise(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
      _fmpq_poly_canonicalise(QNF_ELEM_NUMREF(a), QNF_ELEM_DENREF(a), 2);
   else
      fmpq_poly_canonicalise(NF_ELEM(a));
}

int _nf_elem_invertible_check(nf_elem_t a, const nf_t nf);

/******************************************************************************

    Comparison

******************************************************************************/

int _nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf);

int nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf);

static __inline__
int nf_elem_is_one(const nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      return fmpz_is_one(LNF_ELEM_NUMREF(a)) && fmpz_is_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      return fmpz_is_one(anum) && fmpz_is_zero(anum + 1) 
          && fmpz_is_one(QNF_ELEM_DENREF(a));
   } else
      return fmpq_poly_is_one(a->elem);
}

/******************************************************************************

    I/O

******************************************************************************/

void nf_elem_print(const nf_elem_t a, const nf_t nf);

/******************************************************************************

    Arithmetic

******************************************************************************/

static __inline__
void nf_elem_zero(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_zero(LNF_ELEM_NUMREF(a));
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      fmpz_zero(anum);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_zero(NF_ELEM(a));
}

static __inline__
void nf_elem_one(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_one(LNF_ELEM_NUMREF(a));
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      fmpz_one(anum);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_one(NF_ELEM(a));
}

static __inline__
void nf_elem_set(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b));
      fmpz_set(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      
      fmpz_set(anum, bnum);
      fmpz_set(anum + 1, bnum + 1);
      fmpz_set(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));
   } else
      fmpq_poly_set(NF_ELEM(a), NF_ELEM(b));
}

static __inline__
void nf_elem_neg(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_neg(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b));
      fmpz_set(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      
      fmpz_neg(anum, bnum);
      fmpz_neg(anum + 1, bnum + 1);
      fmpz_set(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));
   } else
      fmpq_poly_neg(NF_ELEM(a), NF_ELEM(b));
}

static __inline__
void nf_elem_swap(nf_elem_t a, nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_swap(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b));
      fmpz_swap(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      fmpz * const bnum = QNF_ELEM_NUMREF(b);
      
      fmpz_swap(anum, bnum);
      fmpz_swap(anum + 1, bnum + 1);
      fmpz_swap(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));
   } else
      fmpq_poly_swap(NF_ELEM(a), NF_ELEM(b));
}

void _nf_elem_add_lf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can);

void _nf_elem_sub_lf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can);

void _nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can);

void _nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can);

void nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, 
                                            const nf_elem_t c, const nf_t nf);

void nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b, 
                                            const nf_elem_t c, const nf_t nf);

static __inline__
void _nf_elem_add(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
      _nf_elem_add_lf(a, b, c, nf, 0);
   else if (nf->flag & NF_QUADRATIC)
      _nf_elem_add_qf(a, b, c, nf, 0);
   else
      fmpq_poly_add_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 0);
}

static __inline__
void _nf_elem_sub(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
      _nf_elem_sub_lf(a, b, c, nf, 0);
   else if (nf->flag & NF_QUADRATIC)
      _nf_elem_sub_qf(a, b, c, nf, 0);
   else
      fmpq_poly_sub_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 0);
}

static __inline__
void nf_elem_add(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
      _nf_elem_add_lf(a, b, c, nf, 1);
   else if (nf->flag & NF_QUADRATIC)
      nf_elem_add_qf(a, b, c, nf);
   else
      fmpq_poly_add_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 1);
}

static __inline__
void nf_elem_sub(nf_elem_t a, const nf_elem_t b, 
                                              const nf_elem_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
      _nf_elem_sub_lf(a, b, c, nf, 1);
   else if (nf->flag & NF_QUADRATIC)
      nf_elem_sub_qf(a, b, c, nf);
   else
      fmpq_poly_sub_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 1);
}

void _nf_elem_mul(nf_elem_t a, const nf_elem_t b, 
                                             const nf_elem_t c, const nf_t nf);

void nf_elem_mul(nf_elem_t a, const nf_elem_t b, 
                                             const nf_elem_t c, const nf_t nf);

void _nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, 
                                    const nf_elem_t c, const nf_t nf, int red);

void nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, 
                                    const nf_elem_t c, const nf_t nf, int red);

void _nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf);

void nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf);

void _nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf);

void nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf);

#ifdef __cplusplus
}
#endif

#endif
