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

#ifdef NF_ELEM_INLINES_C
#define NF_ELEM_INLINE FLINT_DLL
#else
#define NF_ELEM_INLINE static __inline__
#endif

#include "gmp.h"
#include "flint.h"
#include "fmpq_poly.h"
#include "fmpz_mat.h"
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

FLINT_DLL void nf_elem_init(nf_elem_t a, const nf_t nf);

FLINT_DLL void nf_elem_clear(nf_elem_t a, const nf_t nf);

FLINT_DLL void nf_elem_randtest(nf_elem_t a, flint_rand_t state, 
                                              mp_bitcnt_t bits, const nf_t nf);

FLINT_DLL void nf_elem_randtest_not_zero(nf_elem_t a, flint_rand_t state, 
                                              mp_bitcnt_t bits, const nf_t nf);

NF_ELEM_INLINE
void nf_elem_canonicalise(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR) {
      _fmpq_canonicalise(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
      _fmpq_poly_canonicalise(QNF_ELEM_NUMREF(a), QNF_ELEM_DENREF(a), 3);
   else
      fmpq_poly_canonicalise(NF_ELEM(a));
}

FLINT_DLL void _nf_elem_reduce(nf_elem_t a, const nf_t nf);

FLINT_DLL void nf_elem_reduce(nf_elem_t a, const nf_t nf);

FLINT_DLL int _nf_elem_invertible_check(nf_elem_t a, const nf_t nf);

/******************************************************************************

    Comparison

******************************************************************************/

FLINT_DLL int _nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf);

FLINT_DLL int nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf);

NF_ELEM_INLINE
int nf_elem_is_zero(const nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      return fmpz_is_zero(LNF_ELEM_NUMREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      return fmpz_is_zero(anum) && fmpz_is_zero(anum + 1);
   } else
      return fmpq_poly_is_zero(a->elem);
}

NF_ELEM_INLINE
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

NF_ELEM_INLINE
int nf_elem_is_gen(const nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_t t1, t2;
	  int is_gen;
	  
	  /* fast path */
	  if (fmpz_equal(LNF_ELEM_DENREF(a), nf->pol->coeffs + 1))
	    return fmpz_cmpabs(LNF_ELEM_DENREF(a), nf->pol->coeffs) == 0
            && fmpz_sgn(LNF_ELEM_DENREF(a)) == -fmpz_sgn(nf->pol->coeffs);	
			
	  /* slow path */
	  fmpz_init(t1);
	  fmpz_init(t2);
	  
	  fmpz_mul(t1, LNF_ELEM_NUMREF(a), nf->pol->coeffs + 1);
	  fmpz_mul(t2, LNF_ELEM_DENREF(a), nf->pol->coeffs);
	  fmpz_neg(t1, t1);
	  
	  is_gen = fmpz_equal(t1, t2);
	  
	  fmpz_clear(t1);
	  fmpz_clear(t2);
	  
	  return is_gen;
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      return fmpz_equal(anum + 1, QNF_ELEM_DENREF(a)) 
	      && fmpz_is_zero(anum);
   } else
      return fmpq_poly_length(NF_ELEM(a)) == 2
	      && fmpz_equal(NF_ELEM(a)->coeffs + 1, NF_ELEM(a)->den)
		  && fmpz_is_zero(NF_ELEM(a)->coeffs);
}

/******************************************************************************

    I/O

******************************************************************************/

FLINT_DLL void nf_elem_print_pretty(const nf_elem_t a, 
                             const nf_t nf, const char * var);

NF_ELEM_INLINE
char * nf_elem_get_str_pretty(const nf_elem_t a, 
                              const char * var, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      const fmpz * const den = LNF_ELEM_DENREF(a);
      const fmpz * const num = LNF_ELEM_NUMREF(a);
      slong len = 1 - fmpz_is_zero(num);

      return _fmpq_poly_get_str_pretty(num, den, len, var);
   }
   else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const den = QNF_ELEM_DENREF(a);
      const fmpz * const num = QNF_ELEM_NUMREF(a);
      slong len = 3;

      while (len != 0 && fmpz_is_zero(num + len - 1))
         len--;

      return _fmpq_poly_get_str_pretty(num, den, len, var);
   } else
   {
      return fmpq_poly_get_str_pretty(NF_ELEM(a), var);
   }
}

/******************************************************************************

    Element creation

******************************************************************************/

NF_ELEM_INLINE
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

NF_ELEM_INLINE
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

NF_ELEM_INLINE
void nf_elem_gen(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_neg(LNF_ELEM_NUMREF(a), nf->pol->coeffs);
      fmpz_set(LNF_ELEM_DENREF(a), nf->pol->coeffs + 1);
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      fmpz_one(anum + 1);
      fmpz_zero(anum);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_set_coeff_ui(NF_ELEM(a), 1, 1);
}

NF_ELEM_INLINE
void nf_elem_set_si(nf_elem_t a, slong c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set_si(LNF_ELEM_NUMREF(a), c);
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      fmpz_set_si(anum, c);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_set_si(NF_ELEM(a), c);
}

NF_ELEM_INLINE
void nf_elem_set_fmpz(nf_elem_t a, const fmpz_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(LNF_ELEM_NUMREF(a), c);
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      fmpz_set(anum, c);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_set_fmpz(NF_ELEM(a), c);
}

NF_ELEM_INLINE
void nf_elem_set_fmpq(nf_elem_t a, const fmpq_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(LNF_ELEM_NUMREF(a), fmpq_numref(c));
      fmpz_set(LNF_ELEM_DENREF(a), fmpq_denref(c));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      fmpz_set(anum, fmpq_numref(c));
      fmpz_zero(anum + 1);
      fmpz_set(QNF_ELEM_DENREF(a), fmpq_denref(c));
   } else
      fmpq_poly_set_fmpq(NF_ELEM(a), c);
}

NF_ELEM_INLINE
void nf_elem_set_fmpq_poly(nf_elem_t a, const fmpq_poly_t pol, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      if (pol->length == 0)
	  {
	     fmpz_zero(LNF_ELEM_NUMREF(a));
		 fmpz_one(LNF_ELEM_DENREF(a));
	  } else
	  {
	     fmpz_set(LNF_ELEM_NUMREF(a), fmpq_poly_numref(pol));
         fmpz_set(LNF_ELEM_DENREF(a), fmpq_poly_denref(pol));
	  }
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      
      if (pol->length == 0)
	  {
	     fmpz_zero(anum);
		 fmpz_zero(anum + 1);
		 fmpz_one(QNF_ELEM_DENREF(a));
	  } else if (pol->length == 1)
	  {
		 fmpz_zero(anum + 1);
		 fmpz_set(anum, fmpq_poly_numref(pol));
		 fmpz_set(QNF_ELEM_DENREF(a), fmpq_poly_denref(pol));
	  } else
	  {
	     fmpz_set(anum, fmpq_poly_numref(pol));
	     fmpz_set(anum + 1, fmpq_poly_numref(pol) + 1);
         fmpz_set(QNF_ELEM_DENREF(a), fmpq_poly_denref(pol));
	  }
   } else
      fmpq_poly_set(NF_ELEM(a), pol);
}

/******************************************************************************

    Conversion

******************************************************************************/

FLINT_DLL
void nf_elem_set_fmpz_mat_row(nf_elem_t b, const fmpz_mat_t M, 
                                     const slong i, fmpz_t den, const nf_t nf);

FLINT_DLL 
void nf_elem_get_fmpz_mat_row(fmpz_mat_t M, const slong i, fmpz_t den, 
                                             const nf_elem_t b, const nf_t nf);

NF_ELEM_INLINE
void nf_elem_get_fmpq_poly(fmpq_poly_t pol, const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpq_poly_set_fmpz(pol, LNF_ELEM_NUMREF(a));
        fmpz_set(fmpq_poly_denref(pol), LNF_ELEM_DENREF(a));
    } else if (nf->flag & NF_QUADRATIC)
    {
        const fmpz * const anum = QNF_ELEM_NUMREF(a);

        fmpq_poly_fit_length(pol, 2);
        _fmpq_poly_set_length(pol, 2);
        _fmpz_vec_set(pol->coeffs, anum, 2);
        _fmpq_poly_normalise(pol);
        fmpz_set(pol->den, QNF_ELEM_DENREF(a));
    } else
    {
        fmpq_poly_set(pol, NF_ELEM(a));
    }
}

/******************************************************************************
 
    Basic manipulation 

******************************************************************************/

NF_ELEM_INLINE
void nf_elem_get_den(fmpz_t d, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
     fmpz_set(d, LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
     fmpz_set(d, QNF_ELEM_DENREF(b));
   } else
   {
     fmpz_set(d, NF_ELEM_DENREF(b));
   }
}

NF_ELEM_INLINE
void nf_elem_set_den(nf_elem_t b, fmpz_t d, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
     fmpz_set(LNF_ELEM_DENREF(b), d);
   } else if (nf->flag & NF_QUADRATIC)
   {
     fmpz_set(QNF_ELEM_DENREF(b), d);
   } else
   {
     fmpz_set(NF_ELEM_DENREF(b), d);
   }
}

NF_ELEM_INLINE
void nf_elem_get_coeff_fmpq(fmpq_t a, const nf_elem_t b, 
                                                        slong i, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(fmpq_numref(a), LNF_ELEM_NUMREF(b));
      fmpz_set(fmpq_denref(a), LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      
      fmpz_set(fmpq_numref(a), bnum + i);
      fmpz_set(fmpq_denref(a), QNF_ELEM_DENREF(b));

      fmpq_canonicalise(a);
   } else
      fmpq_poly_get_coeff_fmpq(a, NF_ELEM(b), i);
}

NF_ELEM_INLINE
void nf_elem_get_coeff_fmpz(fmpz_t a, const nf_elem_t b, 
                                                        slong i, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(a, LNF_ELEM_NUMREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      
      fmpz_set(a, bnum + i);
   } else
      fmpq_poly_get_coeff_fmpz(a, NF_ELEM(b), i);
}



/******************************************************************************

    Arithmetic

******************************************************************************/

NF_ELEM_INLINE
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

NF_ELEM_INLINE
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

NF_ELEM_INLINE
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

FLINT_DLL void nf_elem_add_si(nf_elem_t a, 
                                   const nf_elem_t b, slong c, const nf_t nf);
								   
FLINT_DLL void nf_elem_add_fmpz(nf_elem_t a,
                                  const nf_elem_t b, const fmpz_t c, const nf_t nf);

FLINT_DLL void nf_elem_add_fmpq(nf_elem_t a,
                                  const nf_elem_t b, const fmpq_t c, const nf_t nf);

FLINT_DLL void nf_elem_sub_si(nf_elem_t a, 
                                   const nf_elem_t b, slong c, const nf_t nf);
								   
FLINT_DLL void nf_elem_sub_fmpz(nf_elem_t a,
                                  const nf_elem_t b, const fmpz_t c, const nf_t nf);

FLINT_DLL void nf_elem_sub_fmpq(nf_elem_t a,
                                  const nf_elem_t b, const fmpq_t c, const nf_t nf);

FLINT_DLL void nf_elem_si_sub(nf_elem_t a, 
                                   slong c, const nf_elem_t b, const nf_t nf);
								   
FLINT_DLL void nf_elem_fmpz_sub(nf_elem_t a,
                                  const fmpz_t c, const nf_elem_t b, const nf_t nf);

FLINT_DLL void nf_elem_fmpq_sub(nf_elem_t a,
                                  const fmpq_t c, const nf_elem_t b, const nf_t nf);

FLINT_DLL void nf_elem_scalar_mul_si(nf_elem_t a, const nf_elem_t b, 
                                                      slong c, const nf_t nf);

FLINT_DLL void nf_elem_scalar_mul_fmpz(nf_elem_t a, const nf_elem_t b, 
                                                     const fmpz_t c, const nf_t nf);

FLINT_DLL void nf_elem_scalar_mul_fmpq(nf_elem_t a, const nf_elem_t b, 
                                                     const fmpq_t c, const nf_t nf);
									
FLINT_DLL void nf_elem_scalar_div_si(nf_elem_t a, const nf_elem_t b, 
                                                      slong c, const nf_t nf);

FLINT_DLL void nf_elem_scalar_div_fmpz(nf_elem_t a, const nf_elem_t b, 
                                                     const fmpz_t c, const nf_t nf);

FLINT_DLL void nf_elem_scalar_div_fmpq(nf_elem_t a, const nf_elem_t b, 
                                                     const fmpq_t c, const nf_t nf);
									
FLINT_DLL void _nf_elem_add_lf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can);

FLINT_DLL void _nf_elem_sub_lf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can);

FLINT_DLL void _nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can);

FLINT_DLL void _nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b, 
                                   const nf_elem_t c, const nf_t nf, int can);

FLINT_DLL void nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, 
                                            const nf_elem_t c, const nf_t nf);

FLINT_DLL void nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b, 
                                            const nf_elem_t c, const nf_t nf);

NF_ELEM_INLINE
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

NF_ELEM_INLINE
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

NF_ELEM_INLINE
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

NF_ELEM_INLINE
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

FLINT_DLL void _nf_elem_mul(nf_elem_t a, const nf_elem_t b, 
                                             const nf_elem_t c, const nf_t nf);

FLINT_DLL void nf_elem_mul(nf_elem_t a, const nf_elem_t b, 
                                             const nf_elem_t c, const nf_t nf);

FLINT_DLL void _nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, 
                                    const nf_elem_t c, const nf_t nf, int red);

FLINT_DLL void nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, 
                                    const nf_elem_t c, const nf_t nf, int red);

FLINT_DLL void _nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf);

FLINT_DLL void nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf);

FLINT_DLL void _nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf);

FLINT_DLL void nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf);

FLINT_DLL void _nf_elem_pow(nf_elem_t res, const nf_elem_t b, ulong e, const nf_t nf);

FLINT_DLL void nf_elem_pow(nf_elem_t res, const nf_elem_t a, ulong e, const nf_t nf);

FLINT_DLL void _nf_elem_norm(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf);

FLINT_DLL void nf_elem_norm(fmpq_t res, const nf_elem_t a, const nf_t nf);

FLINT_DLL void _nf_elem_norm_div(fmpz_t rnum, fmpz_t rden, const nf_elem_t a,
                             const nf_t nf, const fmpz_t divisor, slong nbits);

FLINT_DLL void nf_elem_norm_div(fmpq_t res, const nf_elem_t a, const nf_t nf,
                                            const fmpz_t divisor, slong nbits);

FLINT_DLL void _nf_elem_trace(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, 
                                                                const nf_t nf);

FLINT_DLL void nf_elem_trace(fmpq_t res, const nf_elem_t a, const nf_t nf);

#ifdef __cplusplus
}
#endif

#endif
