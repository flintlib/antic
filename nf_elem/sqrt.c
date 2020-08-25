/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2020 William Hart

******************************************************************************/

#include "nf_elem.h"

#define DEBUG 1

/*
   TODO:

     * Handle aliasing
     * try to reuse information from previous failed attempt
     * improve bounds
     * add LM bound termination for nonsquare case
     * deal with algebraic integers with denominators
     * add linear and quadratic cases
*/

slong _fmpz_poly_get_n_adic(fmpz * sqrt, slong len, fmpz_t z, fmpz_t n)
{
   fmpz_t n2, z2;
   slong i, slen = 0;

   fmpz_init(n2);
   fmpz_init(z2);

   fmpz_fdiv_q_2exp(n2, n, 1);
      
   if (fmpz_sgn(z) < 0)
      fmpz_neg(z2, z);
   else
      fmpz_set(z2, z);

   for (i = 0; i < len; i++)
   {
      fmpz_mod(sqrt + i, z2, n);

      fmpz_sub(z2, z2, sqrt + i);

      if (fmpz_cmpabs(sqrt + i, n2) > 0)
      {
         fmpz_sub(sqrt + i, sqrt + i, n);
         fmpz_add(z2, z2, n);
      }

      fmpz_fdiv_q(z2, z2, n);

      if (!fmpz_is_zero(sqrt + i))
         slen = i + 1;
   }

   fmpz_clear(n2);
   fmpz_clear(z2);

   return slen;
}

int nf_elem_sqrt(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      /* const fmpz * const bnum = LNF_ELEM_NUMREF(b);
      const fmpz * const bden = LNF_ELEM_DENREF(b);
      fmpz * const anum = LNF_ELEM_NUMREF(a);
      fmpz * const aden = LNF_ELEM_DENREF(a); */
      
      flint_printf("Sqrt for linear number fields not implemented yet\n");
      flint_abort();
   } else if (nf->flag & NF_QUADRATIC)
   {
      /* const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      const fmpz * const bden = QNF_ELEM_DENREF(b);
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      fmpz * const aden = QNF_ELEM_DENREF(a); */
      
      flint_printf("Sqrt for quadratic number fields not implemented yet\n");
      flint_abort();
   } else /* generic nf_elem */
   {
      const slong lenb = NF_ELEM(b)->length, lenf = fmpq_poly_length(nf->pol);
      slong nbits, bbits;
      fmpq_t bnorm;
      fmpz_factor_t fac;
      flint_rand_t state;
      fmpz_t disc, z, temp, n, m;
      nf_elem_t sqr;
      slong i, j, k;
      fmpz * r, * mr;
      int res = 0, factored;

      if (!fmpz_is_one(NF_ELEM_DENREF(b)))
      {
         flint_printf("Sqrt1 outside Z[alpha] not implemented yet\n");
         flint_abort();
      }

      if (lenb == 0)
      {
         nf_elem_zero(a, nf);
         return 1;
      }

      /* Step 1: compute norm and check it is square */
#if DEBUG
      flint_printf("Step 1\n");
#endif

      fmpq_init(bnorm);

      nf_elem_norm(bnorm, b, nf);
      
      if (!fmpz_is_one(fmpq_denref(bnorm)))
      {
         flint_printf("Sqrt2 outside Z[alpha] not yet implemented yet\n");
         fmpq_clear(bnorm);
         flint_abort();
      }

      if (!fmpz_is_square(fmpq_numref(bnorm)))
      {
         nf_elem_zero(a, nf);
         fmpq_clear(bnorm);
         return 0;
      }

      /* 
         Step 2: compute number of bits for initial n
                 start by assuming sqrt k has coeffs of about b1 bits where
                 b1 = about half the bits of max_i{abs(coeff(b, i))}
                 we need nbits >= b1 + 1 (extra 1 for sign) and
                 bits(m) >= bits(k(n)) where m = f(n) for f = nf->pol
                 Note f has higher degree than k so this should be true
                 if bits(n) = b1 + 1.
      */
#if DEBUG
      flint_printf("Step 2\n");
#endif

      bbits = FLINT_ABS(_fmpz_vec_max_bits(NF_ELEM_NUMREF(b), lenb));
      nbits = (bbits + 1)/2 + 2;

      /*
         Step 3: find a nbits bit prime such that z = f(n) is a product
                 at most five distinct primes p_i, none of which divide the
                 discriminant of f or the norm of b. This guarantees that
                 P_i = (p_i, x - n) are ideals of degree one  in the number
                 field K defined by f. Then O_K/P_i is isomorphic to Z/pZ via
                 the map s(x) mod P_i -> s(n) mod p.
      */
#if DEBUG
      flint_printf("Step 3\n");
#endif

      fmpz_factor_init(fac);
      fmpz_init(z);
      fmpz_init(temp);
      fmpz_init(n);

      do /* continue increasing nbits until square root found */
      {
         fmpz_init(disc);
         flint_randinit(state);

         _fmpz_poly_discriminant(disc, fmpq_poly_numref(nf->pol), lenf);

         factored = 0;
         
         while (!factored || fac->num > 5) /* no bound known for finding such a factorisation */
         {
            fmpz_factor_clear(fac);
            fmpz_factor_init(fac);

            fmpz_randprime(n, state, nbits, 0);
            
            _fmpz_poly_evaluate_fmpz(z, fmpq_poly_numref(nf->pol), lenf, n);

            factored = fmpz_factor_trial(fac, z, 3512);

            if (!factored)
               factored = fmpz_is_probabprime(fac->p + fac->num - 1);
         }

         flint_randclear(state);
         fmpz_clear(disc);

         /*
            Step 4: compute the square roots r_i of z = b(n) mod p_i for each
                    of the primes p_i dividing f(n). This allows us to compute
                    all of the square roots of b(n) mod m where m = prod_i p_i.
                    If m was sufficiently large, this allows us to retrieve the
                    square root of b from its n-adic expansion. (Note that z
                    here is not the same z as above!)
         */
#if DEBUG
         flint_printf("Step 4\n");
#endif

         r = _fmpz_vec_init(fac->num);

         _fmpz_poly_evaluate_fmpz(z, NF_ELEM_NUMREF(b), lenb, n);

         for (i = 0; i < fac->num; i++)
         {
            fmpz_mod(temp, z, fac->p + i);

            if (!fmpz_sqrtmod(r + i, temp, fac->p + i)) /* check it is square mod p_i */
            {
               res = 0;
               nf_elem_zero(a, nf);
               goto cleanup;
            }
         }

         /* 
            Step 5: CRT recombination of the square roots r_i
                    We compute z congruent to r_i mod p_i, which is hopefully
                    k(n) where k is the square root of g. Of course we must
                    go through all possibilities of r_i and -r_i. (Note again
                    that z here has a different meaning to above.)
         */
#if DEBUG
         flint_printf("Step 5\n");
#endif

         nf_elem_init(sqr, nf);

         mr = _fmpz_vec_init(fac->num); /* compute -r_i mod p_i */
         fmpz_init(m);

         for (i = 0; i < fac->num; i++)
         {
            if (!fmpz_is_zero(r + i))
               fmpz_sub(mr + i, fac->p + i, r + i);
         }

         for (j = 0; j < (1<<fac->num) && res != 1; j++)
         {
            k = j;

            fmpz_set_ui(z, 0);
            fmpz_set_ui(m, 1);
            nf_elem_zero(a, nf);
   
            for (i = 0; i < fac->num; i++)
            {
               if (k & 1)
                  fmpz_CRT(z, z, m, r + i, fac->p + i, 0);
               else
                  fmpz_CRT(z, z, m, mr + i, fac->p + i, 0);

               fmpz_mul(m, m, fac->p + i);

               k >>= 1;
            }

            /* Step 6: retrieve sqrt from n-adic expansion, check sqrt */
#if DEBUG
            flint_printf("Step 6\n");
#endif
            NF_ELEM(a)->length = _fmpz_poly_get_n_adic(NF_ELEM_NUMREF(a),
                                                               lenf - 1, z, n);

            nf_elem_mul(sqr, a, a, nf);

            res = nf_elem_equal(sqr, b, nf);
         }
         
         fmpz_clear(m);
         _fmpz_vec_clear(mr, fac->num);

         nf_elem_clear(sqr, nf);

         nbits = 1.1*nbits + 1; /* increase nbits and try again if sqrt not found */
      } while (res != 1);

cleanup:
      fmpz_clear(n);
      fmpz_clear(temp);
      fmpz_clear(z);

      _fmpz_vec_clear(r, fac->num);

      fmpz_factor_clear(fac);

      fmpq_clear(bnorm);

      return res;
   }
}

