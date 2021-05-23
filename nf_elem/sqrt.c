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

#define DEBUG 0

/*
   TODO:

     * Cache factorisation of f(n) on number field for future square roots
     * fix bug in fmpz_factor_trial #843
*/

#define ROT(u,v,t)   \
    do { fmpz _t = *u; *u = *v; *v = *t; *t = _t; } while (0);

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

int _nf_elem_sqrt(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      const fmpz * const bnum = LNF_ELEM_NUMREF(b);
      const fmpz * const bden = LNF_ELEM_DENREF(b);
      fmpz * const anum = LNF_ELEM_NUMREF(a);
      fmpz * const aden = LNF_ELEM_DENREF(a);
      fmpz_t r;
      int res = 1;

      if (fmpz_sgn(bnum) < 0)
         return 0;

      fmpz_init(r);
      
      fmpz_sqrtrem(anum, r, bnum);

      if (!fmpz_is_zero(r))
         res = 0;
         
      fmpz_sqrtrem(aden, r, bden);

      if (!fmpz_is_zero(r))
         res = 0;
         
      if (!res)
         nf_elem_zero(a, nf);

      fmpz_clear(r);

      return res;
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);
      const fmpz * const bden = QNF_ELEM_DENREF(b);
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      fmpz * const aden = QNF_ELEM_DENREF(a);
      fmpz_t d, t, t2, t3;
      fmpq_t r, s, r1, s1, tq, tq2;
      nf_elem_t sqr;
      int sq, ret = 1;

      /*
         We use that (m + n*sqrt(d))^2 = m^2 + n^2*d + 2*m*n*sqrt(d).
         Thus, finding (m^2)*(n^2*d) and (m^2) + (n^2*d) allows us to
         retrieve m^2 and n^2*d as the roots of (Y - m^2)(Y - n^2*d),
         from which we can retrieve m and n.

         Here d is the discriminant of the defining polynomial of the
         number field.
      */
      fmpq_init(r);
      fmpq_init(s);
      fmpq_init(r1);
      fmpq_init(s1);
      fmpq_init(tq);
      fmpq_init(tq2);
      fmpz_init(t);
      fmpz_init(t2);
      fmpz_init(t3);
      fmpz_init(d);
      nf_elem_init(sqr, nf);

      /* compute discriminant d of the defining polynomial */
      fmpz_mul(d, fmpq_poly_numref(nf->pol) + 1, fmpq_poly_numref(nf->pol) + 1);
      fmpz_mul(t, fmpq_poly_numref(nf->pol) + 0, fmpq_poly_numref(nf->pol) + 2);
      fmpz_mul_2exp(t, t, 2);
      fmpz_sub(d, d, t);

      /* write x mod nf->pol in the form r1 + s1*sqrt(d) */
      fmpz_neg(fmpq_numref(r1), fmpq_poly_numref(nf->pol) + 1);
      fmpz_mul_2exp(fmpq_denref(r1), fmpq_poly_numref(nf->pol) + 2, 1);
      fmpz_set(fmpq_denref(s1), fmpq_denref(r1));
      fmpz_set_ui(fmpq_numref(s1), 1);      
      fmpq_canonicalise(r1);
      fmpq_canonicalise(s1);

      /* write b in the form r + s*sqrt(d) */
      fmpq_set_fmpz_frac(tq, bnum + 1, bden);
      fmpq_mul(r, r1, tq);
      fmpq_mul(s, s1, tq);
      fmpq_set_fmpz_frac(tq, bnum + 0, bden);
      fmpq_add(r, r, tq);

      /* compute m^2*n^2*d and m^2 + n^2*d as above */
      fmpq_div_2exp(s, s, 1);
      fmpq_mul(s, s, s);
      fmpq_mul_fmpz(s, s, d);

      /* compute m^2 and n^2*d as above */
      fmpq_mul(tq, r, r);
      fmpq_mul_2exp(tq2, s, 2);
      fmpq_sub(tq, tq, tq2);
      if (fmpz_sgn(fmpq_numref(tq)) >= 0)
         fmpz_sqrtrem(fmpq_numref(s), t, fmpq_numref(tq));
      if (fmpz_sgn(fmpq_numref(tq)) < 0 || !fmpz_is_zero(t))
      {
         ret = 0;
         nf_elem_zero(a, nf);

         goto quadratic_cleanup;
      }

      fmpz_sqrtrem(fmpq_denref(s), t, fmpq_denref(tq));
      if (!fmpz_is_zero(t))
      {
         ret = 0;
         nf_elem_zero(a, nf);

         goto quadratic_cleanup;
      }
      fmpq_div_2exp(r, r, 1);
      fmpq_div_2exp(s, s, 1);
      
      /*
         we now have m^2 = r + s and n^2*d = r - s, or vice versa, so
         compute m^2 and n^2
      */
      fmpq_add(tq, r, s);
      fmpq_sub(s, r, s);
      fmpq_swap(r, tq);
      fmpq_div_fmpz(s, s, d);

      /* compute m and +/- n */
      if (fmpz_sgn(fmpq_numref(r)) >= 0)
      {
         fmpz_sqrtrem(t2, t, fmpq_numref(r));
         sq = fmpz_is_zero(t);
      } else
         sq = 0;
      if (sq)
      {
         fmpz_sqrtrem(t3, t, fmpq_denref(r));
         sq = fmpz_is_zero(t);
      }
      if (!sq)
      {
         fmpq_mul_fmpz(s, s, d);
         fmpq_div_fmpz(r, r, d);
         fmpq_swap(r, s);
         if (fmpz_sgn(fmpq_numref(r)) >= 0)
         {
            fmpz_sqrtrem(t2, t, fmpq_numref(r));
            sq = fmpz_is_zero(t);
         } else
            sq = 0;
         if (sq)
         {
            fmpz_sqrtrem(t3, t, fmpq_denref(r));
            sq &= fmpz_is_zero(t);
         }
         if (!sq)
         {
            ret = 0;
            nf_elem_zero(a, nf);

            goto quadratic_cleanup;
         }
      }
      fmpz_swap(fmpq_numref(r), t2);
      fmpz_swap(fmpq_denref(r), t3);
      if (fmpz_sgn(fmpq_numref(s)) >= 0)
         fmpz_sqrtrem(fmpq_numref(s), t, fmpq_numref(s));
      if (fmpz_sgn(fmpq_numref(s)) < 0 || !fmpz_is_zero(t))
      {
         ret = 0;
         nf_elem_zero(a, nf);

         goto quadratic_cleanup;
      }
      fmpz_sqrtrem(fmpq_denref(s), t, fmpq_denref(s));
      if (!fmpz_is_zero(t))
      {
         ret = 0;
         nf_elem_zero(a, nf);

         goto quadratic_cleanup;
      }

      /* write m + n*sqrt(d) as alpha*x + beta */
      fmpq_div(tq, s, s1);
      fmpq_mul(tq2, tq, r1);
      fmpq_sub(tq2, r, tq2);

      fmpz_gcd(t, fmpq_denref(tq), fmpq_denref(tq2));
      fmpz_mul(aden, fmpq_denref(tq), fmpq_denref(tq2));
      fmpz_tdiv_q(aden, aden, t);
      fmpz_tdiv_q(anum + 1, fmpq_denref(tq2), t);
      fmpz_tdiv_q(anum + 0, fmpq_denref(tq), t);
      fmpz_mul(anum + 1, anum + 1, fmpq_numref(tq));
      fmpz_mul(anum + 0, anum + 0, fmpq_numref(tq2));

      nf_elem_mul(sqr, a, a, nf);

      if (!nf_elem_equal(sqr, b, nf))
      {
         fmpq_mul_2exp(r, r, 1);
         fmpq_sub(tq2, tq2, r);

         fmpz_gcd(t, fmpq_denref(tq), fmpq_denref(tq2));
         fmpz_mul(aden, fmpq_denref(tq), fmpq_denref(tq2));
         fmpz_tdiv_q(aden, aden, t);
         fmpz_tdiv_q(anum + 1, fmpq_denref(tq2), t);
         fmpz_tdiv_q(anum + 0, fmpq_denref(tq), t);
         fmpz_mul(anum + 1, anum + 1, fmpq_numref(tq));
         fmpz_mul(anum + 0, anum + 0, fmpq_numref(tq2));
      }

quadratic_cleanup:

      fmpq_clear(r);
      fmpq_clear(s);
      fmpq_clear(r1);
      fmpq_clear(s1);
      fmpq_clear(tq);
      fmpq_clear(tq2);
      fmpz_clear(t);
      fmpz_clear(t2);
      fmpz_clear(t3);
      fmpz_clear(d);
      nf_elem_clear(sqr, nf);

      return ret;
   } else /* generic nf_elem */
   {
      const slong lenb = NF_ELEM(b)->length, lenf = fmpq_poly_length(nf->pol);
      slong nbits, bbits, primes;
      fmpq_t bnorm;
      fmpz_factor_t fac;
      flint_rand_t state;
      fmpz_t disc, z, temp, n, m, az, d, bden, g, maxd;
      nf_elem_t sqr;
      slong i, j, k;
      fmpz * r, * mr, * bz, * bz1, * modulus;
      int res = 0, factored, iters;

      if (lenb == 0)
      {
         nf_elem_zero(a, nf);
         return 1;
      }

      /* Step 1: compute norm and check it is square and rationalise denominator */
#if DEBUG
      flint_printf("Step 1\n");
#endif

      fmpq_init(bnorm);

      nf_elem_norm(bnorm, b, nf);

      if (!fmpz_is_square(fmpq_numref(bnorm)))
      {
         nf_elem_zero(a, nf);
         fmpq_clear(bnorm);
#if DEBUG
         flint_printf("Norm does not have square numerator\n");
#endif
         return 0;
      }

      if (!fmpz_is_square(fmpq_denref(bnorm)))
      {
         nf_elem_zero(a, nf);
         fmpq_clear(bnorm);
#if DEBUG
         flint_printf("Norm does not have square denominator\n");
#endif
         return 0;
      }

      fmpz_init(temp);
      fmpz_init(bden);

      fmpz_set_ui(bden, 1);

      /* get rid of denominator */
      bz1 = _fmpz_vec_init(NF_ELEM(b)->length);
      bz = _fmpz_vec_init(NF_ELEM(b)->length);

      fmpz_factor_init(fac);
      nbits = FLINT_ABS(fmpz_bits(NF_ELEM_DENREF(b)));
      primes = FLINT_MIN((slong) 4*log(nbits)*nbits + 2, 3512);
      fmpz_factor_trial(fac, NF_ELEM_DENREF(b), primes);
      
      for (i = 0; i < fac->num; i++)
      {
         if (fac->exp[i] > 1)
         {
            fmpz_pow_ui(temp, fac->p + i, fac->exp[i]/2);
            fmpz_mul(bden, bden, temp);
         }
      }

      fmpz_factor_clear(fac);

      fmpz_fdiv_q(temp, NF_ELEM_DENREF(b), bden);
      fmpz_fdiv_q(temp, temp, bden);

      _fmpz_vec_scalar_mul_fmpz(bz1, NF_ELEM_NUMREF(b), NF_ELEM(b)->length, temp);

      fmpz_mul(bden, temp, bden);

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

      bbits = FLINT_ABS(_fmpz_vec_max_bits(bz1, lenb));
      bbits -= (lenf - 1)*(FLINT_ABS(fmpz_bits(fmpq_poly_numref(nf->pol) + lenf - 1)) - 1);
      bbits = FLINT_MAX(bbits, 1);
      nbits = (2*bbits + 1)/(lenf) + 2;

      /*
         Step 3: find a nbits bit prime such that z = f(n) is a product
                 of powers of at most 14 distinct primes p_i. This ensures the
                 P_i = (p_i, x - n) are ideals of degree one in the number
                 field K defined by f. Then O_K/P_i is isomorphic to Z/pZ via
                 the map s(x) mod P_i -> s(n) mod p.
      */

      fmpz_factor_init(fac);
      fmpz_init(z);
      fmpz_init(n);
      fmpz_init(d);
      fmpz_set_ui(d, 1);

      fmpz_init(disc);
      fmpz_init(g);
      fmpz_init(maxd);
      flint_randinit(state);

      _fmpz_poly_discriminant(disc, fmpq_poly_numref(nf->pol), lenf);

      /* denominator only has primes dividing the discriminant and gcd(lc(f), tc(f)) */
      fmpz_gcd(g, fmpq_poly_numref(nf->pol) + 0, fmpq_poly_numref(nf->pol) + lenf - 1);
      fmpz_pow_ui(temp, g, FLINT_MAX((lenf - 2)*(lenf - 3) - lenf/2, 0));
      fmpz_mul(maxd, disc, temp);
      fmpz_mul(disc, disc, g);
      
      do /* continue increasing nbits until square root found */
      {
         fmpz_t fac1;
         
         fmpz_init(fac1);

         _fmpz_vec_set(bz, bz1, NF_ELEM(b)->length);

         factored = 0;
         iters = 0;

#if DEBUG
      flint_printf("Step 3\n");
#endif

         while (!factored || fac->num > 14) /* no bound known for finding such a factorisation */
         {            
            fmpz_set_ui(d, 1);
            fmpz_factor_clear(fac);
            fmpz_factor_init(fac);

            /* ensure we don't exhaust all primes of the given size */
            if (nbits < 20 && iters >= (1<<(nbits-1)))
            {
               iters = 0;
               nbits++;
            }
            
            fmpz_randbits(n, state, nbits);
            if (fmpz_sgn(n) < 0)
               fmpz_neg(n, n);

            _fmpz_poly_evaluate_fmpz(z, fmpq_poly_numref(nf->pol), lenf, n);

            primes = FLINT_MIN((slong) log(nbits)*nbits + 2, 3512);
            factored = fmpz_factor_trial(fac, z, primes);

            if (!factored)
            {
               fmpz_factor_t fac2;

               fmpz_set(fac1, fac->p + fac->num - 1);
               fac->num--;

               fmpz_factor_init(fac2);
               factored = fmpz_factor_smooth(fac2, fac1, FLINT_MIN(20, fmpz_bits(z)/100 + 1), 0);
               _fmpz_factor_concat(fac, fac2, fac->exp[fac->num]);
               fmpz_factor_clear(fac2);
            }

            if (factored)
            {
               /* check for primes dividing the discriminant */
               j = 0;
               for (i = 0; i < fac->num; i++)
               {
                  fmpz_mod(temp, disc, fac->p + i);

                  if (fmpz_is_zero(temp))
                  {
                     if (fmpz_cmp_ui(fac->p + i, 2*lenf) > 0)
                     {
                        factored = 0;

                        break;
                     } else
                     {
                        if (fac->exp[i] & 1)
                           fac->exp[i]++;
                           
                        fac->exp[i] >>= 1;
                        
                        fmpz_pow_ui(temp, fac->p + i, fac->exp[i]);
                        fmpz_mul(d, d, temp);
                        fmpz_mul(temp, temp, temp);

                        fmpz_zero(fac->p + i);
                     }
                  } else
                  {
                     if (i != j)
                     {
                        fmpz_set(fac->p + j, fac->p + i);
                        fac->exp[j] = fac->exp[i];
                     }
                     j++;
                  }
               }

               if (i == fac->num)
                  fac->num = j;
            }

            iters++;
         }

         if (!fmpz_is_one(d))
         {
            _fmpz_vec_scalar_mul_fmpz(bz, bz, NF_ELEM(b)->length, d);
            _fmpz_vec_scalar_mul_fmpz(bz, bz, NF_ELEM(b)->length, d);
         }

         fmpz_clear(fac1);

         /*
            Step 4: compute the square roots r_i of z = b(n) mod p_i^k_i for
                    each of the prime powers p_i^k_i dividing f(n). This allows
                    us to compute all of the square roots of b(n) mod m where
                    m = prod_i p_i^k_i. If m was sufficiently large, this
                    allows us to retrieve the square root of b from its n-adic
                    expansion. (Note that z here is not the same z as above!)
         */
#if DEBUG
         flint_printf("Step 4\n");
#endif

         r = _fmpz_vec_init(fac->num);

         _fmpz_poly_evaluate_fmpz(z, bz, lenb, n);

         fmpz_init(m);
         fmpz_init(az);
         modulus = _fmpz_vec_init(fac->num);

         for (i = 0; i < fac->num; i++)
         {
            fmpz_mod(temp, z, fac->p + i);

            if (!fmpz_sqrtmod(r + i, temp, fac->p + i)) /* check it is square mod p_i */
            {
               /* we must prove primality of prime used to reject as nonsquare */
               if (fmpz_is_prime(fac->p + i))
               {
                  res = 0;
                  nf_elem_zero(a, nf);
                  _fmpz_vec_clear(r, fac->num);
                  _fmpz_vec_clear(modulus, fac->num);
                  goto cleanup;
               }
            }

            fmpz_set(modulus + i, fac->p + i); /* initial modulus m = p_i */

            /*
               Step 4a: lift the square roots mod p_i to roots mod p_i^k_i
                        We require the lift to be unique so must exclude 2
                        and the case where the root is 0 mod p_i from lifting
             */
            if (fac->exp[i] != 1 && fmpz_cmp_ui(fac->p + i, 2) != 0 && !fmpz_is_zero(r + i))
            {
               for (j = 1; j < fac->exp[i]; j++)
               {
                  fmpz_add(az, r + i, r + i);
                  if (fmpz_cmp(az, modulus + i) >= 0)
                     fmpz_sub(az, az, modulus + i);
               
                  fmpz_mul(m, r + i, r + i);
                  fmpz_sub(m, m, z);
                  fmpz_mul(modulus + i, modulus + i, fac->p + i);
                  fmpz_mod(m, m, modulus + i);
                  fmpz_mul(m, m, az);
                  fmpz_mod(m, m, modulus + i);
                  fmpz_sub(r + i, r + i, m);
                  if (fmpz_sgn(r + i) < 0)
                     fmpz_add(r + i, r + i, modulus + i);
               }
            }
         }

         fmpz_clear(az);

         /* 
            Step 5: CRT recombination of the square roots r_i
                    We compute z congruent to r_i mod p_i^k_i, which is
                    hopefully k(n) where k is the square root of g. Of course
                    we must go through all possibilities of r_i and -r_i. (Note
                    again that z here has a different meaning to above.)
         */
#if DEBUG
         flint_printf("Step 5\n");
#endif

         nf_elem_init(sqr, nf);

         mr = _fmpz_vec_init(fac->num); /* compute -r_i mod p_i^k_i */

         for (i = 0; i < fac->num; i++)
         {
            if (!fmpz_is_zero(r + i))
               fmpz_sub(mr + i, modulus + i, r + i);
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
                  fmpz_CRT(z, z, m, r + i, modulus + i, 0);
               else
                  fmpz_CRT(z, z, m, mr + i, modulus + i, 0);

               fmpz_mul(m, m, modulus + i);

               k >>= 1;
            }

            /* 
               Step 6: retrieve sqrt from n-adic expansion, check sqrt.
                       First try assuming coeffs are integers.
             */
#if DEBUG
            flint_printf("Step 6\n");
#endif
            NF_ELEM(a)->length = _fmpz_poly_get_n_adic(NF_ELEM_NUMREF(a),
                                                               lenf - 1, z, n);

            fmpz_set(NF_ELEM_DENREF(a), bden);
            
            if (!fmpz_is_one(d))
            {
               fmpz_mul(NF_ELEM_DENREF(a), NF_ELEM_DENREF(a), d);
               fmpq_poly_canonicalise(NF_ELEM(a));
            }

            nf_elem_mul(sqr, a, a, nf);

            res = nf_elem_equal(sqr, b, nf);

            /*
               Step 6a. Now try rational coeffs. We loop through all
                        possible denominators of the *evaluation* k(n),
                        and try to compute the numerator as an integer
                        polynomial.
            */
            if (!res)
            {
               fmpz_t q, r, s, t, g, num, den;
               int sgn;

               fmpz_init(q);
               fmpz_init(r);
               fmpz_init(s);
               fmpz_init(t);
               fmpz_init(g);
               fmpz_init(num);
               fmpz_init(den);
               
               fmpz_set(r, m); fmpz_zero(s);
               fmpz_set(num, z); fmpz_one(den);

               while (!res && !fmpz_is_zero(num) && fmpz_cmpabs(den, maxd) <= 0) /* try rational coeffs */
               {
                  fmpz_fdiv_q(q, r, num);
                  fmpz_mul(t, q, num); fmpz_sub(t, r, t); ROT(r, num, t);
                  fmpz_mul(t, q, den); fmpz_sub(t, s, t); ROT(s, den, t);
        
                  sgn = fmpz_sgn(den) < 0;
                  
                  if (sgn)
                  {
                     fmpz_neg(num, num);
                     fmpz_neg(den, den);
                  }
                  
                  fmpz_gcd(g, num, den);
                  res = fmpz_is_one(g);

                  if (res)
                  {
                     NF_ELEM(a)->length = _fmpz_poly_get_n_adic(NF_ELEM_NUMREF(a),
                                                                     lenf - 1, num, n);
            
                     fmpz_mul(NF_ELEM_DENREF(a), bden, den);

                     if (!fmpz_is_one(d))
                     {
                        fmpz_mul(NF_ELEM_DENREF(a), NF_ELEM_DENREF(a), d);
                        fmpq_poly_canonicalise(NF_ELEM(a));
                     }

                     nf_elem_mul(sqr, a, a, nf);

                     res = nf_elem_equal(sqr, b, nf);
                  }

                  if (sgn)
                  {
                     fmpz_neg(num, num);
                     fmpz_neg(den, den);
                  }
               }

               fmpz_clear(q);
               fmpz_clear(r);
               fmpz_clear(s);
               fmpz_clear(t);
               fmpz_clear(g);
               fmpz_clear(num);
               fmpz_clear(den);
            }
         }
         
         fmpz_clear(m);
         _fmpz_vec_clear(r, fac->num);
         _fmpz_vec_clear(mr, fac->num);
         _fmpz_vec_clear(modulus, fac->num);
         nf_elem_clear(sqr, nf);

         nbits = 1.1*nbits + 1; /* increase nbits and try again if sqrt not found */
      } while (res != 1);

cleanup:
      flint_randclear(state);
      fmpz_clear(disc);
      fmpz_clear(g);
      fmpz_clear(maxd);

      fmpz_clear(n);
      fmpz_clear(d);
      fmpz_clear(temp);
      fmpz_clear(z);
      fmpz_clear(bden);

      _fmpz_vec_clear(bz, NF_ELEM(b)->length);
      _fmpz_vec_clear(bz1, NF_ELEM(b)->length);

      fmpz_factor_clear(fac);

      fmpq_clear(bnorm);

      return res;
   }
}

int nf_elem_sqrt(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   nf_elem_t t;

   if (a == b)
   {
      int ret;

      nf_elem_init(t, nf);

      ret = _nf_elem_sqrt(t, b, nf);
      nf_elem_swap(t, a, nf);

      nf_elem_clear(t, nf);

      return ret;
   }
   else
      return _nf_elem_sqrt(a, b, nf);
}