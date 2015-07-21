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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"


static mp_limb_t * prime_array = (void*)0;
static slong num_primes = 0, space = 0;

static void assure_primes(slong np)
{
  slong i;
  mp_limb_t p;
  mp_bitcnt_t pbits;

  np += 5; 
  if (num_primes >= np)
    return;
  
  if (np >= space) {
    mp_limb_t * new_array = malloc(2*np*sizeof(mp_limb_t));
    if (!new_array) {
      printf("no memory\n");
      exit(-1);
    }
    if (num_primes) {
      for(i=0; i<num_primes; i++)
        new_array[i] = prime_array[i];
      free(prime_array);
      prime_array = new_array;
    } else {
      /* set size of first prime */
      prime_array = new_array;
      pbits = FLINT_BITS - 1;
      p = (UWORD(1)<<pbits);
      prime_array[0] = n_nextprime(p, 0);
      num_primes = 1;
    }
    prime_array = new_array;
    space = 2*np;
  }

  for(i=num_primes; i<np; i++)
    prime_array[i] = n_nextprime(prime_array[i-1], 0);

  num_primes = np;
}

static mp_limb_t nth_prime(slong i)
{
  if (i>=num_primes)
    assure_primes(i+5);
  return prime_array[i];
}

void _fmpz_poly_resultant_modular_div(fmpz_t res, const fmpz * poly1, slong len1, 
                                        const fmpz * poly2, slong len2, const fmpz_t divisor, slong num_primes)
{
    slong i, pc;
    fmpz_comb_t comb;
    fmpz_comb_temp_t comb_temp;
    fmpz_t ac, bc, l, modulus, div;
    fmpz * A, * B, * lead_A, * lead_B;
    mp_ptr a, b, rarr, parr;
    mp_limb_t p, d;
    nmod_t mod;
    
    /* special case, one of the polys is a constant */
    if (len2 == 1) /* if len1 == 1 then so does len2 */
    {
        fmpz_pow_ui(res, poly2, len1 - 1);
        fmpz_divexact(res, res, divisor);

        return;
    }
    
    assure_primes(num_primes);

    fmpz_init(ac);
    fmpz_init(bc);
    
    /* compute content of poly1 and poly2 */
    _fmpz_vec_content(ac, poly1, len1);
    _fmpz_vec_content(bc, poly2, len2);
    
    /* divide poly1 and poly2 by their content */
    A = _fmpz_vec_init(len1);
    B = _fmpz_vec_init(len2);
    _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, ac);
    _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, bc);


    fmpz_init(l);

    if (!fmpz_is_one(ac))
    {
       fmpz_pow_ui(l, ac, len2 - 1);
       fmpz_init(div);
       fmpz_gcd(div, l, divisor);
       fmpz_divexact(div, divisor, div);
    } else {
       fmpz_init_set(div, divisor);
    }
    
    if (!fmpz_is_one(bc))
    {
       fmpz_pow_ui(l, bc, len1 - 1);
       fmpz_gcd(l, l, div);
       fmpz_divexact(div, div, l);
    }

    
    /* get product of leading coefficients */
    lead_A = A + len1 - 1;
    lead_B = B + len2 - 1;
    fmpz_mul(l, lead_A, lead_B);

    parr = _nmod_vec_init(num_primes);
    rarr = _nmod_vec_init(num_primes);

    fmpz_init(modulus);
    fmpz_set_ui(modulus, 1);
    fmpz_zero(res);

    /* make space for polynomials mod p */
    a = _nmod_vec_init(len1);
    b = _nmod_vec_init(len2);

    pc = 0;
    
    for (i = 0; i< num_primes; )
    {
        /* get new prime and initialise modulus */
        p = nth_prime(pc++);
        if (fmpz_fdiv_ui(l, p) == 0) 
            continue;
        d = fmpz_fdiv_ui(div, p);
        if (d==0)
            continue;
        d = n_invmod(d, p);

        nmod_init(&mod, p);

        /* reduce polynomials modulo p */
        _fmpz_vec_get_nmod_vec(a, A, len1, mod);
        _fmpz_vec_get_nmod_vec(b, B, len2, mod);

        /* compute resultant over Z/pZ */
        rarr[i] = _nmod_poly_resultant(a, len1, b, len2, mod);
        rarr[i] = n_mulmod2_preinv(rarr[i], d, mod.n, mod.ninv);
        parr[i++] = p;
    }

    fmpz_comb_init(comb, parr, num_primes);
    fmpz_comb_temp_init(comb_temp, comb);
    
    fmpz_multi_CRT_ui(res, rarr, comb, comb_temp, 1);
        
    fmpz_clear(modulus);
    fmpz_comb_temp_clear(comb_temp);
    fmpz_comb_clear(comb);
        
    _nmod_vec_clear(a);
    _nmod_vec_clear(b);

    _nmod_vec_clear(parr);
    _nmod_vec_clear(rarr);
    
    /* finally multiply by powers of content */
    if (!fmpz_is_one(ac))
    {
       fmpz_pow_ui(l, ac, len2 - 1);
       fmpz_mul(res, res, l);
    }
    
    if (!fmpz_is_one(bc))
    {
       fmpz_pow_ui(l, bc, len1 - 1);
       fmpz_mul(res, res, l);
    }

    fmpz_clear(l); 
    fmpz_clear(div); 
    
    _fmpz_vec_clear(A, len1);
    _fmpz_vec_clear(B, len2);

    fmpz_clear(ac);
    fmpz_clear(bc);
}

void
fmpz_poly_resultant_modular_div(fmpz_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2, const fmpz_t divisor, slong num_primes)
{
   slong len1 = poly1->length;
   slong len2 = poly2->length;
   
   if (len1 == 0 || len2 == 0)
     fmpz_zero(res);
   else if (len1 >= len2)
        _fmpz_poly_resultant_modular_div(res, poly1->coeffs, len1, poly2->coeffs, len2, divisor, num_primes);
   else
   {
        _fmpz_poly_resultant_modular_div(res, poly2->coeffs, len2, poly1->coeffs, len1, divisor, num_primes);  
        if ((len1 > 1) && (!(len1 & WORD(1)) & !(len2 & WORD(1))))
            fmpz_neg(res, res);
   }
}

