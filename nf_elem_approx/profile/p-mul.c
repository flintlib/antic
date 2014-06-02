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

    Copyright 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "nf.h"
#include "nf_elem.h"
#include "nf_elem_approx.h"

#define BITS 230

typedef struct
{
   slong length;
   int monic;
} info_t;

void random_fmpq_poly(fmpq_poly_t pol, flint_rand_t state, slong length)
{
   fmpz * arr;
   slong i;

   fmpq_poly_fit_length(pol, length);

   arr = fmpq_poly_numref(pol);

   for (i = 0; i < length; i++)
      fmpz_randbits(arr + i, state, BITS);

   fmpz_randbits(fmpq_poly_denref(pol), state, BITS);

   _fmpq_poly_set_length(pol, length);
   _fmpq_poly_normalise(pol);
   fmpq_poly_canonicalise(pol);
}

void random_nf_elem(nf_elem_t a, flint_rand_t state, nf_t nf)
{
   slong len = nf->pol->length - 1;
   slong i;

   random_fmpq_poly(NF_ELEM(a), state, len);
}

void sample(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   slong length = info->length, i, j;
   int monic = info->monic;
   int scale;
   
   scale = 100;
   if (length >= 50) scale = 10;
   if (length >= 500) scale = 4;
   
   flint_rand_t state;
   flint_randinit(state);

   fmpq_poly_t pol;
   nf_t nf;
   nf_elem_t a, b;
   nf_elem_approx_t a2, b2, c2;

   fmpq_poly_init(pol);
        
   for (i = 0; i < count; i++)
   {
      random_fmpq_poly(pol, state, length);
      if (monic)
      {
         fmpz_one(fmpq_poly_denref(pol));
         fmpq_poly_set_coeff_ui(pol, length - 1, 1);
      }
	
      nf_init(nf, pol);
       
      nf_elem_init(a, nf);
      nf_elem_init(b, nf);
 
      nf_elem_approx_init(a2, nf);
      nf_elem_approx_init(b2, nf);
      nf_elem_approx_init(c2, nf);

      random_nf_elem(a, state, nf);
      random_nf_elem(b, state, nf);
      if (monic)
      {
         fmpz_one(fmpq_poly_denref(NF_ELEM(a)));
         fmpz_one(fmpq_poly_denref(NF_ELEM(b)));
      }
	
      nf_elem_approx_set_nf_elem(a2, a, nf);
      nf_elem_approx_set_nf_elem(b2, b, nf);

      prof_start();
      for (j = 0; j < scale; j++)
      {
         nf_elem_approx_mul(c2, a2, b2, nf);
      }
	   prof_stop();
   }
  
   nf_elem_approx_clear(a2, nf);
   nf_elem_approx_clear(b2, nf);
   nf_elem_approx_clear(c2, nf);

   nf_elem_clear(a, nf);
   nf_elem_clear(b, nf);
        
   nf_clear(nf);

   fmpq_poly_clear(pol);

   flint_randclear(state);
}

int main(void)
{
   double min, max;
   info_t info;
   slong k, scale;

   printf("Number field element approx multiplication\n");
   flint_printf("bits = %ld\n", BITS);

   /*for (k = 4; k <= 1000; k = (slong) ceil(1.1*k))*/
   k = 4;
   {
      info.length = k;
      info.monic = 0;

      scale = 100;
      if (k >= 50) scale = 10;
      if (k >= 500) scale = 4;
      
      prof_repeat(&min, &max, sample, (void *) &info);
      
      flint_printf("generic: length %wd, min %.3e ms, max %.3e ms\n", 
           info.length,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400.0
	     );

     info.monic = 1;
     
     prof_repeat(&min, &max, sample, (void *) &info);
         
     flint_printf("monic: length %wd, min %.3e ms, max %.3e ms\n", 
           info.length,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400.0
	     );
   }

   return 0;
}
