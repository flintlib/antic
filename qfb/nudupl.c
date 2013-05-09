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

    Copyright (C) 2012 William Hart

******************************************************************************/

#undef ulong /* prevent clash with stdlib */
#include <stdlib.h>
#define ulong unsigned long
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

/*void qfb_nudupl(qfb_t r, const qfb_t f, fmpz_t L)
{
   fmpz_t G, dx, dy, 
      Ax, Bx, By, Dy, x, y, 
      z, t1, t2, bx, by, 
      ax, ay, q, t, Q1;
   
   fmpz_init(G);
   fmpz_init(Ax); fmpz_init(Bx); fmpz_init(By); fmpz_init(Dy);
   fmpz_init(x); fmpz_init(y); fmpz_init(z);  
   fmpz_init(t1); fmpz_init(t2);
   fmpz_init(ax); fmpz_init(ay); fmpz_init(bx); fmpz_init(by); 
   fmpz_init(dx); fmpz_init(dy);
   fmpz_init(q); fmpz_init(t);
   fmpz_init(Q1);

   / Step 1: /
   fmpz_mod(t, f->b, f->a);
   fmpz_gcdinv(G, y, t, f->a);
   
   fmpz_set(Ax, G);
   fmpz_divexact(By, f->a, G);
   fmpz_divexact(Dy, f->b, G);
   
   / Step 2: /
   fmpz_mul(Bx, y, f->c);
   fmpz_mod(Bx, Bx, By);

   / Step 3: /
   fmpz_set(bx, Bx);
   fmpz_set(by, By);
   
   fmpz_set_ui(x, 1);
   fmpz_set_ui(y, 0);
   fmpz_set_ui(z, 0);

   while (1)
   {
      fmpz_abs(t, by);
      if (fmpz_cmp(t, L) <= 0 || fmpz_is_zero(bx))
      {
         if (fmpz_is_odd(z))
         {
            fmpz_neg(by, by);
            fmpz_neg(y, y);
         }

         fmpz_mul(ax, G, x);
         fmpz_mul(ay, G, y);

         break;
      }

      fmpz_fdiv_qr(q, t, by, bx);
      fmpz_set(by, bx);
      fmpz_set(bx, t);
      fmpz_mul(t, q, x);
      fmpz_sub(t, y, t);
      fmpz_set(y, x);
      fmpz_set(x, t);
      fmpz_add_ui(z, z, 1);
   }
   
   if (fmpz_is_zero(z))
   {
      / Step 6: /
      fmpz_mul(dx, bx, Dy);
      fmpz_sub(dx, dx, f->c);
      fmpz_divexact(dx, dx, By);
      fmpz_mul(r->a, by, by);
      fmpz_mul(r->c, bx, bx);
      fmpz_add(t, bx, by);
      fmpz_mul(t, t, t);
      fmpz_sub(r->b, f->b, t);
      fmpz_add(r->b, r->b, r->a);
      fmpz_add(r->b, r->b, r->c);
      fmpz_mul(t, G, dx);
      fmpz_sub(r->c, r->c, t);
   } else
   {
      / Step 7: /
      fmpz_mul(t1, Dy, bx);
      fmpz_mul(t2, f->c, x);
      fmpz_sub(t1, t1, t2);
      fmpz_divexact(dx, t1, By);
      fmpz_mul(Q1, y, dx);
      fmpz_add(dy, Q1, Dy);
      fmpz_add(r->b, dy, Q1);
      fmpz_mul(r->b, r->b, G);
      fmpz_divexact(dy, dy, x);
      fmpz_mul(r->a, by, by);
      fmpz_mul(r->c, bx, bx);
      fmpz_add(t, bx, by);
      fmpz_mul(t, t, t);
      fmpz_sub(r->b, r->b, t);
      fmpz_add(r->b, r->b, r->a);
      fmpz_add(r->b, r->b, r->c);
      fmpz_mul(t, ay, dy);
      fmpz_sub(r->a, r->a, t);
      fmpz_mul(t, ax, dx);
      fmpz_sub(r->c, r->c, t);
   }

   fmpz_clear(G);
   fmpz_clear(Ax); fmpz_clear(Bx); fmpz_clear(By); fmpz_clear(Dy);
   fmpz_clear(x); fmpz_clear(y); fmpz_clear(z);  
   fmpz_clear(t1); fmpz_clear(t2);
   fmpz_clear(ax); fmpz_clear(ay); fmpz_clear(bx); fmpz_clear(by); 
   fmpz_clear(dx); fmpz_clear(dy);
   fmpz_clear(q); fmpz_clear(t);
   fmpz_clear(Q1);
}*/

void qfb_nudupl(qfb_t r, const qfb_t f, fmpz_t D, fmpz_t L)
{
   fmpz_t a1, b1, c1, ca, cb, cc, k, s, t, u2, v1, v2;

   fmpz_init(a1); fmpz_init(b1); fmpz_init(c1);
   fmpz_init(ca); fmpz_init(cb); fmpz_init(cc); 
   fmpz_init(k);
   fmpz_init(s); 
   fmpz_init(t); fmpz_init(u2); fmpz_init(v1); fmpz_init(v2); 
   
   /* nucomp calculation */

   fmpz_set(a1, f->a);
   fmpz_set(c1, f->c);

   fmpz_zero(k);
   
   if (fmpz_cmpabs(b1, a1) == 0)
   {
      fmpz_set(s, a1);
      fmpz_zero(v2);
   } else if (fmpz_sgn(f->b) < 0)
   {
      fmpz_neg(b1, f->b);
      fmpz_gcdinv(s, v2, b1, a1);
      fmpz_neg(v2, v2);
   } else
      fmpz_gcdinv(s, v2, f->b, a1);

   fmpz_mul(t, v2, c1);
   fmpz_neg(k, t);

   if (!fmpz_is_one(s))
   {
      fmpz_fdiv_q(a1, a1, s);
      fmpz_mul(c1, c1, s);
   }

   fmpz_fdiv_r(k, k, a1);

   if (fmpz_cmp(a1, L) < 0)
   {
      fmpz_mul(t, a1, k);

      fmpz_mul(ca, a1, a1);

      fmpz_mul_2exp(cb, t, 1);
      fmpz_add(cb, cb, f->b);

      fmpz_add(cc, f->b, t);
      fmpz_mul(cc, cc, k);
      fmpz_add(cc, cc, c1);
      
      fmpz_fdiv_q(cc, cc, a1);
   } else
   {
      fmpz_t m2, r1, r2, co1, co2, temp;

      fmpz_init(m2); fmpz_init(r1); fmpz_init(r2);
      fmpz_init(co1); fmpz_init(co2); fmpz_init(temp);

      fmpz_set(r2, a1);
      fmpz_set(r1, k);

      fmpz_xgcd_partial(co2, co1, r2, r1, L);
      
      fmpz_mul(t, a1, r1);
      
      fmpz_mul(m2, f->b, r1);
      fmpz_mul(temp, c1, co1);
      fmpz_sub(m2, m2, temp);
      fmpz_tdiv_q(m2, m2, a1);
      
      fmpz_mul(ca, r1, r1);
      fmpz_mul(temp, co1, m2);
      if (fmpz_sgn(co1) < 0)
         fmpz_sub(ca, ca, temp);
      else
         fmpz_sub(ca, temp, ca);

      fmpz_mul(cb, ca, co2);
      fmpz_sub(cb, t, cb);
      fmpz_mul_2exp(cb, cb, 1);
      fmpz_fdiv_q(cb, cb, co1);
      fmpz_sub(cb, cb, f->b);
      fmpz_mul_2exp(temp, ca, 1);
      fmpz_fdiv_r(cb, cb, temp);
      
      fmpz_mul(cc, cb, cb);
      fmpz_sub(cc, cc, D);
      fmpz_fdiv_q(cc, cc, ca);
      fmpz_fdiv_q_2exp(cc, cc, 2);

      if (fmpz_sgn(ca) < 0)
      {
         fmpz_neg(ca, ca);
         fmpz_neg(cc, cc);
      }

      fmpz_clear(m2); fmpz_clear(r1); fmpz_clear(r2);
      fmpz_clear(co1); fmpz_clear(co2); fmpz_clear(temp);
   }

   fmpz_set(r->a, ca);
   fmpz_set(r->b, cb);
   fmpz_set(r->c, cc);

   fmpz_clear(ca); fmpz_clear(cb); fmpz_clear(cc); 
   fmpz_clear(k);  
   fmpz_clear(s); 
   fmpz_clear(t); fmpz_clear(u2); fmpz_clear(v2);
   fmpz_clear(a1); fmpz_clear(b1); fmpz_clear(c1);
}
