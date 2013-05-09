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

/*void qfb_nucomp(qfb_t r, qfb_t f, qfb_t g, fmpz_t L)
{
   int Fdivs;
   
   fmpz_t s, m, a, b, c, F, G, H, 
      Ax, Bx, By, Cy, Dy, x, y, 
      z, l, t1, t2, bx, by, 
      ax, ay, q, t, Q1, Q2, Q3, Q4,
      dx, dy, cx, cy;
   
   fmpz_init(s); fmpz_init(l); fmpz_init(m);
   fmpz_init(a); fmpz_init(c); fmpz_init(b);
   fmpz_init(F); fmpz_init(G); fmpz_init(H);
   fmpz_init(Ax); fmpz_init(Bx); fmpz_init(By); 
   fmpz_init(Cy); fmpz_init(Dy);
   fmpz_init(x); fmpz_init(y); fmpz_init(z);  
   fmpz_init(t1); fmpz_init(t2);
   fmpz_init(ax); fmpz_init(ay); fmpz_init(bx); fmpz_init(by); 
   fmpz_init(cx); fmpz_init(cy); fmpz_init(dx); fmpz_init(dy);
   fmpz_init(q); fmpz_init(t);
   fmpz_init(Q1); fmpz_init(Q2); fmpz_init(Q3); fmpz_init(Q4);
   
   / Step 1: /
   if (fmpz_cmp(f->c, g->c) < 0)
   {
      qfb * temp = f;
      f = g;
      g = temp;
   }

   fmpz_add(s, f->b, g->b);
   fmpz_fdiv_q_2exp(s, s, 1);
   fmpz_sub(m, g->b, s);

   / Step 2: /
   fmpz_xgcd(F, b, c, g->a, f->a);
   if ((Fdivs = fmpz_divisible(s, F)))
   {
      fmpz_set(G, F);
      fmpz_set(Ax, G);
      fmpz_mul(Bx, m, b);
   } else
   {
      / Step 3: /
      fmpz_mod(t1, s, F);
      fmpz_gcdinv(G, y, t1, F);
      fmpz_divexact(H, F, G);
   }

   fmpz_divexact(By, f->a, G);
   fmpz_divexact(Cy, g->a, G);
   fmpz_divexact(Dy, s, G);
   
   / Step 4 /
   if (!Fdivs)
   {
      fmpz_mod(t1, f->c, H);
      fmpz_mul(t1, t1, b);
      fmpz_mod(t2, g->c, H);
      fmpz_mul(t2, t2, c);
      fmpz_add(t1, t1, t2);
      fmpz_mod(t1, t1, H);
      fmpz_mul(t1, t1, y);
      fmpz_mod(l, t1, H);

      fmpz_divexact(t1, m, H);
      fmpz_mul(t1, t1, b);
      fmpz_divexact(t2, By, H);
      fmpz_mul(t2, t2, l);
      fmpz_add(Bx, t1, t2);
   }

   / Step 5: /
   fmpz_mod(bx, Bx, By);
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
      fmpz_mul(Q1, Cy, bx);
      fmpz_sub(cx, Q1, m);
      fmpz_divexact(cx, cx, By);
      fmpz_mul(dx, bx, Dy);
      fmpz_sub(dx, dx, g->c);
      fmpz_divexact(dx, dx, By);
      fmpz_mul(r->a, by, Cy);
      fmpz_mul(t1, bx, cx);
      fmpz_mul(t2, G, dx);
      fmpz_sub(r->c, t1, t2);
      fmpz_set(t1, g->b);
      fmpz_mul_2exp(r->b, Q1, 1);
      fmpz_sub(r->b, t1, r->b);
   } else
   {
      / Step 7: /
      fmpz_mul(t1, Cy, bx);
      fmpz_mul(t2, m, x);
      fmpz_sub(t1, t1, t2);
      fmpz_divexact(cx, t1, By);
      fmpz_mul(Q1, by, cx);
      fmpz_add(Q2, Q1, m);
      fmpz_mul(t1, Dy, bx);
      fmpz_mul(t2, g->c, x);
      fmpz_sub(t1, t1, t2);
      fmpz_divexact(dx, t1, By);
      fmpz_mul(Q3, y, dx);
      fmpz_add(Q4, Q3, Dy);
      fmpz_divexact(dy, Q4, x);
      if (!fmpz_is_zero(bx))
         fmpz_divexact(cy, Q2, bx);
      else
      {
         fmpz_mul(cy, cx, dy);
         fmpz_sub(cy, cy, f->c);
         fmpz_divexact(cy, cy, dx);
      }
      fmpz_mul(t1, by, cy);
      fmpz_mul(t2, ay, dy);
      fmpz_sub(r->a, t1, t2);
      fmpz_mul(t1, bx, cx);
      fmpz_mul(t2, ax, dx);
      fmpz_sub(r->c, t1, t2);
      fmpz_add(r->b, Q3, Q4);
      fmpz_mul(r->b, r->b, G);
      fmpz_sub(r->b, r->b, Q1);
      fmpz_sub(r->b, r->b, Q2);
   }

   fmpz_clear(s); fmpz_clear(l); fmpz_clear(m);
   fmpz_clear(a); fmpz_clear(b); fmpz_clear(c);
   fmpz_clear(F); fmpz_clear(G); fmpz_clear(H);
   fmpz_clear(Ax); fmpz_clear(Bx); fmpz_clear(By); 
   fmpz_clear(Cy); fmpz_clear(Dy);
   fmpz_clear(x); fmpz_clear(y); fmpz_clear(z);  
   fmpz_clear(t1); fmpz_clear(t2);
   fmpz_clear(ax); fmpz_clear(ay); fmpz_clear(bx); fmpz_clear(by); 
   fmpz_clear(cx); fmpz_clear(cy); fmpz_clear(dx); fmpz_clear(dy);
   fmpz_clear(q); fmpz_clear(t);
   fmpz_clear(Q1); fmpz_clear(Q2); fmpz_clear(Q3); fmpz_clear(Q4);
}*/

void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, fmpz_t r2, fmpz_t r1, fmpz_t L)
{
   fmpz_t q, r;
   long aa2, aa1, bb2, bb1, rr1, rr2, qq, bb, t1, t2, t3, i;
   long bits, bits1, bits2;

   fmpz_init(q); fmpz_init(r);
   
   fmpz_zero(co2);
   fmpz_set_si(co1, -1);
  
   while (!fmpz_is_zero(r1) && fmpz_cmp(r1, L) > 0)
   {
      bits2 = fmpz_bits(r2);
      bits1 = fmpz_bits(r1);
      bits = FLINT_MAX(bits2, bits1) - FLINT_BITS + 1;
      if (bits < 0) bits = 0;
      
      fmpz_tdiv_q_2exp(r, r2, bits);
      rr2 = fmpz_get_ui(r);
      fmpz_tdiv_q_2exp(r, r1, bits);
      rr1 = fmpz_get_ui(r);
      fmpz_tdiv_q_2exp(r, L, bits);
      bb = fmpz_get_ui(r);

      aa2 = 0; aa1 = 1;
      bb2 = 1; bb1 = 0;

      for (i = 0; rr1 != 0 && rr1 > bb; i++)
      {
         qq = rr2 / rr1;

         t1 = rr2 - qq*rr1; 
         t2 = aa2 - qq*aa1; 
         t3 = bb2 - qq*bb1; 

         if (i & 1)
         {
            if (t1 < -t3 || rr1 - t1 < t2 - aa1) break;
         } else 
         {
            if (t1 < -t2 || rr1 - t1 < t3 - bb1) break;
         }

         rr2 = rr1; rr1 = t1;
         aa2 = aa1; aa1 = t2;
         bb2 = bb1; bb1 = t3;
      }

      if (i == 0)
      {
         fmpz_fdiv_qr(q, r2, r2, r1);
         fmpz_swap(r2, r1);

         fmpz_submul(co2, co1, q);
         fmpz_swap(co2, co1);
      } else
      {
         fmpz_mul_si(r, r2, bb2);
         if (aa2 >= 0)
            fmpz_addmul_ui(r, r1, aa2);
         else
            fmpz_submul_ui(r, r1, -aa2);
         fmpz_mul_si(r1, r1, aa1);
         if (bb1 >= 0)
            fmpz_addmul_ui(r1, r2, bb1);
         else
            fmpz_submul_ui(r1, r2, -bb1);
         fmpz_set(r2, r);

         fmpz_mul_si(r, co2, bb2);
         if (aa2 >= 0)
            fmpz_addmul_ui(r, co1, aa2);
         else
            fmpz_submul_ui(r, co1, -aa2);
         fmpz_mul_si(co1, co1, aa1);
         if (bb1 >= 0)
            fmpz_addmul_ui(co1, co2, bb1);
         else
            fmpz_submul_ui(co1, co2, -bb1);
         fmpz_set(co2, r);

         if (fmpz_sgn(r1) < 0) { fmpz_neg(co1, co1); fmpz_neg(r1, r1); }
         if (fmpz_sgn(r2) < 0) { fmpz_neg(co2, co2); fmpz_neg(r2, r2); }
      }
   }

   if (fmpz_sgn(r2) < 0)
   { 
      fmpz_neg(co2, co2); fmpz_neg(co1, co1);
      fmpz_neg(r2, r2);
   }

   fmpz_clear(q); fmpz_clear(r);
}

/*
   This is based (with permission) on an implementation of NUCOMP 
   by Michael Jacobson which is an adaption to binary quadratic 
   forms of the algorithm presented in Appendix A of "Solving the 
   Pell Equation", Michael Jacobson and High Williams, CMS Books in 
   Mathematics, Springer 2009.
*/
void qfb_nucomp(qfb_t r, const qfb_t f, const qfb_t g, fmpz_t L)
{
   fmpz_t n;

   fmpz_t a1, a2, c2, ca, cb, cc, k, s, sp, ss, m, t, u2, v1, v2;

   if (fmpz_cmp(f->a, g->a) > 0)
   {
      qfb_nucomp(r, g, f, L);
      return;
   }

   fmpz_init(a1); fmpz_init(a2); fmpz_init(c2);
   fmpz_init(ca); fmpz_init(cb); fmpz_init(cc); 
   fmpz_init(k); fmpz_init(m); 
   fmpz_init(s); fmpz_init(sp); fmpz_init(ss); 
   fmpz_init(t); fmpz_init(u2); fmpz_init(v1); fmpz_init(v2); 
   fmpz_init(n);

   /* compute n = b^2 - 4ac */
   fmpz_mul(n, f->b, f->b);
   fmpz_mul(t, f->a, f->c);
   fmpz_mul_2exp(t, t, 2);
   fmpz_sub(n, n, t);

   /* nucomp calculation */

   fmpz_set(a1, f->a);
   fmpz_set(a2, g->a);
   fmpz_set(c2, g->c);

   fmpz_add(ss, f->b, g->b);
   fmpz_fdiv_q_2exp(ss, ss, 1);

   fmpz_sub(m, f->b, g->b);
   fmpz_fdiv_q_2exp(m, m, 1);

   fmpz_fdiv_r(t, a2, a1);
   if (fmpz_is_zero(t))
   {
      fmpz_set_ui(v1, 0);
      fmpz_set(sp, a1);
   } else
      fmpz_gcdinv(sp, v1, t, a1);

   fmpz_mul(k, m, v1);
   fmpz_fdiv_r(k, k, a1);
   
   if (!fmpz_is_one(sp))
   {
      fmpz_xgcd(s, v2, u2, ss, sp);
      
      fmpz_mul(k, k, u2);
      fmpz_mul(t, v2, c2);
      fmpz_sub(k, k, t);

      if (!fmpz_is_one(s))
      {
         fmpz_fdiv_q(a1, a1, s);
         fmpz_fdiv_q(a2, a2, s);
         fmpz_mul(c2, c2, s);
      }

      fmpz_fdiv_r(k, k, a1);
   }

   if (fmpz_cmp(a1, L) < 0)
   {
      fmpz_mul(t, a2, k);

      fmpz_mul(ca, a2, a1);

      fmpz_mul_2exp(cb, t, 1);
      fmpz_add(cb, cb, g->b);

      fmpz_add(cc, g->b, t);
      fmpz_mul(cc, cc, k);
      fmpz_add(cc, cc, c2);
      
      fmpz_fdiv_q(cc, cc, a1);
   } else
   {
      fmpz_t m1, m2, r1, r2, co1, co2, temp;

      fmpz_init(m1); fmpz_init(m2); fmpz_init(r1); fmpz_init(r2);
      fmpz_init(co1); fmpz_init(co2); fmpz_init(temp);

      fmpz_set(r2, a1);
      fmpz_set(r1, k);

      fmpz_xgcd_partial(co2, co1, r2, r1, L);
      
      fmpz_mul(t, a2, r1);
      fmpz_mul(m1, m, co1);
      fmpz_add(m1, m1, t);
      fmpz_tdiv_q(m1, m1, a1);
      
      fmpz_mul(m2, ss, r1);
      fmpz_mul(temp, c2, co1);
      fmpz_sub(m2, m2, temp);
      fmpz_tdiv_q(m2, m2, a1);
      
      fmpz_mul(ca, r1, m1);
      fmpz_mul(temp, co1, m2);
      if (fmpz_sgn(co1) < 0)
         fmpz_sub(ca, ca, temp);
      else
         fmpz_sub(ca, temp, ca);

      fmpz_mul(cb, ca, co2);
      fmpz_sub(cb, t, cb);
      fmpz_mul_2exp(cb, cb, 1);
      fmpz_fdiv_q(cb, cb, co1);
      fmpz_sub(cb, cb, g->b);
      fmpz_mul_2exp(temp, ca, 1);
      fmpz_fdiv_r(cb, cb, temp);
      
      fmpz_mul(cc, cb, cb);
      fmpz_sub(cc, cc, n);
      fmpz_fdiv_q(cc, cc, ca);
      fmpz_fdiv_q_2exp(cc, cc, 2);

      if (fmpz_sgn(ca) < 0)
      {
         fmpz_neg(ca, ca);
         fmpz_neg(cc, cc);
      }

      fmpz_clear(m1); fmpz_clear(m2); fmpz_clear(r1); fmpz_clear(r2);
      fmpz_clear(co1); fmpz_clear(co2); fmpz_clear(temp);
   }

   fmpz_set(r->a, ca);
   fmpz_set(r->b, cb);
   fmpz_set(r->c, cc);

   fmpz_clear(ca); fmpz_clear(cb); fmpz_clear(cc); 
   fmpz_clear(k); fmpz_init(m); 
   fmpz_clear(s); fmpz_clear(sp); fmpz_clear(ss); 
   fmpz_clear(t); fmpz_init(u2); fmpz_init(v1); fmpz_init(v2);
   fmpz_clear(a1); fmpz_clear(a2); fmpz_clear(c2);
   fmpz_clear(n);
}