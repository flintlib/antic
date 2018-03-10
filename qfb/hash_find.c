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

#include <stdlib.h>
#include <gmp.h>
#include "flint/fmpz.h"
#include "qfb.h"

slong qfb_hash_find(qfb_hash_t * qhash, qfb_t q, slong depth)
{
   slong size = (1L<<depth), i;
   fmpz_t r;

   fmpz_init(r);

   fmpz_fdiv_r_2exp(r, q->a, depth);
   i = fmpz_get_ui(r);

   while (!fmpz_is_zero(qhash[i].q->a))
   {
      if (fmpz_cmp(qhash[i].q->a, q->a) == 0)
      {
         if (fmpz_cmpabs(qhash[i].q->b, q->b) == 0)
         {
            fmpz_clear(r);
            return i;
         }
      }
      
      i++;
      if (i == size)
         i = 0;
   }

   fmpz_clear(r);
   return -1;
}
