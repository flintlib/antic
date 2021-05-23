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

int nf_elem_is_square(const nf_elem_t b, const nf_t nf)
{
   nf_elem_t t;

   int ret;

   nf_elem_init(t, nf);

   ret = _nf_elem_sqrt(t, b, nf);

   nf_elem_clear(t, nf);

   return ret;
}
