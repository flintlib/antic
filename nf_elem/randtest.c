/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013, 2014 William Hart

******************************************************************************/

#include "nf_elem.h"
#include "flint/ulong_extras.h"

void nf_elem_randtest(nf_elem_t a, flint_rand_t state, 
                                               mp_bitcnt_t bits, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_randtest(LNF_ELEM_NUMREF(a), state, bits);

        if (n_randint(state, 2))
        {
           fmpz_randtest_not_zero(LNF_ELEM_DENREF(a), state, bits);
           fmpz_abs(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(a));

           _fmpq_canonicalise(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a));
        } else
           fmpz_one(LNF_ELEM_DENREF(a));
    } else if (nf->flag & NF_QUADRATIC)
    {
        fmpz_randtest(QNF_ELEM_NUMREF(a), state, bits);
        fmpz_randtest(QNF_ELEM_NUMREF(a) + 1, state, bits);

        if (n_randint(state, 2))
        {
           fmpz_t d;
           
           fmpz_randtest_not_zero(QNF_ELEM_DENREF(a), state, bits);
           fmpz_abs(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(a));

           fmpz_init(d);
           fmpz_gcd(d, QNF_ELEM_NUMREF(a), QNF_ELEM_NUMREF(a) + 1);
           if (!fmpz_is_one(d))
           {
              fmpz_gcd(d, d, QNF_ELEM_DENREF(a));

              if (!fmpz_is_one(d))
              {
                 _fmpz_vec_scalar_divexact_fmpz(QNF_ELEM_NUMREF(a), QNF_ELEM_NUMREF(a), 2, d);
                 fmpz_divexact(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(a), d);
              }
           }
        } else
           fmpz_one(QNF_ELEM_DENREF(a));
    }
    else
    {
        fmpq_poly_randtest(NF_ELEM(a), state, nf->pol->length - 1, bits);
    }
}

void nf_elem_randtest_bounded(nf_elem_t a, flint_rand_t state, 
                                               mp_bitcnt_t bits, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        fmpz_randtest(LNF_ELEM_NUMREF(a), state, bits);

        if (n_randint(state, 2))
        {
           fmpz_randtest_not_zero(LNF_ELEM_DENREF(a), state, bits);
           fmpz_abs(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(a));

           _fmpq_canonicalise(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a));
        } else
           fmpz_one(LNF_ELEM_DENREF(a));
    } else if (nf->flag & NF_QUADRATIC)
    {
        fmpz_randtest(QNF_ELEM_NUMREF(a), state, bits);
        fmpz_randtest(QNF_ELEM_NUMREF(a) + 1, state, bits);

        if (n_randint(state, 2))
        {
           fmpz_t d;
           
           fmpz_randtest_not_zero(QNF_ELEM_DENREF(a), state, bits);
           fmpz_abs(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(a));

           fmpz_init(d);
           fmpz_gcd(d, QNF_ELEM_NUMREF(a), QNF_ELEM_NUMREF(a) + 1);
           if (!fmpz_is_one(d))
           {
              fmpz_gcd(d, d, QNF_ELEM_DENREF(a));

              if (!fmpz_is_one(d))
              {
                 _fmpz_vec_scalar_divexact_fmpz(QNF_ELEM_NUMREF(a), QNF_ELEM_NUMREF(a), 2, d);
                 fmpz_divexact(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(a), d);
              }
           }
        } else
           fmpz_one(QNF_ELEM_DENREF(a));
    }
    else
    {
        slong i, lenf = nf->pol->length;
        fmpz_t d;

        for (i = 0; i < lenf - 1; i++)
            fmpz_randtest(NF_ELEM_NUMREF(a) + i, state, bits);

        if (n_randint(state, 2))                                                                                                {
            fmpz_init(d);

            fmpz_randtest_not_zero(NF_ELEM_DENREF(a), state, bits);
            fmpz_abs(NF_ELEM_DENREF(a), NF_ELEM_DENREF(a));

            _fmpz_vec_content(d, NF_ELEM_NUMREF(a), lenf - 1);

            if (!fmpz_is_one(d))
            {
                fmpz_gcd(d, d, NF_ELEM_DENREF(a));

                if (!fmpz_is_one(d))
                {
                    _fmpz_vec_scalar_divexact_fmpz(NF_ELEM_NUMREF(a), NF_ELEM_NUMREF(a), lenf - 1, d);
                    fmpz_divexact(NF_ELEM_DENREF(a), NF_ELEM_DENREF(a), d);
                }
            }

            fmpz_clear(d);
        } else
            fmpz_one(NF_ELEM_DENREF(a));

        _fmpq_poly_set_length(NF_ELEM(a), lenf - 1);
        _fmpq_poly_normalise(NF_ELEM(a));
    }
}

void nf_elem_randtest_not_zero(nf_elem_t a, flint_rand_t state, 
                                               mp_bitcnt_t bits, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
       do {
          nf_elem_randtest(a, state, bits, nf);
       } while (fmpz_is_zero(QNF_ELEM_NUMREF(a)));
   } else if (nf->flag & NF_QUADRATIC)
   {
       do {
          nf_elem_randtest(a, state, bits, nf);
       } while (fmpz_is_zero(QNF_ELEM_NUMREF(a)) && fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1));
   } else
   {
      do {
          nf_elem_randtest(a, state, bits, nf);
       } while (fmpq_poly_is_zero(NF_ELEM(a)));
   }
}
