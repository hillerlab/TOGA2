/**                                                     
 * Logodd Tests
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <criterion/criterion.h>

#include <stdio.h>

#include "Logodd.h"


Test(logodd, exp) {
  cr_assert_eq(Logodd__exp(0), 1);
  cr_assert_eq(Logodd__exp(LOGODD_NEGINF), 0);
}


Test(logodd, log) {
  cr_assert_eq(Logodd__log(1), 0);
  cr_assert_eq(Logodd__log(0), LOGODD_NEGINF);
}


Test(logodd, add) {
  cr_assert_eq(Logodd__add(Logodd__log(0), Logodd__log(0)), LOGODD_NEGINF);
  cr_assert_eq(Logodd__add(Logodd__log(1), Logodd__log(0)), LOGODD_NEGINF);
  cr_assert_eq(Logodd__add(Logodd__log(0), Logodd__log(1)), LOGODD_NEGINF);
  cr_assert_eq(Logodd__add(LOGODD_NEGINF, .1234), LOGODD_NEGINF);
}
