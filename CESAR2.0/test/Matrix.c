/**                                                     
 * Matrix Test
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <criterion/criterion.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Logodd.h"

#include "Matrix.h"

struct LogoddMatrix* TestMatrix__dummy;
struct PathMatrix* TestMatrix__puppet;

Test(matrix, LogoddMatrix__str) {
  LOGODD_T dephault = .5;
  TestMatrix__dummy = LogoddMatrix__create(3, 3, dephault);

  char buffer[1024] = "";
  LogoddMatrix__str(TestMatrix__dummy, buffer);

  char cmp[1024] = "0\t+5.000000E-01\t+5.000000E-01\t+5.000000E-01\t\n1\t+5.000000E-01\t+5.000000E-01\t+5.000000E-01\t\n2\t+5.000000E-01\t+5.000000E-01\t+5.000000E-01\t\n";

  cr_assert_str_eq(buffer, cmp);
}

Test(matrix, PathMatrix__str) {
  STATE_ID_T dephault = 42;
  TestMatrix__puppet = PathMatrix__create(3, 3, dephault);

  char buffer[1024] = "";
  PathMatrix__str(TestMatrix__puppet, buffer);

  char cmp[1024] = "0\t42\t42\t42\t\n1\t42\t42\t42\t\n2\t42\t42\t42\t\n";

  cr_assert_str_eq(buffer, cmp);
}
