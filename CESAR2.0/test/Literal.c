/**                                                     
 * Literal Tests
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <criterion/criterion.h>
#include <stdio.h>

#include "EmissionTable.h"
#include "HMM.h"

#include "Literal.h"


Test(literal, from_char) {
  cr_assert_eq(Literal__from_char('A'), LITERAL_A);
  cr_assert_eq(Literal__from_char('C'), LITERAL_C);
  cr_assert_eq(Literal__from_char('G'), LITERAL_G);
  cr_assert_eq(Literal__from_char('T'), LITERAL_T);
}


Test(literal, _char) {
  cr_assert_eq(Literal__char(LITERAL_A), 'A');
  cr_assert_eq(Literal__char(LITERAL_C), 'C');
  cr_assert_eq(Literal__char(LITERAL_G), 'G');
  cr_assert_eq(Literal__char(LITERAL_T), 'T');
}


Test(literal, str) {
  Literal codon[3] = { LITERAL_A, LITERAL_C, LITERAL_G };
  char buffer[4] = "\0";
  Literal__str(3, codon, buffer);
  cr_assert_eq(buffer[0], 'A');
  cr_assert_eq(buffer[1], 'C');
  cr_assert_eq(buffer[2], 'G');
  cr_assert_str_eq(buffer, "ACG\0");
}


Test(literal, _uint) {
  Literal sequence[3] = { LITERAL_A, LITERAL_C, LITERAL_G };
  cr_assert_eq(Literal__uint(3, sequence),     0x06);
  cr_assert_eq(Literal__uint(2, sequence),     0x01);
  cr_assert_eq(Literal__uint(1, sequence),     0x00);
  cr_assert_eq(Literal__uint(2, &sequence[1]), 0x06);
  cr_assert_eq(Literal__uint(1, &sequence[2]), 0x02);
  cr_assert_eq(Literal__uint(0, &sequence[3]), 0x00);
}
