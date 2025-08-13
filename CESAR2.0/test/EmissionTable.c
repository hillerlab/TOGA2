/**                                                     
 * EmissionTable Tests
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <criterion/criterion.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Distribution.h"

#include "EmissionTable.h"

struct EmissionTable* EmissionTable__dummy;

void EmissionTable__setup(void) {
  EmissionTable__dummy = (struct EmissionTable*) malloc(sizeof(EmissionTable));
  EmissionTable__init(EmissionTable__dummy, 1, UNIFORM_DISTRIBUTION);
}

void EmissionTable__teardown(void) {
  EmissionTable__destroy(EmissionTable__dummy);
  free(EmissionTable__dummy);
}

Test(emissiontable, by_literals, .init=EmissionTable__setup, .fini=EmissionTable__teardown) {
  EmissionTable__init(EmissionTable__dummy, 1, LAMBDA_DISTRIBUTION);
  Literal reference[1] = { LITERAL_A };
  Literal query[1] = { LITERAL_A };
  EmissionTable__set(EmissionTable__dummy, query, (LOGODD_T) .1);
  cr_assert_eq(EmissionTable__by_literals(EmissionTable__dummy, reference, query), .1);
}

Test(emissiontable, by_literals__too_many, .init=EmissionTable__setup, .fini=EmissionTable__teardown) {
  Literal literals[3] = { LITERAL_C, LITERAL_A, LITERAL_A };
  LOGODD_T logodd = EmissionTable__by_literals(EmissionTable__dummy, literals, literals);
  cr_assert_eq(log(0.25), logodd);
}

Test(emissiontable, read64, .fini=EmissionTable__teardown) {
  EmissionTable__dummy = (struct EmissionTable*) malloc(sizeof(EmissionTable));
  EmissionTable__init(EmissionTable__dummy, 3, LAMBDA_DISTRIBUTION);
  EmissionTable__read(EmissionTable__dummy, "../extra/tables/eth_codon_sub.txt");
  cr_assert_eq(EmissionTable__dummy->num_literals, (uint8_t) 3);
  cr_assert_eq(EmissionTable__get(EmissionTable__dummy, 0, 0), Logodd__log(0.40849));
  cr_assert_eq(EmissionTable__get(EmissionTable__dummy, 63, 0), Logodd__log(0.00037));
  cr_assert_eq(EmissionTable__get(EmissionTable__dummy, 0, 63), Logodd__log(0.00054));
}

Test(emissiontable, different_sizes) {
  for (int size=1; size < 3; size++) {
    struct EmissionTable* dummy = (struct EmissionTable*) malloc(sizeof(EmissionTable));
    EmissionTable__init(dummy, size, LAMBDA_DISTRIBUTION);

    /*
    char tmp[1000] = "\0";
    EmissionTable__str(dummy, tmp);
    printf("etable dummy:\n%s\n", tmp);
    */

    EmissionTable__destroy(dummy);
  }
}
