/**
 * Test reading profiles.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <stdio.h>
#include <criterion/criterion.h>

#include "Logodd.h"

#include "Profile.h"


struct Profile* TestProfile__dummy;


void TestProfile__setup() {
  TestProfile__dummy = Profile__create("dummy");
  Profile__read(TestProfile__dummy, "../extra/tables/human/do_profile.txt");
}


void TestProfile__teardown() {
  Profile__destroy(TestProfile__dummy);
}


Test(profile, create_and_destroy, .init=TestProfile__setup, .fini=TestProfile__teardown) {
}


Test(profile, read_sequence, .init=TestProfile__setup, .fini=TestProfile__teardown) {
  Literal profile[5] = {LITERAL_A, LITERAL_G, LITERAL_G, LITERAL_G, LITERAL_G};
  cr_assert_eq(Profile__by_literals(TestProfile__dummy, profile), LOGODD_NEGINF);
}


Test(profile, read, .init=TestProfile__setup, .fini=TestProfile__teardown) {
  cr_assert_eq(EmissionTable__get(&TestProfile__dummy->emission_tables[0], LITERAL_A, 0), LOGODD_NEGINF);
  cr_assert_eq(EmissionTable__get(&TestProfile__dummy->emission_tables[0], LITERAL_G, 0), 0);
  cr_assert_eq(EmissionTable__get(&TestProfile__dummy->emission_tables[5], LITERAL_T, 0), Logodd__log(0.483));
}
