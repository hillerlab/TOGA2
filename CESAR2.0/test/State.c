/**                                                     
 * State Tests
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <criterion/criterion.h>

#include "stdio.h"

#include "EmissionTable.h"
#include "State.h"

struct State* TestState__dummy;

void State__setup(void) {
  TestState__dummy = malloc(sizeof(struct State));
}

void State__teardown(void) {
  free(TestState__dummy);
}

Test(state, namelength, .init=State__setup, .fini=State__teardown) {
  char name[STATE_NAME_LENGTH] = "0123456789012345678\0";
  State__init_silent(TestState__dummy, name);

  cr_assert_str_eq(TestState__dummy->name, name);
}

Test(state, string_representation, .init=State__setup, .fini=State__teardown) {
  State__init_silent(TestState__dummy, "Foobar");
  cr_assert_str_eq(TestState__dummy->name, "Foobar");

  TestState__dummy->id = 4711;
  cr_expect_eq(TestState__dummy->id, 4711);
  cr_expect_eq(TestState__dummy->num_incoming, 0);

  State__add_incoming(TestState__dummy, 0.5, TestState__dummy);
  cr_expect_eq(TestState__dummy->num_incoming, 1);
  cr_assert_eq(TestState__dummy->incoming[0].origin, 4711);
  cr_assert_eq(TestState__dummy->incoming[0].logodd, 0.5);
}

Test(state, namelength_toolong, .init=State__setup, .fini=State__teardown) {
  char name[41] = "0123456789012345678";
  State__init_silent(TestState__dummy, name);

  cr_assert_eq(strlen(TestState__dummy->name), STATE_NAME_LENGTH-1);
  cr_assert_str_eq(TestState__dummy->name, "0123456789012345678");
}

Test(state, string_single_emission, .init=State__setup, .fini=State__teardown) {
  State__init_silent(TestState__dummy, "Foobar");
  TestState__dummy->id = 4711;

  char buffer[255];
  State__str(TestState__dummy, buffer);

  cr_assert_str_eq(buffer, "State {.name=\"Foobar\", .id=4711, .num_emissions=0, .reference=\"\", .num_incoming=0}");
}

Test(state, string_triple_emission, .init=State__setup, .fini=State__teardown) {
  struct EmissionTable* etable = (struct EmissionTable*) malloc(sizeof(EmissionTable));
  TestState__dummy->id = 4711;
  EmissionTable__init(etable, 3, UNIFORM_DISTRIBUTION);

  State__init_uniform(TestState__dummy, "Foobarbaz", 3, etable);
  cr_assert_eq(TestState__dummy->reference[0], LITERAL_A);
  cr_assert_eq(TestState__dummy->reference[1], LITERAL_A);
  cr_assert_eq(TestState__dummy->reference[2], LITERAL_A);

  char buffer[255];
  State__str(TestState__dummy, buffer);
  cr_assert_str_eq(buffer, "State {.name=\"Foobarbaz\", .id=4711, .num_emissions=3, .reference=\"AAA\", .num_incoming=0}");

  EmissionTable__destroy(etable);
}
