/**                                                     
 * Viterbi Tests
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <criterion/criterion.h>
#include <stdio.h>
#include <math.h>

#include "EmissionTable.h"
#include "HMM.h"

#include "Viterbi.h"

struct HMM* TestViterbi__hmm;
struct EmissionTable TestViterbi__et[5];
size_t TestViterbi__num_emissiontables = 0;


void TestViterbi__simpliest_setup(void) {
  EmissionTable__init(&TestViterbi__et[0], 1, UNIFORM_DISTRIBUTION);
  TestViterbi__num_emissiontables = 1;

  TestViterbi__hmm = HMM__create(1, 1, 1);

  struct State* s = HMM__new_state(TestViterbi__hmm);
  State__init_uniform(s, "1nt_uni", 1, &TestViterbi__et[0]);
  State__add_incoming(s, 0.5, s);

  struct Transition t = {.origin=s->id, .logodd=0};
  HMM__set_start(TestViterbi__hmm, t);
  t.logodd = log(0.5);
  HMM__set_end(TestViterbi__hmm, t);
}


void TestViterbi__complex_setup(void) {
  EmissionTable__init(&TestViterbi__et[0], 1, UNIFORM_DISTRIBUTION);
  EmissionTable__init(&TestViterbi__et[1], 3, UNIFORM_DISTRIBUTION);
  TestViterbi__num_emissiontables = 2;

  TestViterbi__hmm = HMM__create(2, 1, 1);

  struct State* s0 = HMM__new_state(TestViterbi__hmm);
  struct State* s1 = HMM__new_state(TestViterbi__hmm);

  State__init_uniform(s0, "1nt_uni", 1, &TestViterbi__et[0]);
  State__init_uniform(s1, "3nt_uni", 3, &TestViterbi__et[1]);

  State__add_incoming(s0, 0.5, s0);
  State__add_incoming(s0, 0.5, s1);
  State__add_incoming(s1, 0.5, s1);

  struct Transition t = {.origin=s0->id, .logodd=log(0.5)};
  HMM__set_start(TestViterbi__hmm, t);

  t.origin=s1->id;
  HMM__set_end(TestViterbi__hmm, t);
}


void TestViterbi__simpliest_teardown(void) {
  HMM__destroy(TestViterbi__hmm);

  for (size_t i=0; i < TestViterbi__num_emissiontables; i++) {
    EmissionTable__destroy(&TestViterbi__et[i]);
  }
  TestViterbi__num_emissiontables = 0;
}


void TestViterbi__validate_simple_path(void) {
  size_t path_length;
  struct State* path[3];
  Literal observations[3] = { LITERAL_A, LITERAL_A, LITERAL_A };

  Viterbi(TestViterbi__hmm, 3, observations, &path_length, path);

  for (size_t i=0; i < 3; i++) {
    struct State* state = path[i];
    /*
    char tmp[255] = "\0";
    State__str(state, tmp);
    printf("%c -> %s\n", Literal__char(observations[i]), tmp);
    */
    cr_assert_neq(state, NULL);
    cr_assert_eq(state->id, 0);
  }
}


Test(viterbi, build_and_breakdown, .init=TestViterbi__simpliest_setup, .fini=TestViterbi__simpliest_teardown) {
  cr_assert(1);
}


Test(viterbi, single_emission, .init=TestViterbi__simpliest_setup, .fini=TestViterbi__simpliest_teardown) {
  //TestViterbi__validate_simple_path();
}


Test(viterbi, complex__HMM, .init=TestViterbi__complex_setup, .fini=TestViterbi__simpliest_teardown) {

}
