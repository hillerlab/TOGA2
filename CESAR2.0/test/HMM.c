/**
 * HMM Tests
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include <criterion/criterion.h>
#include "stdio.h"
#include "math.h"

#include "State.h"
#include "Transition.h"
#include "EmissionTable.h"

#include "HMM.h"

struct HMM* TestHMM__hmm;
struct EmissionTable TestHMM__et[1];


void HMM__setup(void) {
  EmissionTable__init(&TestHMM__et[0], 1, UNIFORM_DISTRIBUTION);

  TestHMM__hmm = HMM__create(3, 3, 3);

  struct State* rock = HMM__new_state(TestHMM__hmm);
  struct State* paper = HMM__new_state(TestHMM__hmm);
  struct State* scissors = HMM__new_state(TestHMM__hmm);

  State__init_uniform(rock,     "Rock",     1, &TestHMM__et[0]);
  State__init_uniform(paper,    "Paper",    1, &TestHMM__et[0]);
  State__init_uniform(scissors, "Scissors", 1, &TestHMM__et[0]);

  State__add_incoming(rock, log(.33), rock);
  State__add_incoming(rock, log(.33), paper);
  State__add_incoming(rock, log(.33), scissors);
  State__add_incoming(paper, log(.33), rock);
  State__add_incoming(paper, log(.33), paper);
  State__add_incoming(paper, log(.33), scissors);
  State__add_incoming(scissors, log(.33), rock);
  State__add_incoming(scissors, log(.33), paper);
  State__add_incoming(scissors, log(.33), scissors);

  struct Transition t = {.origin=rock->id, .logodd=log(.33)};
  HMM__set_start(TestHMM__hmm, t);
  HMM__set_end(TestHMM__hmm, t);
  t.origin = paper->id;
  HMM__set_start(TestHMM__hmm, t);
  HMM__set_end(TestHMM__hmm, t);
  t.origin = scissors->id;
  HMM__set_start(TestHMM__hmm, t);
  HMM__set_end(TestHMM__hmm, t);
}


void HMM__teardown(void) {
  HMM__destroy(TestHMM__hmm);
  EmissionTable__destroy(&TestHMM__et[0]);
}


Test(hmm, HMM__new_state) {
  TestHMM__hmm = HMM__create(3, 0, 0);

  struct State* one = HMM__new_state(TestHMM__hmm);
  struct State* two = HMM__new_state(TestHMM__hmm);
  struct State* three = HMM__new_state(TestHMM__hmm);

  cr_assert_eq(one->id, 0);
  cr_assert_eq(two->id, 1);
  cr_assert_eq(three->id, 2);

  cr_assert_eq(TestHMM__hmm->states[0].id, 0);
  cr_assert_eq(TestHMM__hmm->states[1].id, 1);
  cr_assert_eq(TestHMM__hmm->states[2].id, 2);

  cr_assert_eq(one->num_incoming, 0);
  cr_assert_eq(two->num_incoming, 0);
}


Test(hmm, names_and_origins, .init=HMM__setup, .fini=HMM__teardown) {
  cr_assert_str_eq(TestHMM__hmm->states[0].name, "Rock");
  cr_assert_str_eq(TestHMM__hmm->states[1].name, "Paper");
  cr_assert_str_eq(TestHMM__hmm->states[2].name, "Scissors");

  char buffer[255];
  State__str(&TestHMM__hmm->states[0], buffer);

  cr_assert_eq(TestHMM__hmm->states[0].id, 0);
  cr_assert_eq(TestHMM__hmm->states[1].id, 1);
  cr_assert_eq(TestHMM__hmm->states[2].id, 2);

  cr_assert_eq(TestHMM__hmm->states[0].num_incoming, 3);
  cr_assert_eq(TestHMM__hmm->states[0].incoming[0].origin, 0);
  cr_assert_eq(TestHMM__hmm->states[0].incoming[0].logodd, log(0.33));

  cr_assert_str_eq(TestHMM__hmm->states[TestHMM__hmm->states[0].incoming[0].origin].name, "Rock");
}
