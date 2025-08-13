/**
 * HMM definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#ifndef HMM_H_
#define HMM_H_

#include <stdio.h>
#include <stdbool.h>

#include "State.h"
#include "Transition.h"

#define HMM_MAX_STRING_LENGTH 4096

typedef struct HMM {
  struct Transition* starts;
  size_t num_starts;
  size_t max_starts;

  struct Transition* ends;
  size_t num_ends;
  size_t max_ends;

  struct State* states;
  size_t num_states;
  size_t max_states;
} HMM;

struct HMM* HMM__create(size_t states, size_t num_starts, size_t num_ends);
bool HMM__destroy(struct HMM* hmm);
struct State* HMM__new_state(struct HMM* self);
bool HMM__set_start(struct HMM* self, struct Transition transition);
bool HMM__set_end(struct HMM* self, struct Transition transition);
bool HMM__dot(struct HMM* self, FILE* file);
bool HMM__normalize(struct HMM* self, struct State* state);

#endif  // HMM_H_
