/**
 * HMM implementation.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Logging.h"
#include "SafeAlloc.h"

#include "HMM.h"

/**
 * Create an empty Hidden Markov Model.
 * @param num_states the number of states.
 * @param num_starts the number of possible first states.
 * @param num_ends the number of final states.
 * @return pointer to the new HMM.
 */
struct HMM* HMM__create(size_t num_states, size_t num_starts, size_t num_ends) {
  struct HMM* self = (struct HMM*) SAFEMALLOC(sizeof(struct HMM));

  self->starts = (struct Transition*) SAFEMALLOC(sizeof(Transition) * num_starts);
  self->num_states = 0;
  self->max_states = num_states;

  self->ends = (struct Transition*) SAFEMALLOC(sizeof(Transition) * num_ends);
  self->num_ends = 0;
  self->max_ends = num_ends;

  self->states = (struct State*) SAFECALLOC(num_states, sizeof(struct State));
  self->num_starts = 0;
  self->max_starts = num_states;

  return self;
}

/**
 * Remove HMM from memory.
 * @param self a HMM.
 * @return success boolean.
 */
bool HMM__destroy(struct HMM* self) {
  free(self->starts);
  free(self->ends);
  free(self->states);
  free(self);
  return true;
}

/**
 * Add state.
 * @param self a HMM.
 * @return pointer to the next used state.
 */
struct State* HMM__new_state(struct HMM* self) {
  if (self->num_states >= self->max_states) {
    die("More states requested than reserved: %lu of %lu", self->num_states+1, self->max_states);
  }

  STATE_ID_T id = self->num_states++;
  self->states[id].id = id;

  return &self->states[id];
}

/**
 * Set state to be a start state.
 * @param self a HMM.
 * @param transition the new transition.
 * @return success boolean.
 */
bool HMM__set_start(struct HMM* self, struct Transition transition) {
  if (self->num_starts >= self->max_starts) {
    return false;
  }

  self->starts[self->num_starts++] = transition;
  return true;
}

/**
 * Set end.
 * @param self a HMM.
 * @param transition the new transition.
 * @return success boolean.
 */
bool HMM__set_end(struct HMM* self, struct Transition transition) {
  if (self->num_ends >= self->max_ends) {
    return false;
  }

  self->ends[self->num_ends++] = transition;
  return true;
}

/**
 * Compose a string representation of a HMM.
 * @param self a HMM.
 * @param file a path to a file where the string representation will be stored in.
 * @return success boolean.
 */
bool HMM__dot(struct HMM* self, FILE * file) {
  fprintf(file, "digraph {\n");
  fprintf(file, "  rankdir = \"LR\"\n");
  for (STATE_ID_T i=0; i < self->num_states; i++) {
    struct State* state = &self->states[i];

    for (int j=0; j < state->num_incoming; j++) {
      struct Transition* incoming = &state->incoming[j];
      struct State* origin = &self->states[incoming->origin];

      if (strlen(origin->name) == 0) {
        origin->name[0] = 'x';
      }
      fprintf(file, "  %s_"SID" -> %s_"SID" [label=\"%le\"];\n", origin->name, origin->id, state->name, state->id, Logodd__exp(incoming->logodd));
    }
  }
  for (uint16_t i=0; i < self->num_starts; i++) {
    struct Transition t = self->starts[i];
    struct State* state = &self->states[t.origin];

    fprintf(file, "  START -> %s_"SID" [label=\"%le\"];\n", state->name, state->id, Logodd__exp(t.logodd));
  }
  for (uint16_t i=0; i < self->num_ends; i++) {
    struct Transition t = self->ends[i];
    struct State* state = &self->states[t.origin];

    fprintf(file, "  %s_"SID" -> END [label=\"%le\"];\n", state->name, state->id, Logodd__exp(t.logodd));
  }
  fprintf(file, "}\n");
  return true;
}

/**
 * Normalize transitions
 * @param hmm a HMM.
 * @param left a state of which the  outgoing transitions are to be normalized.
 * @return success boolean.
 */
bool HMM__normalize(struct HMM* hmm, struct State* left) {
  // gather all outgoing edges that leave the left state
  uint16_t num_outgoing = 0;
  struct Transition* outgoing[15];
  struct State* rights[15];

  for (size_t j=0; j < hmm->num_states; j++) {
    struct State* right = &hmm->states[j];

    for (size_t k=0; k < right->num_incoming; k++) {
      struct Transition* transition = &right->incoming[k];

      if (transition->origin == left->id) {
        rights[num_outgoing] = right;
        outgoing[num_outgoing++] = transition;
      }
    }  // O(+deg)
  }  // O(+deg * S)
  for (size_t k=0; k < hmm->num_ends; k++) {
    struct Transition* transition = &hmm->ends[k];
    if (transition->origin == left->id) {
      rights[num_outgoing] = NULL;
      outgoing[num_outgoing++] = transition;
    }
  }

  // perform sumexp(outgoings) (This is how yahmm does it.)
  LOGODD_T total = LOGODD_NEGINF;
  for (uint16_t j=0; j < num_outgoing; j++) {
    struct Transition* transition = outgoing[j];
    total = Logodd__sumexp(total, transition->logodd);
  }

  if (total == 0.0) {
    logv(5, "No normalization\t%s ==%E==> X", left->name, total);
    return true;
  }

  // substract log(sumexp) from all outgoing transition log odds
  for (uint16_t j=0; j < num_outgoing; j++) {
    LOGODD_T old = outgoing[j]->logodd;
    outgoing[j]->logodd = Logodd__add(outgoing[j]->logodd, -total);
    if (rights[j] != NULL) {
      logv(5, "Normalization:\t%s ==(%f to %f)==> %s", left->name, Logodd__exp(old), Logodd__exp(outgoing[j]->logodd), rights[j]->name);
    } else {
      logv(5, "Normalization:\t%s ==(%f to %f)==> END", left->name, Logodd__exp(old), Logodd__exp(outgoing[j]->logodd));
    }
  }
  return true;
}  // O(+deg * S)
