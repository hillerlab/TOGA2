/**
 * State implementation.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Logging.h"
#include "EmissionTable.h"
#include "SafeAlloc.h"

#include "State.h"

/**
 * Initialize a State.
 * @param self a state.
 * @param name the state name.
 * @param num_emissions the number of literals the state will emit.
 * @param emission_reference an array of the reference sequence.
 * @param emission_table a pointer to an emission table.
 * @return success boolean.
 */
bool State__init(struct State* self, char name[STATE_NAME_LENGTH], uint16_t num_emissions, Literal emission_reference[num_emissions], struct EmissionTable* emission_table) {
  logv(5, "State init:\t%s", name);
  strncpy(self->name, name, STATE_NAME_LENGTH-1);
  self->name[STATE_NAME_LENGTH-1] = '\0';

  if (g_loglevel >= 5) {
    char tmp[4] = "";
    Literal__str(num_emissions, emission_reference, tmp);
    logv(5, SID"\tref=%s\t%s(%i)", self->id, tmp, name, num_emissions);
  }

  self->num_incoming = 0;
  self->num_emissions = num_emissions;
  for (uint16_t i=0; i < num_emissions; i++) {
    self->reference[i] = emission_reference[i];
  }
  self->emission_table = emission_table;

  // Die if emission table does not allow the state to emit its reference.
  // This originally was thought to avoid stop codon emissions.
  /*
  if (num_emissions) {
    uint16_t reference_id = Literal__uint(num_emissions, emission_reference);
    if (!EmissionTable__emittable(emission_table, reference_id)) {
      char tmp[4] = "";
      Literal__str(num_emissions, emission_reference, tmp);
      die("State cannot emit reference %s: %s ", tmp, self->name);
    }
  }
  */

  return true;
}

/**
 * Fewer arguments for irrelevant reference.
 * @param self a state.
 * @param name the name of a state.
 * @param num_emissions number of emissitted literals from this state.
 * @param emission_table the table containing num_emissions entries.
 * @return success boolean.
 */
bool State__init_uniform(struct State* self, char name[STATE_NAME_LENGTH], uint16_t num_emissions, struct EmissionTable* emission_table) {
  Literal emission_reference[STATE_MAX_REF_LEN] = { LITERAL_A, LITERAL_A, LITERAL_A };
  return State__init(self, name, num_emissions, emission_reference, emission_table);
}

/**
 * Silent State
 * @param self a state.
 * @param name a string.
 * @return success boolean.
 */
bool State__init_silent(struct State* self, char name[STATE_NAME_LENGTH]) {
  return State__init(self, name, 0, NULL, NULL);
}

/**
 * Add an incoming transition to the target state.
 * @param self the origin of the transition.
 * @param target the target state of the incoming transition.
 * @param logodd the log of the probability to transide from self to target.
 * @return success boolean.
 */
bool State__add_incoming(struct State* self, LOGODD_T logodd, struct State* target) {
  if(target->num_incoming > STATE_MAX_NUM_INCOMING) {
    die("Too many incoming transitions for state %s ("SID").", target->name, target->id);
  }
  logv(5, "State incoming:\t%s--(%f)-->%s", self->name, Logodd__exp(logodd), target->name);
  struct Transition t = {.origin = self->id, .logodd = logodd};
  target->incoming[target->num_incoming++] = t;
  return true;
}

/**
 * Compose a string representation of the state.
 * @param self the state.
 * @param buffer a string that can be filled.
 * @return success boolean.
 */
bool State__str(struct State* self, char* buffer) {
  char* tmp = SAFEMALLOC(sizeof(char) * self->num_emissions+1);
  Literal__str(self->num_emissions, self->reference, tmp);
  sprintf(buffer, "State {.name=\"%s\", .id=%Lu, .num_emissions=%u, .reference=\"%s\", .num_incoming=%u}", self->name, (long long unsigned int) self->id, (unsigned int) self->num_emissions, tmp, (unsigned int) self->num_incoming);
  free(tmp);
  return true;
}
