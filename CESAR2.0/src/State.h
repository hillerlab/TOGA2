/**
 * State definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#ifndef STATE_H_
#define STATE_H_

#include <stdint.h>
#include <stdbool.h>

#include "Literal.h"
#include "Transition.h"
#include "Logodd.h"
#include "Stateid.h"

#define STATE_NAME_LENGTH 20
#define STATE_MAX_REF_LEN 3
#define STATE_MAX_NUM_INCOMING 15
#define STATE_CUSTOM_T uint16_t

typedef struct State {
  STATE_ID_T id;
  char name[STATE_NAME_LENGTH];
  struct EmissionTable* emission_table;
  uint16_t num_emissions;
  Literal reference[STATE_MAX_REF_LEN];
  uint16_t num_incoming;
  struct Transition incoming[STATE_MAX_NUM_INCOMING];
  STATE_CUSTOM_T custom;
} State;

bool State__init(struct State* self, char name[STATE_NAME_LENGTH], uint16_t num_emissions, Literal emission_reference[num_emissions], struct EmissionTable* emission_table);
bool State__init_uniform(struct State* self, char name[STATE_NAME_LENGTH], uint16_t num_emissions, struct EmissionTable* emission_table);
bool State__init_silent(struct State* self, char name[STATE_NAME_LENGTH]);
bool State__add_incoming(struct State* self, LOGODD_T logodd, struct State* target);
bool State__str(struct State* self, char* buffer);
LOGODD_T State__remainder(struct State* self);

#endif  // STATE_H_
