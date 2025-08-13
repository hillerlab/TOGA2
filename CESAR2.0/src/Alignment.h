/**
 * Alignment type and function definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "Fasta.h"
#include "Params.h"
#include "State.h"
#include "Literal.h"

typedef struct Alignment {
  size_t length;
  char* reference;
  char* query;
} Alignment;

struct Alignment* Alignment__create(struct Fasta* fasta, uint16_t query_id, struct Params* params, size_t path_length, struct State** sequence);
bool Alignment__destroy(struct Alignment* self);

#endif  // ALIGNMENT_H_
