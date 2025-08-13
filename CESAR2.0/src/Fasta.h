/**
 * Fasta parser declarations
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#ifndef FASTA_H_
#define FASTA_H_

#include <stdbool.h>

#include "Params.h"
#include "Sequence.h"

typedef struct Fasta {
  uint16_t num_references;
  uint16_t num_queries;
  struct Sequence** references;
  struct Sequence** queries;
} Fasta;

bool Fasta__init(struct Fasta* self);
bool Fasta__destroy(struct Fasta* self);
bool Fasta__read(struct Fasta* self, char* filename);
bool Fasta__add_reference(struct Fasta* self, struct Sequence* reference);
bool Fasta__add_query(struct Fasta* self, struct Sequence* query);

#endif  // FASTA_H_
