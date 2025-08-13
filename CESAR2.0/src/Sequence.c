/**
 * Sequence management
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "Literal.h"
#include "Logging.h"
#include "SafeAlloc.h"

#include "Sequence.h"

bool Sequence__init(struct Sequence* self) {
  self->length = 0;
  self->length_reserved = 1000;
  self->sequence = SAFECALLOC(sizeof(Literal), self->length_reserved);
  self->name[0] = '\0';
  self->start_split_length = 0;
  self->end_split_length = 0;
  self->codons_offset = 0;
  self->num_codons = 0;
  self->num_codon_bases = 0;
  self->num_states = 0;
  self->num_align_spaces = 0;
  self->hasStartCodonAsSplitCodonState = 0;

  self->acceptor[0] = '\0';
  self->donor[0] = '\0';

  logv(2, "init: length reserved: %lu", self->length_reserved);
  return true;
}

bool Sequence__destroy(struct Sequence* self) {
  free(self->sequence);
  return true;
}

bool Sequence__set_profiles(struct Sequence* self, char acceptor[], char donor[]) {
  strcpy(self->acceptor, acceptor);
  strcpy(self->donor, donor);
  return true;
}

bool Sequence__append(struct Sequence* self, Literal literal) {
  // increase memory space if necessary
  if (self->length+1 >= self->length_reserved) {
    self->length_reserved *= 2;
    logv(2, "length reserved: %lu", self->length_reserved);
    self->sequence = SAFEREALLOC(self->sequence, sizeof(Literal) * self->length_reserved);
    if (self->sequence == NULL) {
      return false;
    }
  }
  self->sequence[self->length++] = literal;
  return true;
}
