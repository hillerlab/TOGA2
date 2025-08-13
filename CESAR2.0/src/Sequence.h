/**
 * Fasta parser declarations
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <stdint.h>

#include "Literal.h"
#include "Profile.h"

#define SEQUENCENAMELENGTH 255

typedef struct Sequence {
  size_t length;
  Literal* sequence;
  char name[SEQUENCENAMELENGTH],
       acceptor[PROFILE_FILENAME_LENGTH],
       donor[PROFILE_FILENAME_LENGTH];
  uint16_t start_split_length,
          end_split_length;
  size_t num_states,
         codons_offset,
         num_codon_bases,
         num_align_spaces,
         num_codons,
         length_reserved,
         genome_location_start,
         hasStartCodonAsSplitCodonState;
} Sequence;

bool Sequence__init(struct Sequence* self);
bool Sequence__destroy(struct Sequence* self);
bool Sequence__append(struct Sequence* self, Literal literal);

#endif  // SEQUENCE_H_
