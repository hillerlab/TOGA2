/**
 * Fasta parser implementation.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Logging.h"
#include "SafeAlloc.h"
#include "Literal.h"
#include "Params.h"

#include "Fasta.h"

#ifndef g_loglevel
#define g_loglevel 0
#endif

/**
 * Initialize an empty Fasta struct.
 * @param self a Fasta struct.
 * @return success boolean.
 */
bool Fasta__init(struct Fasta* self) {
  self->num_references = 0;
  self->num_queries = 0;
  self->queries = SAFEMALLOC(0);
  self->references = SAFEMALLOC(0);
  return self->references != NULL && self->queries != NULL;
}


/**
 * Destroy a Fasta struct.
 * @param self a Fasta struct.
 * @return success boolean.
 */
bool Fasta__destroy(struct Fasta* self) {
  for(uint16_t i = 0; i < self->num_references; i++) {
    Sequence__destroy(self->references[i]);
  }
  free(self->references);

  for(uint16_t i = 0; i < self->num_queries; i++) {
    Sequence__destroy(self->queries[i]);
  }
  free(self->queries);

  return true;
}


/**
 * Add a new reference sequence.
 * @param self a Fasta struct.
 * @param reference a reference sequence.
 * @return success boolean.
 */
bool Fasta__add_reference(struct Fasta* self, struct Sequence* reference) {
  self->num_references++;
  self->references = SAFEREALLOC(self->references, sizeof(struct Sequence*) * self->num_references);
  if (self->references == NULL) {
    return false;
  }
  self->references[self->num_references-1] = reference;
  return true;
}


/**
 * Add a new query sequence.
 * @param self a Fasta struct.
 * @param query a query sequence.
 * @return success boolean.
 */
bool Fasta__add_query(struct Fasta* self, struct Sequence* query) {
  self->num_queries++;
  self->queries = SAFEREALLOC(self->queries, sizeof(struct Sequence*) * self->num_queries);
  if (self->queries == NULL) {
    return false;
  }
  self->queries[self->num_queries-1] = query;
  return true;
}


/**
 * Read a Fasta file.
 * @param self a Fasta struct.
 * @param filename 
 * @return success boolean.
 */
bool Fasta__read(struct Fasta* self, char* filename) {
  FILE* file_descriptor = fopen(filename, "r");
  if (file_descriptor == NULL) {
    die("Cannot open file: %s", filename);
  }
  
  struct Sequence* sequence = SAFECALLOC(sizeof(Sequence), 1);
  Sequence__init(sequence);

  uint16_t name_length = 0;
  uint16_t acc_length = 0;
  uint16_t do_length = 0;
  uint16_t state = 0;
  size_t lineno = 0;
  bool reached_codons = false;
  bool reached_queries = false;

  // initialize new sequence
  while (1) {
    char c = fgetc(file_descriptor);

    if (c == '\n') {
      lineno++;
    }

    logv(3, "Fasta-read-in sequence: %p\treferences: %u\tqueries: %u\tline: %lu\tstate: %u\tc: '%c'", sequence, self->num_references, self->num_queries, lineno, state, c);
    switch (state) {
      case 0:  // awaiting '>' for sequence name
        if (c == EOF) {
          state = 3;
          break;
        }

        if (c == '#') {
          reached_queries = true;
          state = 13;
          break;
        }

        if (c == '>') {
          state = 1;
          break;
        }

        break;

      case 1:  // reading sequence name
        if (c == '\t') {
          sequence->name[name_length] = '\0';
          state = 11;
          break;
        }
        if (c == '\n') {
          sequence->name[name_length] = '\0';
          state = 2;
          break;
        }
        sequence->name[name_length++] = c;
        break;

      case 11:  // reading acceptor profile file name
        if (c == '\t') {
          sequence->acceptor[acc_length] = '\0';
          acc_length = 0;
          logv(3, "Full acceptor filename: %s", sequence->acceptor);
          state = 12;
          break;
        }
        sequence->acceptor[acc_length++] = c;
        break;

      case 12:  // reading donor profile file name
        if (c == '\n') {
          sequence->donor[do_length] = '\0';
          do_length = 0;
          logv(3, "Full donor filename: %s", sequence->donor);
          state = 2;
          break;
        }
        sequence->donor[do_length++] = c;
        break;

      case 13:  // between references and query
        if (c == '\n') {
          state = 0;
          break;
        }
        break;

      case 2:  // read sequence
        if (c == '\n') {
          state = 0;

          if (sequence != NULL) {
            // save current sequence
            if (reached_queries) {
              Fasta__add_query(self, sequence);
            } else {
              sequence->num_codons = floor(sequence->num_codon_bases / 3);
              Fasta__add_reference(self, sequence);
            }
          }

          // initialize new sequence
          sequence = SAFECALLOC(sizeof(Sequence), 1);
          Sequence__init(sequence);
          name_length = 0,
          reached_codons = false;

          break;
        }
        // handle aligned fasta
        if (c == ' ') {
          sequence->num_align_spaces++;
          break;  // skip spaces
        }
        if (c == '-') {
          // ignore deletes
          break;
        }
        if (c == '|') {
          sequence->start_split_length += sequence->num_codon_bases;
          sequence->codons_offset = sequence->start_split_length;
          sequence->end_split_length = 0;
          sequence->num_codon_bases = 0;
          reached_codons = false;
          state = 21;
          break;
        }
        if ('a' <= c && c <= 'z') {  // lower case
          // lower case indicate split codon
          if(!reached_codons) {
            sequence->codons_offset++;
          }
          if (sequence->length < 2) {
            sequence->start_split_length++;  // 3 = codon length
          } else {
            sequence->end_split_length++;
          }
        } else if ('A' <= c && c <= 'Z') {  // upper case
          reached_codons = true;
          sequence->num_codon_bases++;
        }
        Sequence__append(sequence, Literal__from_char(c));
        break;

      case 21:  // full codons following pipes
        if (c == '\n') {
          state = 0;

          if (sequence != NULL) {
            // save current sequence
            if (reached_queries) {
              Fasta__add_query(self, sequence);
            } else {
              sequence->num_codons = floor(sequence->num_codon_bases / 3);
              Fasta__add_reference(self, sequence);
            }
          }

          // initialize new sequence
          sequence = SAFECALLOC(sizeof(Sequence), 1);
          Sequence__init(sequence);
          name_length = 0,
          reached_codons = false;

          break;
        }
        if (c == '-') {
          // ignore deletes
          break;
        }
        if (c == '|') {
          state = 22;
          break;
        }
        reached_codons = true;
        sequence->num_codon_bases++;
        Sequence__append(sequence, Literal__from_char(c));
        break;

      case 22:  // split codon donor end
        if (c == '\n') {
          state = 0;

          if (sequence != NULL) {
            // save current sequence
            if (reached_queries) {
              Fasta__add_query(self, sequence);
            } else {
              sequence->num_codons = floor(sequence->num_codon_bases / 3);
              Fasta__add_reference(self, sequence);
            }
          }

          // initialize new sequence
          sequence = SAFECALLOC(sizeof(Sequence), 1);
          Sequence__init(sequence);
          name_length = 0,
          reached_codons = false;

          break;
        }
        if (c == '-') {
          // ignore deletes
          break;
        }
        sequence->end_split_length++;
        Sequence__append(sequence, Literal__from_char(c));
        break;

      case 3:  // stop
        if (sequence != NULL && !reached_queries) {
          logv(1, "No separator line beginning with '#' found. Interpreting last sequence as query sequence.");
          logv(3, "References: %u, Queries: %u", self->num_references, self->num_queries);
          // remove from references, add to queries
          Fasta__add_query(self, self->references[self->num_references-1]);
          self->references[self->num_references-1] = NULL;
          self->num_references--;
          logv(3, "References: %u, Queries: %u", self->num_references, self->num_queries);
        }
        fclose(file_descriptor);
        return true;

      default:
        fclose(file_descriptor);
        die("Unknown state: %u", state);
    }  // switch
  }  // while

  fclose(file_descriptor);
  return false;
}
