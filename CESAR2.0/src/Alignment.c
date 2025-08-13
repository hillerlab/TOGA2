/**
 * Alignment composition
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "Params.h"
#include "Alignment.h"
#include "Logging.h"
#include "SafeAlloc.h"

/**
 * Write the most probable deletion to a string buffer.
 * @param deletion a partial codon match state, e.g. 1nt-Deletion or 2nt-Deletion.
 * @param emission_table the emission table where to pick the probabilties from.
 * @param reference array of Literals.
 * @param query array of Literals.
 * @param result pre-allocated string buffer.
 * @return success boolean.
 */
bool find_best_deletion(struct State* deletion, struct EmissionTable* emission_table, Literal* reference, Literal* query, char* result) {

  double prob=LOGODD_NEGINF, highest_prob=LOGODD_NEGINF;
  uint_fast16_t best_pos=0, pos=0;
  Literal lookup_query[3];

  if (deletion->num_emissions == 2) {
    for(; pos < 3; pos++) {
      for (uint_fast16_t i=0; i < 3; i++) {
        if (pos == i) {
          lookup_query[i] = LITERAL_N;
          continue;
        }
        if (pos > i) {
          lookup_query[i] = query[i];
          continue;
        }
        if (pos < i) {
          lookup_query[i] = query[i-1];
          continue;
        }
      }
      prob = EmissionTable__by_literals(emission_table, reference, lookup_query);
      if (prob >= highest_prob) {
        highest_prob = prob;
        for (uint_fast16_t j=0; j < 3; j++) {
          result[j] = Literal__char(lookup_query[j]);
        }
        result[pos] = '-';
      }
      logv(6, "lookup_query: %c%c%c (%c%c/%c%c)\tprob (highest): %E (%E)\tbest deletion: '%s' at pos=%u\n",
          Literal__char(lookup_query[0]), Literal__char(lookup_query[1]), Literal__char(lookup_query[2]),
          Literal__char(reference[0]), Literal__char(reference[1]),
          Literal__char(query[0]), Literal__char(query[1]),
          prob, highest_prob, result, best_pos
          );
    }
    return true;
  } else if (deletion->num_emissions == 1) {
    strncpy(result, "---", 3);
    for(; pos < 3; pos++) {
      for (uint_fast16_t i=0; i < 3; i++) {
        if (pos == i) {
          lookup_query[i] = query[0];
        } else {
          lookup_query[i] = LITERAL_N;
        }
      }
      prob = EmissionTable__by_literals(deletion->emission_table, reference, lookup_query);
      if (prob >= highest_prob) {
        highest_prob = prob;
        best_pos = pos;
      }
      logv(6, "lookup_query: %c%c%c (%c%c%c/%c)\tprob (highest): %E (%E)\tbest deletion: '%s' at pos=%u\n",
          Literal__char(lookup_query[0]), Literal__char(lookup_query[1]), Literal__char(lookup_query[2]),
          Literal__char(reference[0]),Literal__char(reference[1]),Literal__char(reference[2]),
          Literal__char(query[0]),
          prob, highest_prob, result, best_pos
          );
    }
    result[best_pos] = Literal__char(query[0]);
    return true;
  }
  return false;
}

/**
 * Create a new alignment.
 * @param fasta the Fasta that contains reference and query sequences.
 * @param query_id the query index of which the alignment shall be created.
 * @param params the Params.
 * @param path_length the length of the Viterbi path.
 * @param path the Viterbi path.
 * @return the alignment.
 */
struct Alignment* Alignment__create(struct Fasta* fasta, uint16_t query_id, struct Params* params, size_t path_length, struct State** path) {
  uint16_t reference_id = 0;
  const char lower = 'a' - 'A';
  int numAlignedRefChars = 0;		/* for sanity check that the entire reference seq is contained in the alignment */

  struct Alignment* self = (struct Alignment*) SAFEMALLOC(sizeof(struct Alignment));
  size_t length = fasta->queries[query_id]->length;
  for (uint16_t ref_id=0; ref_id < fasta->num_references; ref_id++) {
    length += fasta->references[ref_id]->length;
  }
  self->reference = (char*) SAFECALLOC(sizeof(char), length+1+20*fasta->num_references);
  self->query = (char*) SAFECALLOC(sizeof(char), length+1+20*fasta->num_references);

  // assemble the alignment
  char* deletion;
  char bases[4] = "";

  uint_fast16_t pending_deletion=0;
  size_t q = 0, r = 0, t = 0;
  for (size_t i=1; i < path_length; i++) {
    uint16_t j=0, emissions = path[i]->num_emissions;

    // deleting 1nt and 2nt will always emit 3 bases/dashes to maintain reading frame intact.
    if (!strncmp("delete_1nt", path[i]->name, 10) || !strncmp("delete_2nt", path[i]->name, 10)) {
      emissions = 3;
    }

    if (emissions == 0) {
      if (i > 0) {
        logv(3, "i=%lu\tt=%lu\tj=%u\tq=%lu\tr=%lu\tref=''\tqry=''\t%s -> %s (%s)", i, t, j, q, r, path[i-1]->name, path[i]->name, "");
      } else {                                                       
        logv(3, "i=%lu\tt=%lu\tj=%u\tq=%lu\tr=%lu\tref=''\tqry=''\t[] -> %s (%s)", i, t, j, q, r, path[i]->name, "");
      }
    }
    for (; j < emissions; j++) {
      if (r >= fasta->references[reference_id]->length && reference_id < fasta->num_references-1) {
        reference_id++;
        logv(3, "reference id: %u/%u", reference_id, fasta->num_references);
        r = 0;
      }
      if (!strncmp("intron", path[i]->name, 6)) {
        self->reference[t+j] = ' ';
        self->query[t+j] = Literal__char(fasta->queries[query_id]->sequence[q++]) + lower;
      } else if (!strncmp("split_codon", path[i]->name, 11)) {
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]) + lower;
        numAlignedRefChars ++;
        self->query[t+j] = Literal__char(fasta->queries[query_id]->sequence[q++]);
      } else if (!strncmp("start_codon", path[i]->name, 11) || !strncmp("stop_codon", path[i]->name, 10) || !strncmp("match_codon", path[i]->name, 11)) {
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]);
        numAlignedRefChars ++;
        self->query[t+j] = Literal__char(fasta->queries[query_id]->sequence[q++]);
      } else if (!strncmp("match", path[i]->name, 5)) {  // donor, acceptor
        self->reference[t+j] = ' ';
        self->query[t+j] = Literal__char(fasta->queries[query_id]->sequence[q++]) + lower;
      } else if (!strncmp("insert", path[i]->name, 6)) {
        self->reference[t+j] = '-';
        self->query[t+j] = Literal__char(fasta->queries[query_id]->sequence[q++]);
      } else if (!strncmp("delete", path[i]->name, 6)) {
        if (pending_deletion == 0) {
          deletion = (char*) SAFECALLOC(sizeof(char), 3);
          find_best_deletion(
              path[i],
              params->emission_table_61_LAMBDA,
              &fasta->references[reference_id]->sequence[r],
              &fasta->queries[query_id]->sequence[q],
              deletion
              );
          pending_deletion = 3;
          q += path[i]->num_emissions;
        }
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]);
        numAlignedRefChars ++;
        self->query[t+j] = deletion[3-pending_deletion--];

        if (pending_deletion == 0) {
          free(deletion);
        }
      }
      Literal__str(path[i]->num_emissions, path[i]->reference, bases);
      if (i > 0) {
        logv(3, "i=%lu\tt=%lu\tj=%u\tq=%lu\tr=%lu\tref='%c'\tqry='%c'\t%s -> %s (%s)", i, t, j, q, r, self->reference[t+j], self->query[t+j], path[i-1]->name, path[i]->name, bases);
      } else {                                             
        logv(3, "i=%lu\tt=%lu\tj=%u\tq=%lu\tr=%lu\tref='%c'\tqry='%c'\t[] -> %s (%s)", i, t, j, q, r, self->reference[t+j], self->query[t+j], path[i]->name, bases);
      }
    }

    if (!strncmp("start_first_codon", path[i-1]->name, 17) && !strncmp("between_match_donor", path[i]->name, 19)) {
		/* special case that should only exist if the stop codon is the only codon */
      logv(3, "add the deletion of the only existing codon");
      for (; j < 3; j++) {
  	     self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]);
        numAlignedRefChars ++;
     	  self->query[t+j] = '-';
      }
    }

    if (!strncmp("between_split_donor", path[i-1]->name, 19) && !strncmp("between_acc_split", path[i]->name, 17)) {
      for (j=0; j < 19; j++) {
        self->reference[t+j] = '>';
        self->query[t+j] = '-';
      }
    }

    if (i > 0 && !strncmp("between_acc_split", path[i-1]->name, 17) &&
        !strncmp("between_split_ins", path[i]->name, 17)) {
      for (; j < fasta->references[reference_id]->start_split_length; j++) {
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]) + lower;
        numAlignedRefChars ++;
        self->query[t+j] = '-';
        logv(3, "i=%lu\tt=%lu\tj=%u\tq=%lu\tr=%lu\tref='%c'\tqry='%c'\t%s -> %s (%s)", i, t, j, q, r, self->reference[t+j], self->query[t+j], path[i-1]->name, path[i]->name, bases);
      }
    }
	 
    if (i > 0 && !strncmp("start_first_codon", path[i-1]->name, 17) && !strncmp("between_split_donor", path[i]->name, 19)) {
      for (; j < fasta->references[reference_id]->end_split_length; j++) {
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]) + lower;
        numAlignedRefChars ++;
        self->query[t+j] = '-';
        logv(3, "i=%lu\tt=%lu\tj=%u\tq=%lu\tr=%lu\tref='%c'\tqry='%c'\t%s -> %s (%s)", i, t, j, q, r, self->reference[t+j], self->query[t+j], path[i-1]->name, path[i]->name, bases);
      }
    }

    //  between_acc_split -> between_match_ins
    if (!strncmp("between_acc_split", path[i-1]->name, 17) &&
        !strncmp("between_match_ins", path[i]->name, 17)) {
      for (; j < 3; j++) {
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]);
        numAlignedRefChars ++;
        self->query[t+j] = '-';
      }
    }

    // end_codon -> between_split_donor
    if (!strncmp("end_codon", path[i-1]->name, 9) &&
        !strncmp("between_match_donor", path[i]->name, 19)) {
      for (; j < 3; j++) {
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]);
        numAlignedRefChars ++;
        self->query[t+j] = '-';
      }
    }



    if (i > 0 && !strncmp("between_split_donor", path[i]->name, 19) && 
        !strncmp("end_codon", path[i-1]->name, 9) && path[i-1]->custom == 
        fasta->references[reference_id]->num_codons-1-fasta->references[reference_id]->hasStartCodonAsSplitCodonState) {
      for (; j < fasta->references[reference_id]->end_split_length; j++) {
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r++]) + lower;
        numAlignedRefChars ++;
        self->query[t+j] = '-';
      }
    }
    if (i > 0 && (!strncmp("start_first_codon", path[i-1]->name, 17) ||
          !strncmp("end_codon", path[i-1]->name, 9)) && !strncmp("end_codon",
          path[i]->name, 9)) {
      // determine the number of deleted codons from the state names
      size_t deleted_codons = path[i]->custom;
      if (!strncmp("end_codon", path[i-1]->name, 9)) {
        deleted_codons -= path[i-1]->custom;
      } else if (!strncmp("start_first_codon", path[i-1]->name, 17)) {
        deleted_codons -= -1;
      }

      char bases[4] = "";
      Literal__str(path[i]->num_emissions, path[i]->reference, bases);
      for (; j < 3*deleted_codons; j++) {
        if (r >= fasta->references[reference_id]->length && reference_id <= fasta->num_references) {
          r = 0;
          reference_id++;
          logv(3, "next reference: %u", reference_id)
        }
        self->reference[t+j] = Literal__char(fasta->references[reference_id]->sequence[r]);
        numAlignedRefChars ++;
        self->query[t+j] = '-';
        logv(3, "i=%lu\tt=%lu\tj=%u\tq=%lu\tr=%lu\tref='%c'\tqry='%c'\t%s -> %s (%s)", i, t, j, q, r, self->reference[t+j], self->query[t+j], path[i-1]->name, path[i]->name, bases);
        r++;
      }
    }
    t += j;
  }

  /* sanity check that the entire reference seq is contained in the alignment */
  int totalRefLen = 0;
  for (uint16_t i=0; i < fasta->num_references; i++) {
      totalRefLen += fasta->references[i]->length;
  }
  if (numAlignedRefChars != totalRefLen) {
		if (params->sanityChecks) {
	     	fprintf(stderr, "ERROR: there are %d bases in the reference but only %d ref bases are in the final alignment.\n", totalRefLen, numAlignedRefChars);
		   exit(-1);
		}
  }

  return self;
}


/**
 * Destroy an alignment.
 * @param self the alignment.
 * @return success boolean.
 */
bool Alignment__destroy(struct Alignment* self) {
  free(self->reference);
  free(self->query);
  free(self);
  return true;
}
