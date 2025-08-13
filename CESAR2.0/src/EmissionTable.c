/**
 * EmissionTable implementation.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "Logging.h"
#include "SafeAlloc.h"

#include "EmissionTable.h"

/**
 * Read distribution from file.
 * @param self a pointer to the emission table.
 * @param filename the path to the file.
 * @return success boolean.
 */
bool EmissionTable__read(struct EmissionTable* self, char* filename) {
  if (self->distribution != LAMBDA_DISTRIBUTION) {
    die("Cannot load file into non-lambda distributed emission table.");
  }

  // read distribution from file here.
  FILE* file_descriptor = fopen(filename, "r");
  if (file_descriptor == NULL) {
    die("Cannot open file: %s", filename);
  }

  #define LINELENGTH 1024
  #define DELIMITERS " \t"

  bool done = false;
  size_t lineno=0;
  while (!done) {
    // fill the line
    char line[LINELENGTH];
    for (size_t i=0; i < LINELENGTH; i++) {
      assert(i < LINELENGTH);

      char c = fgetc(file_descriptor);

      bool stop = false;
      switch (c) {
        case EOF:
          line[i] = '\0';
          done = true;
          stop = true;
          break;
        case '\n':
          line[i] = '\0';
          stop = true;
          break;
        default:
          line[i] = c;
      }

      if (stop) {
        break;
      }
    }

    // set matrix entries using the tokens
    size_t i=0;
    size_t line_offset = 0;
    char* token = strtok(line, DELIMITERS);
    while (token != NULL) {
      // skip commented lines
      if (token[0] == '#') {
        if (i == 0) {
          line_offset++;
        }
        token = NULL;
        break;
      }

      // check the first token: it should be the codon of which the array index
      // corresponds to the line number.
      if (i==0) {
        Literal* codon = (Literal*) SAFEMALLOC(sizeof(Literal) * self->num_literals);
        for (uint16_t j=0; j < self->num_literals; j++) {
          codon[j] = Literal__from_char(token[j]);
        }

        size_t expected = lineno - line_offset;
        EMISSION_ID_T index = Literal__uint(self->num_literals, codon);
        free(codon);

        if (index != expected) {
          die("Unsupported order of oligomers found in %s:%lu: Expected %lu, got %u (%s)", filename, lineno+1, expected, index, token);
        }
      } else {
        double prob;
        sscanf(token, "%lf", &prob);
        LogoddMatrix__set(self->values, lineno-line_offset, i-1, Logodd__log(prob));
      }

      token = strtok(NULL, DELIMITERS);
      i++;
    }
    lineno++;
  }

  fclose(file_descriptor);

  return true;
}

/**
 * Set individual entries in the EmissionTable.
 * @param self emission table.
 * @param sequence Literal sequence whose matches to all other sequences are set to logodd.
 * @param logodd the logodd value.
 * @return success boolean.
 */
bool EmissionTable__set(struct EmissionTable* self, Literal sequence[], LOGODD_T logodd) {
  bool result = true;
  // for each query, set logodd concerning 
  EMISSION_ID_T row = Literal__uint(self->num_literals, sequence);
  for (uint16_t column=0; column < pow(4, self->num_literals); column++) {
    result &= LogoddMatrix__set(self->values, column, row, logodd);
    if (!result) {
      break;
    }
  }
  return result;
}


/**
 * Set the emission of given literal sequence to zero probability.
 * @param self emission table.
 * @param sequence sequence of literals that will be forbidden.
 */
bool EmissionTable__forbid(struct EmissionTable* self, Literal sequence[]) {
  return EmissionTable__set(self, sequence, LOGODD_NEGINF);
}

/**
 * Recursively check for Ns in query and replace them by A,C,T,G while recording those in visited.
 */
void EmissionTable__variants(uint16_t num_literals, Literal query[num_literals], uint16_t position, bool visited[]) {
  if (g_loglevel >= 6) {
    char tmp[4] = "";
    Literal__str(num_literals, query, tmp);
    logv(7, "Variants? %s[%u/%u]", tmp, position, num_literals);
  }

  if (position >= num_literals) {
    visited[Literal__uint(num_literals, query)] = true;
    return;
  }

  if (query[position] != LITERAL_N) {
    return EmissionTable__variants(num_literals, query, position+1, visited);
  }

  Literal copy[4];
  for (uint16_t i=0; i < num_literals; i++) {
    copy[i] = query[i];
  }

  for (Literal l=LITERAL_A; l < LITERAL_N; l++) {
    copy[position] = l;
    EmissionTable__variants(num_literals, copy, position+1, visited);
  }
}

/**
 * Look up an emission table for two arrays of literals.
 * @param self the emission table.
 * @param reference will be used to look up the row.
 * @param query will be used to look up the column.
 * @return the probability to emit query.
 */
LOGODD_T EmissionTable__by_literals(struct EmissionTable* self, Literal reference[], Literal query[]) {
  uint16_t reference_Ns = Literal__Ns(self->num_literals, reference);
  uint16_t query_Ns = Literal__Ns(self->num_literals, query);

  if (query_Ns == 0 && reference_Ns == 0) {
    EMISSION_ID_T row = Literal__uint(self->num_literals, reference);
    EMISSION_ID_T column = Literal__uint(self->num_literals, query);

    if (g_loglevel >= 6) {
      char ref[5] = "";
      char qry[5] = "";
      Literal__str(self->num_literals, reference, ref);
      Literal__str(self->num_literals, query, qry);
      logv(7, "    by literals: qry=%s=%i x ref=%s=%i", qry, column, ref, row);
    }

    return EmissionTable__get(self, column, row);
  }

  bool* visited_reference = (bool*) SAFECALLOC(sizeof(bool), self->values->num_columns);
  EmissionTable__variants(self->num_literals, reference, 0, visited_reference);

  uint16_t reference_visits = 0;
  double total_sum = 0;
  for (EMISSION_ID_T row = 0; row < self->values->num_rows; row++) {

    if (!visited_reference[row]) {
      continue;
    }
    reference_visits++;

    bool* visited_query = (bool*) SAFECALLOC(sizeof(bool), self->values->num_columns);
    EmissionTable__variants(self->num_literals, query, 0, visited_query);

    uint16_t query_visits = 0;
    double sum = 0;
    for (EMISSION_ID_T column = 0; column < self->values->num_columns; column++) {
      if (!visited_query[column]) {
        continue;
      }

      query_visits++;
      sum += Logodd__exp(EmissionTable__get(self, column, row));
      logv(7, "Visit: %02x", column);
    }
    
    free(visited_query);
    logv(7, "query_visits=%u\tsum=%f", query_visits, sum);

    total_sum += sum / query_visits;
  }

  logv(7, "reference_visits=%u\ttotalsum=%f", reference_visits, total_sum);
  free(visited_reference);

  return Logodd__log(total_sum / reference_visits);
}

/**
 * Initialize an Emissiontable.
 * @param self the EmissionTable.
 * @param num_emissions the number of literals the table should provide for.
 * @param distribution the type of distribution.
 * @return A pointer to the created emission table.
 */
bool EmissionTable__init(struct EmissionTable* self, EMISSION_ID_T num_emissions, Distribution distribution) {
  self->distribution = distribution;
  self->num_literals = num_emissions;

  // the value for each of N literals will be 1/N
  LOGODD_T value;

  switch (distribution) {
    case UNIFORM_DISTRIBUTION:
      value = Logodd__log((LOGODD_T) pow(.25, num_emissions));
      break;
    case LAMBDA_DISTRIBUTION:
      value = LOGODD_NEGINF;
      break;
    default:
      die("Unknown distribution: %x", distribution);
  }
  self->values = LogoddMatrix__create(pow(4, num_emissions), pow(4, num_emissions), value);

  return true;
}

/**
 * Set the emission probability to one n-th for each of n codons.
 * @param self the EmissionTable.
 * @param num_codons the number of codons (n).
 * @param codons an array of a series of codons (e.g. [TAATGATAG] for all stop codons).
 * @return success boolean.
 */
bool EmissionTable__init_single_codons(struct EmissionTable* self, EMISSION_ID_T num_codons, Literal codons[num_codons]) {
  EmissionTable__init(self, 3, LAMBDA_DISTRIBUTION);
  LOGODD_T one_nth = 1.0 - Logodd__log((double) num_codons);

  for(uint16_t i = 0; i < num_codons; i++) {
    // for each query, set logodd concerning 
    EMISSION_ID_T row = Literal__uint(self->num_literals, &codons[3*i]);
    for (uint16_t j = 0; j < num_codons; j++) {
      EMISSION_ID_T column = Literal__uint(self->num_literals, &codons[3*j]);
      if (! LogoddMatrix__set(self->values, column, row, one_nth)) {
        return false;
      }
    }
  }
  return true;
}


/**
 * Destroy an emission table's sub structure.
 * Notice: This only frees memory that has been allocated during init, not the given pointer itself.
 * @param self the EmissionTable to deconstruct.
 * @return success boolean.
 */
bool EmissionTable__destroy(struct EmissionTable* self) {
  LogoddMatrix__destroy(self->values);
  return true;
}

/**
 * Lookup function for the emission table.
 * @param self emission table.
 * @param row table column.
 * @param column table row.
 * @return emission log odd.
 */
LOGODD_T EmissionTable__get(struct EmissionTable* self, EMISSION_ID_T column, EMISSION_ID_T row) {
  return LogoddMatrix__get(self->values, (size_t) column, (size_t) row);
}

/**
 * Compose a string representation of the emission table.
 * @param self emission table.
 * @param buffer storage of the resulting string.
 * @return success boolean.
 */
bool EmissionTable__str(struct EmissionTable* self, char buffer[]) {
  return LogoddMatrix__str(self->values, buffer);
}


bool EmissionTable__emittable(struct EmissionTable* self, EMISSION_ID_T row) {
  bool emittable = false;
  uint16_t column = 0;
  while (column < self->values->num_columns && !emittable) {
    emittable = EmissionTable__get(self, column++, row) != LOGODD_NEGINF;
  }
  if (!emittable) return false;
  return true;
}
