/**
 * Profile implementation.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "EmissionTable.h"
#include "Logging.h"
#include "SafeAlloc.h"

#include "Profile.h"

bool Profile__read(struct Profile* self, char filename[]) {
  strcpy(self->filename, filename);

  FILE* file_descriptor = fopen(filename, "r");
  if (file_descriptor == NULL) {
    die("Cannot open file: %s", filename);
  }

  #define LINELENGTH 1000
  #define DELIMITERS " \t"
  
  bool done = false;
  size_t lineno = 0;
  while (!done) {
    // fill the line
    char line[LINELENGTH] = "\0";
    for (size_t i=0; i < LINELENGTH; i++) {
      assert(i < LINELENGTH);

      char c = fgetc(file_descriptor);

      bool stop = false;
      switch (c) {
        case EOF:
          line[i] = '\0';
          done = true;
          stop = true;
        case '\n':
          line[i] = '\0';
          i=0;
          stop = true;
          break;
        default:
          line[i] = c;
      }

      if (stop) {
        break;
      }
    }

    size_t i=0;
    size_t line_offset = 0;
    Literal keys[4];
    char* token = strtok(line, DELIMITERS);
    struct EmissionTable* etable;
    while (token != NULL) {

      // skip commented lines
      if (token[0] == '#') {
        if (i == 0) {
          line_offset++;
        }
        token = NULL;  // fulfill postcondition
        break;
      }

      // the first line looks like `A\tC\tG\tT\n'
      if (lineno == 0) {
        // In a profile, num_literals of each emission_table will be const 1.
        keys[i] = Literal__from_char(token[0]);

        token = strtok(NULL, DELIMITERS);
        i++;
        continue;
      }


      if (i == 0) {
        etable = Profile__add_emission(self);
      }

      assert(self->length == lineno-line_offset);

      double prob;
      sscanf(token, "%lf", &prob);
      LogoddMatrix__set(etable->values, keys[i], 0, Logodd__log(prob));

      token = strtok(NULL, DELIMITERS);
      i++;
    }

    lineno++;
  }

  fclose(file_descriptor);

  return true;
}

LOGODD_T Profile__by_literals(struct Profile* self, Literal query[]) {
  LOGODD_T result = 0;
  Literal reference = LITERAL_A;  // not important, which base it is.
  for (size_t i=0; i < self->length; i++) {
    LOGODD_T summand = EmissionTable__by_literals(&self->emission_tables[i], &query[i], &reference);
    result = Logodd__add(result, summand);
    if (result == LOGODD_NEGINF) {
      return result;
    }
  }
  return result;
}

struct Profile* Profile__create(char name[STATE_NAME_LENGTH]) {
  struct Profile* self = (struct Profile*) SAFEMALLOC(sizeof(Profile));
  self->emission_tables = NULL;  // will be SAFEREALLOCed when tables are added during read
  strncpy(self->name, name, STATE_NAME_LENGTH);
  self->length=0;
  return self;
}

void Profile__destroy(struct Profile* self) {
  for (size_t i=0; i < self->length; i++) {
    EmissionTable__destroy(&self->emission_tables[i]);
  }
  free(self->emission_tables);
  free(self);
}

struct EmissionTable* Profile__add_emission(struct Profile* self) {
  size_t pos = self->length++;

  self->emission_tables = (struct EmissionTable*) SAFEREALLOC(self->emission_tables, sizeof(struct EmissionTable) * self->length);

  EmissionTable__init(&self->emission_tables[pos], 1, LAMBDA_DISTRIBUTION);
  return &self->emission_tables[pos];
}

bool Profile__str(struct Profile* self, char* buffer) {
  char tmp[255] = "";
  sprintf(tmp, "%s(%u)\n", self->name, self->length);
  strcat(buffer, tmp);
  for (uint16_t t=0; t < self->length; t++) {
    sprintf(tmp, "%i\t", t);
    strcat(buffer, tmp);
    for (Literal l=0; l < LITERAL_N; l++) {
      sprintf(tmp, "%E", EmissionTable__get(&self->emission_tables[t], l, 0));
      strcat(buffer, tmp);
      if (l == LITERAL_T) {
        sprintf(tmp, "\n");
      } else {
        sprintf(tmp, "\t");
      }
      strcat(buffer, tmp);
    }
  }
  return true;
}
