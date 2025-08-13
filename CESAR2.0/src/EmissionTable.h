/**
 * EmissionTable definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#ifndef EMISSIONTABLE_H_
#define EMISSIONTABLE_H_

#include <stddef.h>
#include <stdbool.h>

#include "Logodd.h"
#include "Literal.h"
#include "Matrix.h"
#include "Distribution.h"

typedef struct EmissionTable {
  LogoddMatrix* values;
  Distribution distribution;
  uint16_t num_literals;
} EmissionTable;


bool EmissionTable__read(struct EmissionTable* self, char* filename);
LOGODD_T EmissionTable__by_literals(struct EmissionTable* self, Literal reference[], Literal query[]);
bool EmissionTable__init(struct EmissionTable* self, EMISSION_ID_T num_literals, Distribution distribution);
bool EmissionTable__init_single_codons(struct EmissionTable* self, EMISSION_ID_T num_literals, Literal codons[]);
bool EmissionTable__destroy(struct EmissionTable* self);
bool EmissionTable__set(struct EmissionTable* self, Literal sequence[], LOGODD_T logodd);
bool EmissionTable__forbid(struct EmissionTable* self, Literal sequence[]);
LOGODD_T EmissionTable__get(struct EmissionTable* self, uint_fast16_t column, uint_fast16_t row);
bool EmissionTable__str(struct EmissionTable* self, char buffer[]);
bool EmissionTable__emittable(struct EmissionTable* self, EMISSION_ID_T reference);

#endif  // EMISSIONTABLE_H_
