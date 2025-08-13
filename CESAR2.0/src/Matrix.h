/**
 * LogoddMatrix definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#ifndef VMATRIX_H_
#define VMATRIX_H_

#include <stddef.h>

#include "Logodd.h"
#include "State.h"


typedef struct PathMatrix {
  size_t num_rows;
  size_t num_columns;
  STATE_ID_T* v;
} PathMatrix;

struct PathMatrix* PathMatrix__create(size_t columns, size_t rows, STATE_ID_T default_value);

bool PathMatrix__destroy(struct PathMatrix* self);

bool PathMatrix__set(struct PathMatrix* self, size_t column, size_t row, STATE_ID_T value);

STATE_ID_T PathMatrix__get(struct PathMatrix* self, size_t column, size_t row);

bool PathMatrix__str(struct PathMatrix* self, char* buffer);

size_t PathMatrix__bytes(struct PathMatrix* self);

/******************************************/

typedef struct LogoddMatrix {
  size_t num_rows;
  size_t num_columns;
  LOGODD_T* v;
} LogoddMatrix;

struct LogoddMatrix* LogoddMatrix__create(size_t columns, size_t rows, LOGODD_T default_value);

bool LogoddMatrix__destroy(struct LogoddMatrix* self);

bool LogoddMatrix__set(struct LogoddMatrix* self, size_t column, size_t row, LOGODD_T value);

LOGODD_T LogoddMatrix__get(struct LogoddMatrix* self, size_t column, size_t row);

bool LogoddMatrix__str(struct LogoddMatrix* self, char buffer[]);

size_t LogoddMatrix__bytes(struct LogoddMatrix* self);


#endif  // VMATRIX_H_
