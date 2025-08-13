/**
 * LogoddLogoddMatrix implementation.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "Stateid.h"
#include "Logodd.h"
#include "Logging.h"
#include "SafeAlloc.h"

#include "Matrix.h"

/**
 * Create a PathMatrix.
 * @param columns number of columns.
 * @param rows number of rows.
 * @param default_value value, the matrix will be filled with ab initio.
 * @return pointer to the created matrix.
 */
struct PathMatrix* PathMatrix__create(size_t columns, size_t rows, STATE_ID_T default_value) {
  logv(2, "create(%lux%lu, "SID")", columns, rows, default_value);
  PathMatrix* self = (struct PathMatrix*) SAFEMALLOC(sizeof(struct PathMatrix));

  self->num_rows = rows;
  self->num_columns = columns;

  self->v = (STATE_ID_T*) SAFEMALLOC(sizeof(STATE_ID_T) * rows * columns);
  for (size_t i=0; i < rows * columns; i++) {
    self->v[i] = default_value;
  }

  return self;
}

/**
 * Destroy a PathMatrix.
 * @param self a PathMatrix.
 * @return success boolean.
 */
bool PathMatrix__destroy(struct PathMatrix* self) {
  free(self->v);
  free(self);
  return true;
}
 
/**
 * Set a value in the matrix.
 * @param self the PathMatrix.
 * @param column the first value of the coordinates.
 * @param row the second value of the coordinates.
 * @param value the value that will be written.
 * @return success boolean.
 */
bool PathMatrix__set(struct PathMatrix* self, size_t column, size_t row, STATE_ID_T value) {
  self->v[self->num_rows * column + row] = value;
  return true;
}

/**
 * Get a value in the matrix.
 * @param self the PathMatrix.
 * @param column the first value of the coordinates.
 * @param row the second value of the coordinates.
 * @return the value.
 */
STATE_ID_T PathMatrix__get(struct PathMatrix* self, size_t column, size_t row) {
  if (column >= self->num_columns || row >= self->num_rows) {
    die("Invalid matrix access: %lux%lu[%lu][%lu]", self->num_columns, self->num_rows, column, row);
  }

  return self->v[self->num_rows * column + row];
}

/**
 * Compose a string representation for the PathMatrix.
 * @param self a PathMatrix.
 * @param buffer a pre-allocated string that will contain the resulting representation.
 * @return success boolean.
 */
bool PathMatrix__str(struct PathMatrix* self, char* buffer) {
  buffer[0] = '\0';
  char tmp[255] = "";
  for (size_t row=0; row < self->num_rows; row++) {
    sprintf(tmp, "%lu\t", row);
    strcat(buffer, tmp);
    /*
    if (row == 0 && self->num_rows > 1) {
      strcat(buffer, "/ ");
    } else if (row == self->num_rows-1 && self->num_rows > 1) {
      strcat(buffer, "\\ ");
    } else {
      strcat(buffer, "[ ");
    }
    */

    for (size_t column=0; column < self->num_columns; column++) {
      STATE_ID_T id = PathMatrix__get(self, column, row);
      if (id == STATE_MAX_ID) {
        sprintf(tmp, "NA\t");
      } else {
        sprintf(tmp, SID"\t", id);
      }
      strcat(buffer, tmp);
    }

    /*
    if (row == 0 && self->num_rows > 1) {
      strcat(buffer, "\\");
    } else if (row == self->num_rows-1 && self->num_rows > 1) {
      strcat(buffer, "/");
    } else {
      strcat(buffer, "]");
    }
    */
    strcat(buffer, "\n");
  }

  return true;
}

/**
 * Get the size of a PathMatrix in Bytes.
 * @param self a PathMatrix.
 * @return the number of Bytes used by this PathMatrix
 */
size_t PathMatrix__bytes(struct PathMatrix* self) {
  return (self->num_rows * self->num_columns * sizeof(STATE_ID_T));
}


/*****************************************/


/**
 * Create a LogoddMatrix.
 * @param columns number of columns.
 * @param rows number of rows.
 * @param default_value the resulting matrix will be filled with this value.
 * @return pointer to resulting matrix.
 */
struct LogoddMatrix* LogoddMatrix__create(size_t columns, size_t rows, LOGODD_T default_value) {
  logv(2, "create(%lux%lu, %E)", (unsigned long) columns, (unsigned long) rows, default_value);
  LogoddMatrix* self = (struct LogoddMatrix*) SAFEMALLOC(sizeof(struct LogoddMatrix));

  self->num_rows = rows;
  self->num_columns = columns;

  self->v = (LOGODD_T*) SAFEMALLOC(sizeof(LOGODD_T) * rows * columns);
  for (size_t i=0; i < rows * columns; i++) {
    self->v[i] = default_value;
  }

  return self;
}

/**
 * Destroy a LogoddMatrix.
 * @param self a LogoddMatrix.
 * @return success boolean.
 */
bool LogoddMatrix__destroy(LogoddMatrix* self) {
  free(self->v);
  free(self);
  return true;
}
 
/**
 * Set a value in the LogoddMatrix
 * @param self the LogoddMatrix.
 * @param column the first value of the coordinates.
 * @param row the second value of the coordinates.
 * @param value the value that will be written.
 * @return success boolean.
 */
bool LogoddMatrix__set(LogoddMatrix* self, size_t column, size_t row, LOGODD_T value) {
	self->v[self->num_rows * column + row] = value;
  return true;
}

/**
 * Get a value from the LogoddMatrix
 * @param self the LogoddMatrix.
 * @param column the first value of the coordinates.
 * @param row the second value of the coordinates.
 * @return the value.
 */
LOGODD_T LogoddMatrix__get(LogoddMatrix* self, size_t column, size_t row) {
  if (column >= self->num_columns || row >= self->num_rows) {
    die("Invalid matrix access: %lux%lu[%lu][%lu]", self->num_columns, self->num_rows, column, row);
  }

  return self->v[self->num_rows * column + row];
}

/**
 * Write string representation of a LogoddMatrix to buffer.
 * @param self a LogoddMatrix.
 * @param buffer pre-allocated string that will contain the resulting representation.
 * @return success boolean.
 */
bool LogoddMatrix__str(LogoddMatrix* self, char* buffer) {
  buffer[0] = '\0';
  char tmp[255] = "";

  if (self->num_columns * self->num_rows > 64*64) {
    warn("You tried to print a very large table.");
  }

  for (size_t row=0; row < self->num_rows; row++) {
    sprintf(tmp, "%lu\t", row);
    strcat(buffer, tmp);
    /*
    if (row == 0 && self->num_rows > 1) {
      strcat(buffer, "/ ");
    } else if (row == self->num_rows-1 && self->num_rows > 1) {
      strcat(buffer, "\\ ");
    } else {
      strcat(buffer, "[ ");
    }
    */

    for (size_t column=0; column < self->num_columns; column++) {
      LOGODD_T logodd = LogoddMatrix__get(self, column, row);
      if (logodd == LOGODD_NEGINF) {
        sprintf(tmp, "-inf\t");
      } else {
        sprintf(tmp, "%+E\t", logodd);
      }
      strcat(buffer, tmp);
    }

    /*
    if (row == 0 && self->num_rows > 1) {
      strcat(buffer, "\\");
    } else if (row == self->num_rows-1 && self->num_rows > 1) {
      strcat(buffer, "/");
    } else {
      strcat(buffer, "]");
    }
    */
    strcat(buffer, "\n");
  }

  return true;
}

/**
 * Get the size of a LogoddMatrix in Bytes.
 * @param self a LogoddMatrix.
 * @return number of Bytes used by given LogoddMatrix.
 */
size_t LogoddMatrix__bytes(struct LogoddMatrix* self) {
  return (self->num_rows * self->num_columns * sizeof(LOGODD_T));
}
