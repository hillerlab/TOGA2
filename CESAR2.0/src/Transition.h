/**
 * Transition definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#ifndef TRANSITION_H_
#define TRANSITION_H_

#include "Logodd.h"
#include "Stateid.h"

typedef struct Transition {
  STATE_ID_T origin;
  LOGODD_T logodd;
} Transition;

#endif  // TRANSITION_H_
