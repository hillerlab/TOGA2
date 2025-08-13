/**
 * Logodd definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#ifndef LOGODD_H_
#define LOGODD_H_

#include <stddef.h>
#include <stdint.h>
#include <float.h>

//#define LOGODD_T long double
//#define LOGODD_NEGINF -LDBL_MAX
#define LOGODD_T double
#define LOGODD_NEGINF -DBL_MAX

LOGODD_T Logodd__log(LOGODD_T input);
LOGODD_T Logodd__exp(LOGODD_T input);
LOGODD_T Logodd__add(LOGODD_T summand, LOGODD_T addend);
LOGODD_T Logodd__sumexp(LOGODD_T a, LOGODD_T b);
LOGODD_T Logodd__subexp(LOGODD_T a, LOGODD_T b);

#endif  // LOGODD_H_
