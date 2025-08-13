/**
 * Argument parser.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#ifndef ARGUMENTS_H_
#define ARGUMENTS_H_

#include "Params.h"

void print_version();
void print_help();
bool Arguments__read(int argc, char** argv, struct Params* parameters);

#endif  // ARGUMENTS_H_
