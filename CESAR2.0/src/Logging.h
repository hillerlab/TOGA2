/**
 * Logging macros.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#ifndef LOGGING_H_
#define LOGGING_H_

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

// set loglevel
#ifndef LOGLEVEL
#if DEBUG
#define LOGLEVEL 5
#else
#define LOGLEVEL 0
#endif  // DEBUG
#endif  // LOGLEVEL

#define FAIL_EXIT 1

extern char g_loglevel;

#define warn(format, ...)\
  fprintf(stderr, "WARNING! %s:%d %s():\t" format "\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__);\

/** logv prints to stderr if given level exceeds LOGLEVEL
 * http://stackoverflow.com/questions/1644868/c-define-macro-for-debug-printing
 */
#define logv(level, format, ...) if (level <= g_loglevel)\
  fprintf(stderr, "VERBOSE%i %s:%d %s():\t" format "\n", level, __FILE__, __LINE__, __func__, ##__VA_ARGS__);\

/**
 * Kill the program after printing some information to stderr.
 */
#define die(format, ...)\
  fprintf(stderr, "\x1B[31mCRITICAL %s:%d %s():\t\x1B[0m" format "\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__); exit(FAIL_EXIT);



#endif  // LOGGING_H_
