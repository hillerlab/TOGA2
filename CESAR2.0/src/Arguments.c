/**
 * Argument parser implementation.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#include <stdbool.h>
#include <string.h>
#include <stdio.h>

#include "Logging.h"
#include "Params.h"

#include "Arguments.h"

#ifndef VERSION
#define VERSION "0.01 build"
#endif

/**
 * Print information about the program.
 */
void print_version() {
  printf("Cesar %s\n"
         "Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede\n"
         "\n"
         "Align a pair of sequences using HMM.\n",
         VERSION);
}

/**
 * Print a list of available input parameters.
 */
void print_help() {
  print_version();
  printf("\n"
         "Usage: cesar <fasta file>\n"
         "             -c/--clade <human|mouse>\n"
         "             -m/--matrix <matrix file>\n"
         "             -p/--profiles <acc profile> <do profile>\n"
         "             -i/--split_codon_emissions <acc split codon emissions> <do split codon emissions>\n"
         "             -f/--firstexon\n"
         "             -l/--lastexon\n"
         "             -x/--max-memory\n"
         "             -v/--verbosity <0 .. 2>\n"
         "             -V/--version\n"
         "             -s/--set name1=value1 .. nameN=valueN\n"
         "             -S/--sanityChecks  (error-exit if sanity checks fail)\n"
         "             -h/--help\n");
}

/**
 * Read arguments from user input and store them in a Params struct.
 * @param argc number of given arguments.
 * @param argv pointer to argument vector
 * @param parameters a Params struct.
 * @return success boolean.
 */
bool Arguments__read(int argc, char** argv, struct Params* parameters) {
  if (argc < 2) {
    print_help();
    die("Insufficient number of arguments. Please provide at least an input file.");
  }

  uint16_t num_input_files = 0;
  int i = 1;  // [0] is the cesar binary itself.
  while (argv[i] && i < argc) {
    char* argument = argv[i];

    if (!strcmp(argument, "-h") || !strcmp(argument, "-?") || !strcmp(argument, "--help")) {
      print_help();
      return false;
    } else if (!strcmp(argument, "-V") || !strcmp(argument, "--version")) {
      print_version();
      return false;
    }

    if (!strcmp(argument, "-i") || !strcmp(argument, "--split_codon_emissions")) {
      Params__set_via_str(parameters, "split_emissions_acceptor", argv[++i]);
      Params__set_via_str(parameters, "split_emissions_donor", argv[++i]);
    } else if (!strcmp(argument, "-x") || !strcmp(argument, "--max-memory")) {
      Params__set_via_str(parameters, "max_memory", argv[++i]);
    } else if (!strcmp(argument, "-c") || !strcmp(argument, "--clade")) {
      Params__set_via_str(parameters, "clade", argv[++i]);
    } else if (!strcmp(argument, "-m")  || !strcmp(argument, "--matrix")) {
      Params__set_via_str(parameters, "eth_file", argv[++i]);
    } else if (!strcmp(argument, "-v")  || !strcmp(argument, "--verbosity")) {
      g_loglevel = (char) atoi(argv[++i]);
    } else if (!strcmp(argument, "-p")  || !strcmp(argument, "--profiles")) {
      Params__set_via_str(parameters, "acc_profile", argv[++i]);
      Params__set_via_str(parameters, "do_profile", argv[++i]);
		parameters->acc_do_specified = true;
    } else if (!strcmp(argument, "-s")  || !strcmp(argument, "--set")) {
      int j = i+1;
      for (; j < argc; j++) {
        if (argv[j][0] == '-') {
          break;
        }
        char* key = strtok(argv[j], "=");
        char* val = strtok(NULL, "=");
        if(!Params__set_via_str(parameters, key, val)) {
          warn("Ignoring unknown parameter `%s'.", key);
        }
      }
      i = j;
    } else if (!strcmp(argument, "-d") || !strcmp(argument, "--dot")) {
      Params__set_via_str(parameters, "dot", argv[++i]);
    } else if (!strcmp(argument, "-f") || !strcmp(argument, "--firstexon")) {
      parameters->firstexon = true;
    } else if (!strcmp(argument, "-l") || !strcmp(argument, "--lastexon")) {
      parameters->lastexon = true;
    } else if (!strcmp(argument, "-S") || !strcmp(argument, "--sanityChecks")) {
      parameters->sanityChecks = true;
    } else {
      num_input_files++;
      if (num_input_files > 1) {
        die("Too many input files given, e.g.: %s (%u)", argv[i], num_input_files);
      }
      Params__set_via_str(parameters, "fasta_file", argv[i]);
    }

    i++;
  }  // while
 if (num_input_files == 0) {
   die("Too few input files given: %u", num_input_files);
 }
 return true;
}
