/**
 * Params defintion.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#ifndef PARAM_CONTAINER_H_
#define PARAM_CONTAINER_H_

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "Logodd.h"
#include "Literal.h"
#include "Profile.h"
#include "EmissionTable.h"
#include "Sequence.h"

#define PATH_STRING_LENGTH 1024
#define NUM_CS_DELETIONS 10

typedef struct Params {
  char clade[PATH_STRING_LENGTH],
    //blosum_file[PATH_STRING_LENGTH],
    eth_file[PATH_STRING_LENGTH],
    acc_profile[PATH_STRING_LENGTH],
    do_profile[PATH_STRING_LENGTH],
    first_codon_profile[PATH_STRING_LENGTH],
    last_codon_profile[PATH_STRING_LENGTH],
    u12_acc_profile[PATH_STRING_LENGTH],
    u12_donor_profile[PATH_STRING_LENGTH],
    fasta_file[PATH_STRING_LENGTH],
    dot[PATH_STRING_LENGTH];

  struct EmissionTable* emission_table_64_UNIFORM,
    * emission_table_16_UNIFORM,
    * emission_table_4_UNIFORM,
    * emission_table_61_UNIFORM,
    * emission_table_61_LAMBDA,
    * emission_table_64_LAMBDA,
    * emission_table_start_LAMBDA,
    * emission_table_stop_LAMBDA;

  bool dirty;
  
  bool acc_do_specified;  /* set to 1 if -p {acc_profile} {do_profile} is given */

  size_t num_factors;
  long double multiple_cd_factors[10];

  long double no_leading_introns_prob,
           no_trailing_introns_prob;

  size_t split_emissions_acceptor,
          split_emissions_donor,
          max_memory;

  uint16_t num_start_codons;
  Literal* start_codons;

  uint16_t num_stop_codons;
  Literal* stop_codons;

  bool multiexon, lastexon, firstexon, forcelong, sanityChecks;

  LOGODD_T stop_codon_emission_logodd,
           fs_logodd,
           ci_logodd,
           ci2_logodd,
           total_cd_logodd,
           cd_logodd,

           cd_acc,
           cd_do,

           no_leading_introns_logodd,
           no_trailing_introns_logodd,
           intron_del,

           c3_i1_do,
           i3_i1_acc,
           i3_i1_do,
           i3_js_acc,
           i3_js_do,

           js_js,
           bas_sca,
           b1_bas,
           b2_bas,

           b1_acc,
           b1_b2,
           b2_b2,
           b2_acc,
           skip_acc,

           js_scd,
           bsd_do,
           bsd_do_id,
           do2_e1,
           do2_e2,
           skip_do,

           splice_js,

           splice_nti,
           nti_js,
           nti_nti,

           splice_i1,
           i3_i1,
           i3_js,

           js_c1,
           c3_i1,
           c3_js,

           bsd_e2,
           e1_e1,
           e1_e2;
} Params;  // struct params

bool Params__create(struct Params* self, struct EmissionTable emission_tables[6]);
void Params__destroy(struct Params* self);
bool Params__recalculate(struct Params* self);
bool Params__set_paths(struct Params* self, char *BaseDir);
bool Params__set_via_str(struct Params* self, char* string, char* value);

#endif  // PARAM_CONTAINER_H_
