/**
 * Model construction.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "Literal.h"
#include "Profile.h"
#include "Fasta.h"
#include "HMM.h"

struct HMM* multi_exon(struct Params* params, struct Fasta* fasta, struct Profile** acceptors, struct Profile** donors);

#endif  // MODEL_H_
