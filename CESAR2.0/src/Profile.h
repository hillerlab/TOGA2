/**
 * Profile definition
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#ifndef PROFILE_H_
#define PROFILE_H_

#include "State.h"
#include "EmissionTable.h"

#define PROFILE_FILENAME_LENGTH 256

typedef struct Profile {
    char name[STATE_NAME_LENGTH];
    char filename[PROFILE_FILENAME_LENGTH];
    uint16_t length;
    struct EmissionTable* emission_tables;
} Profile;

struct Profile* Profile__create(char name[STATE_NAME_LENGTH]);
void Profile__destroy();
bool Profile__read(struct Profile* self, char filename[]);
struct EmissionTable* Profile__add_emission(struct Profile* self);
LOGODD_T Profile__by_literals(struct Profile* self, Literal query[]);
bool Profile__str(struct Profile* self, char* buffer);

#endif  // PROFILE_H_
