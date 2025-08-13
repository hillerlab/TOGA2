/**
 * State id definition.
 * Copyright 2017 Peter Schwede
 *
 * These lines are separated from State.h to solve dependencies issues.
 */

#include <stddef.h>

#ifndef STATE_ID_H_
#define STATE_ID_H_

//#define STATE_ID_T size_t
//#define STATE_MAX_ID SIZE_MAX

#define STATE_ID_T uint32_t
#define STATE_MAX_ID UINT32_MAX
#define SID "%u" //PRIu32

#endif  // STATE_ID_H_
