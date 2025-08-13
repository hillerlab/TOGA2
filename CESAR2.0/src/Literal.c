/**
 * Literal definition.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */
#include "Logging.h"

#include "Literal.h"

/**
 * Given a character, return its Literal representation.
 * An invalid character causes this function to kill the process (die).
 * @param c a character.
 * @return the according Literal.
 */
Literal Literal__from_char(char c) {
  switch (c) {
    case 'a':
    case 'A':
      return LITERAL_A;
    case 'c':
    case 'C':
      return LITERAL_C;
    case 'g':
    case 'G':
      return LITERAL_G;
    case 't':
    case 'T':
      return LITERAL_T;
    case 'n':
    case 'N':
      return LITERAL_N;
    default:
      die("Unkown literal 0x%x=%c", c, c);
  }
  return LITERAL_N;
}

/**
 * Given a Literal, return its character representation.
 * This function kills the process if the literal is unknown (die).
 * @param literal a Literal.
 * @return character representation.
 */
char Literal__char(Literal literal) {
  switch (literal) {
    case LITERAL_A:
      return 'A';
    case LITERAL_C:
      return 'C';
    case LITERAL_G:
      return 'G';
    case LITERAL_T:
      return 'T';
    case LITERAL_N:
      return 'N';
    default:
      die("Unknown literal 0x%x", literal);
  }
  return 'N';
}

/**
 * Given an array of literals, fill a string with character representations.
 * @param length the number of Literals.
 * @param literals the array of Literals.
 * @param result the string buffer that will contain the string representation.
 */
void Literal__str(size_t length, Literal literals[length], char result[length+1]) {
  for (size_t i=0; i < length; i++) {
    result[i] = Literal__char(literals[i]);
  }
  result[length] = '\0';
}


/**
 * Convert an array of literals to an unsigned int.
 * @param length size of the array.
 * @param array the array of literals.
 * @return an uint8_fastest.
 */
EMISSION_ID_T Literal__uint(uint16_t length, Literal array[length]) {
  EMISSION_ID_T byte = 0;

  for (uint16_t i=0; i < length; i++) {
    if (array[i] == LITERAL_N) {
      warn("N literal found.");
    }
    byte <<= 2;
    byte += array[i];
  }

  return byte;
}


/**
 * Count the number of N-Literals in a given array.
 * @param length length of the given array.
 * @param array an array of Literals.
 * @return the number of N-Literals.
 */
uint16_t Literal__Ns(uint16_t length, Literal array[length]) {
  uint16_t count = 0;
  for (uint16_t i=0; i < length; i++) {
    count += array[i] == LITERAL_N;
  }
  return count;
}
