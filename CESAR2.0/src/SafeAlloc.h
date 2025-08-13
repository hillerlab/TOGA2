/**
 * Makros for safety checks after memory allocations.
 * Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede
 */

#define outofmem(n, line) fprintf(stderr, "CRITICAL %s:%lu %s():\tOut of memory: %lu bytes\n", __FILE__, line, __func__, n);

//https://stackoverflow.com/posts/16298916/revisions
static void* safe_malloc(size_t n, size_t line) {
  void* p = malloc(n);
  if (!p) {
    outofmem(n, line);
    exit(EXIT_FAILURE);
  }
  return p;
}
#define SAFEMALLOC(n) safe_malloc(n, __LINE__)

static void* safe_calloc(size_t s, size_t n, size_t line) {
  void* p = calloc(s, n);
  if (!p) {
    outofmem(n, line);
    exit(EXIT_FAILURE);
  }
  return p;
}
#define SAFECALLOC(s, n) safe_calloc(s, n, __LINE__)

static void* safe_realloc(void* ptr, size_t n, size_t line) {
  void* p = realloc(ptr, n);
  if (!p) {
    outofmem(n, line);
    exit(EXIT_FAILURE);
  }
  return p;
}
#define SAFEREALLOC(ptr, n) safe_realloc(ptr, n, __LINE__)
