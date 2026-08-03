/* Wrapper header for bindgen - includes all primer3 public headers.
 * libprimer3.h includes thal.h, oligotm.h, dpal.h, p3_seq_lib.h, masker.h
 * in the correct order. */
#include "vendor/primer3/src/libprimer3.h"

/* create_empty_seq_lib is defined in p3_seq_lib.c but not declared in
 * the public header. We declare it here for bindgen. */
seq_lib *create_empty_seq_lib(void);

