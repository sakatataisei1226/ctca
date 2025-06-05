#include <stdlib.h>
#include "oh_part.h"

static int psort_comp(const void* x, const void* y);
static int psort_spec(const void* x, const void* y);

void particle_sort_(struct S_particle *p, int *n) {
  qsort(p, *n, sizeof(struct S_particle), psort_comp);
}

void particle_spec_(struct S_particle *p, int *n) {
  qsort(p, *n, sizeof(struct S_particle), psort_spec);
}

static int psort_comp(const void* pa, const void* pb) {
  struct S_particle *a = (struct S_particle*)pa;
  struct S_particle *b = (struct S_particle*)pb;

  return (a->z - b->z);
}

static int psort_spec(const void* pa, const void* pb) {
  struct S_particle *a = (struct S_particle*)pa;
  struct S_particle *b = (struct S_particle*)pb;

  return (a->spec - b->spec);
}
