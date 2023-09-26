#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memory.h"

/* temporary work area */

static unsigned char *mempool = 0, *current_p;	// uint8_t
static size_t mempoolsiz = 0;

#define NWALIM 64

static struct workarea {
  void *begin, *end;
  size_t size;
  int n;
  char name[BUFSIZ];
} wa[NWALIM];

int wa_next_id;

#define ALIGN_UNIT 16
		// sizeof(void *) = sizeof(double) = 8
static unsigned int
roundm(unsigned int x)
/* returns the smallest multiple of `unit` not less than x */
{
  static unsigned int unit = ALIGN_UNIT;
  return (x + unit - 1) & ~(unit - 1);
}

static unsigned char *
align(unsigned char *p)
{
  return (unsigned char *)( ((intptr_t)p + ALIGN_UNIT - 1)
			    & ~((intptr_t)ALIGN_UNIT - 1) );
}

int init_workarea(size_t size)
{
  int i;

  mempoolsiz = roundm(size);
  i = posix_memalign((void **)&mempool, ALIGN_UNIT, mempoolsiz);
  if (i != 0)
    return i;
  current_p = align(mempool + 1);
  wa_next_id = 0;
  return 0;
}

void finalize_workarea()
{
  if (mempool)
    free(mempool);
}

/* allocate an n x size bytes area of the type of size `unit` */
void *get_workarea(char *name, int n, size_t unit)
{
  void *p;
  int current_id;

  if (wa_next_id >= NWALIM) {
    fprintf(stderr, "cannot alloc workarea. too many areas.\n"); return 0;
  }
  current_id = wa_next_id;
  wa[current_id].size = n * unit;
  if (current_p - mempool + wa[current_id].size > mempoolsiz) {
    fprintf(stderr, "cannot alloc workarea. insufficient pool size.\n");
    return 0;
  }
  wa[current_id].size = roundm(n * unit);
  wa[current_id].n = n;
  p = wa[current_id].begin = current_p;
  current_p += wa[current_id].size;
  wa[current_id].end = current_p;
  strcpy(wa[current_id].name, name);
  wa_next_id = current_id + 1;
  return p;
}

int remove_workarea(char *name)
{
  if (wa_next_id <= 0) return 0;
  if (strcmp(name, wa[wa_next_id-1].name) != 0) {
    fprintf(stderr, "id %s is not the area on the top.\n", name); return -1;
  }
  wa_next_id--;
  current_p -= wa[wa_next_id].size;
  wa[wa_next_id].name[0] = 0;
  return 0;
}

void print_workarea_usage()
{
  struct workarea *p;
  int id;
  
  printf("Workareas:\n");
  p = wa;
  for (id = 0; id < wa_next_id; p++, id++)
    printf("\t%-16s: %16lx -- %16lx %8lud bytes (%d items)\n",
	   p->name, (uintptr_t)(p->begin), (uintptr_t)(p->end), p->size, p->n);
  return;
}

int exceed_ulim(char *name, void *p)
{
  int id;
  for (id = wa_next_id - 1; id >= 0; id--)
    if (strcmp(name, wa[id].name) == 0) break;
  if (id < 0) return -1;
  return (unsigned char *)p - (unsigned char *)wa[id].begin > (unsigned int) wa[id].size;
}
