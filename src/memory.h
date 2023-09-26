#ifndef _MEMORY_H
#define _MEMORY_H
#include <stdlib.h> 

int init_workarea(size_t size);
void finalize_workarea(void);
void *get_workarea(char *name, int n, size_t unit);
int remove_workarea(char *name);
void print_workarea_usage(void);
int exceed_ulim(char *name, void *p);

#endif
