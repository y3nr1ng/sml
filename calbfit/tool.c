#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calb.h"

void *memalloc(size_t n_bytes, void *v0)
{
    void *v;

    if (v0 == NULL) {
	if ((v = malloc(n_bytes)) == NULL) {
	    printf("!!! memalloc: not enough memory.\n");
	    exit(1);
	}
	memset(v, 0, n_bytes);
    }
    else {
	if ((v = realloc(v0, n_bytes)) == NULL) {
	    printf("!!! memalloc: not enough memory.\n");
	    exit(1);
	}
    }

    return v;
}

FILE *openfile(const char *fn, const char *mode)
{
    FILE *f;

    if ((f = fopen(fn, mode)) == NULL) {
	printf("!!! openfile: cannot open file: %s\n", fn);
	exit(1);
    }
    return f;
}
