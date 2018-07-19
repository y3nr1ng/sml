#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aimg.h"

/*-------------------------------------------------------------------------
 *
 *	Output the list of spots.
 *
 *------------------------------------------------------------------------*/

void output_spot(para_t *p)
{
    char  fn2[1024], *s;
    FILE *f;
    int   i;

    strncpy(fn2, p->outf, 1023);
    if ((s = strrchr(fn2, '.')) != NULL)
	*s = '\0';
    strncat(fn2, ".tab", 1023);

    if ((f = fopen(fn2, "wt")) == NULL)
	pstop("output_spots: cannot open file: %s\n", fn2);
    fprintf(f, "%4s%11s%11s%11s%10s%13s\n", "idx", "x", "y", "w", "II", "S/N");
    fprintf(f, "------------------------------------------------------------\n");
    for (i=0; i < p->np; i++) {
	fprintf(f, "%04d%11.6f%11.6f%11.6f%10.6f%13.6f\n",
		   i, p->sp[i].x, p->sp[i].y, p->sp[i].w, p->sp[i].II,
		   p->sp[i].II/p->noise);
    }
    fclose(f);
}

/*-------------------------------------------------------------------------
 *
 *	Output the raw values of the image.
 *
 *------------------------------------------------------------------------*/

void output_raw(para_t *p)
{
    char  fn2[1024], *s;
    int   x, y, xx;
    FILE *f;

    strncpy(fn2, p->outf, 1023);
    if ((s = strrchr(fn2, '.')) != NULL)
	*s = '\0';
    strncat(fn2, ".txt", 1023);

    if ((f = fopen(fn2, "wt")) == NULL)
	pstop("output_raw: cannot open file: %s\n", fn2);
    fprintf(f, "%5s%5s%6s\n", "x", "y", "pixel");
    fprintf(f, "----------------\n");

    xx = 0;
    for (y=0; y < p->iH; y++) {
    for (x=0; x < p->iW; x++) {
	fprintf(f, "%5d%5d%6d\n", x, y, (int)(p->img[xx]));
	xx ++;
    }}
    fclose(f);
}
