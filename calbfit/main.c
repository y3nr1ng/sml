#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calb.h"

/*-------------------------------------------------------------------------
 *  Get the input parameters.
 */
static void
inputs(int argc, char **argv, int *verb, char *fnX, char *fnY, char *outfn)
{
    int  idx;

    if (argc < 4) {
	printf("Usage: %s [-v] <calb_rawX_file> <calb_rawY_file> <outfn>\n",
		argv[0]);
	exit(0);
    }
    if (argc == 5 && strcmp(argv[1], "-v") == 0) {
	*verb = 1;
	idx   = 2;
    }
    else {
	*verb = 0;
	idx   = 1;
    }
    strncpy(fnX,   argv[idx],   BUFLEN-1);
    strncpy(fnY,   argv[idx+1], BUFLEN-1);
    strncpy(outfn, argv[idx+2], BUFLEN-1);
}

/*-------------------------------------------------------------------------
 *  Main program
 */
int main(int argc, char **argv)
{
    int     verb;
    char    fnX[BUFLEN], fnY[BUFLEN], outfn[BUFLEN];
    calb_t *cbX, *cbY;

    inputs(argc, argv, &verb, fnX, fnY, outfn);
    cbX = calb_readraw(fnX);
    cbY = calb_readraw(fnY);

    calb_fit(verb, cbX, "Wx");
    calb_fit(verb, cbY, "Wy");
    calb_output(outfn, cbX, cbY);

    return 0;
}
