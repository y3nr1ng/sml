#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "aimg.h"

/*-------------------------------------------------------------------------
 *
 *      Stop code with error messages.
 *
 *------------------------------------------------------------------------*/

void pstop(char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    exit(1);
}

/*-------------------------------------------------------------------------
 *
 *      Parse the arguments of the input file.
 *
 *------------------------------------------------------------------------*/

int skipline(FILE *f)
{
    int c;

    do {
        c = fgetc(f);
    } while (c != EOF && c != '\n');

    return (c == EOF) ? -1 : 0;
}

int inp_cutcmt(char *buf)
{
    char *s;

    if ((s = strrchr(buf, '!')) != NULL) {
        *s = '\0';
        s --;
    }
    s = (s == NULL) ? buf+strlen(buf)-1 : s;
    while (s >= buf && isspace(*s)) {
        *s = '\0';
        s --;
    }
    return (s >= buf) ? 1 : 0;
}

void inp_getINT(char *buf, int *iarg, int n, int nline)
{
    int   i;
    char *s;

    for (i=0; i < n; i++) {
        if (sscanf(buf, "%d", iarg+i) != 1) {
            printf("!!! inp(int): reading input line %d failed.\n", nline);
            exit(1);
        }
        if ((s = strchr(buf, ' ')) != NULL)
            buf = s;
        while (*buf == ' ' && *buf != '\0')
            buf ++;
    }
}

void inp_getUINT(char *buf, unsigned int *iarg, int n, int nline)
{
    int   i;
    char *s;

    for (i=0; i < n; i++) {
        if (sscanf(buf, "%u", iarg+i) != 1) {
            printf("!!! inp(unsigned int): reading input line %d failed.\n",
		    nline);
            exit(1);
        }
        if ((s = strchr(buf, ' ')) != NULL)
            buf = s;
        while (*buf == ' ' && *buf != '\0')
            buf ++;
    }
}

void inp_getDBL(char *buf, double *darg, int n, int nline)
{
    int   i;
    char *s;

    for (i=0; i < n; i++) {
        if (sscanf(buf, "%lf", darg+i) != 1) {
            printf("!!! inp(dbl): reading input line %d failed.\n", nline);
            exit(1);
        }
        if ((s = strchr(buf, ' ')) != NULL)
            buf = s;
        while (*buf == ' ' && *buf != '\0')
            buf ++;
    }
}

void inp_getSTR(char *buf, char **pstr, int nline)
{
    *pstr = strdup(buf);
    if (*pstr == NULL) {
        printf("!!! inp(str): reading input line %d failed.\n", nline);
        exit(1);
    }
}

/*-------------------------------------------------------------------------
 *
 *      Get parameters from input file.
 *
 *------------------------------------------------------------------------*/

static void inputs(int argc, char **argv, para_t *p)
{
    FILE   *f;
    int     nline=0, nlab=0, iarg[256];
    double  darg[256];
    char   *inpf, buf[1024];

    if (argc != 2) {
        printf("Usage: %s <inputf>\n", argv[0]);
        exit(0);
    }
    inpf = argv[1];

    memset(p, 0, sizeof(para_t));
    if ((f = fopen(inpf, "rt")) == NULL)
        pstop("!!! inputs: cannot open file: %s\n", inpf);
    while (fgets(buf, 1023, f) != NULL) {
        nline ++;
        if (inp_cutcmt(buf) == 0) continue;
        nlab  ++;

        switch (nlab) {
        case 1:  inp_getSTR(buf, &(p->outf), nline);
                 break;

	case 2:  inp_getUINT(buf, &(p->seed), 1, nline);
		 break;

	case 3:  inp_getINT(buf, iarg, 2, nline);
		 p->iW = iarg[0];
		 p->iH = iarg[1];
		 break;

	case 4:  inp_getDBL(buf, darg, 2, nline);
		 p->I1 = darg[0];
		 p->I2 = darg[1];
		 break;

	case 5:  inp_getDBL(buf, darg, 2, nline);
		 p->N1 = darg[0];
		 p->N2 = darg[1];
		 break;

	case 6:  inp_getDBL(buf, darg, 2, nline);
		 p->w1 = darg[0];
		 p->w2 = darg[1];
		 break;

	case 7:  inp_getINT(buf, &(p->np), 1, nline);
		 break;

	case 8:  inp_getINT(buf, &(p->spixel), 1, nline);
		 break;

	case 9:  inp_getINT(buf, &(p->smesh), 1, nline);
		 break;
	}
    }
    fclose(f);
}

/*-------------------------------------------------------------------------
 *
 *      Main program.
 *
 *------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
    para_t p;

    inputs(argc, argv, &p);
    srand(p.seed);
    aimg(&p);
    output_JPEG(&p);
    output_spot(&p);
    output_raw(&p);

    return 0;
}
