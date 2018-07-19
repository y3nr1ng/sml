#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "pclst.h"

/*-------------------------------------------------------------------------
 *
 *	Stop the code with error message.
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

static int inp_cutcmt(char *buf)
{
    char *s;

    if ((s = strrchr(buf, '!')) != NULL) {
        *s = '\0';
        s --;
    }
    s = (s == NULL) ? buf+strlen(buf)-1 : s;
    while (s >= buf && isspace((int)(*s))) {
        *s = '\0';
        s --;
    }
    return (s >= buf) ? 1 : 0;
}

static void inp_getINT(char *buf, int *iarg, int n, int nline)
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

static void inp_getDBL(char *buf, double *darg, int n, int nline)
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

static void inp_getSTR(char *buf, char **pstr, int nline)
{
    *pstr = strdup(buf);
    if (*pstr == NULL) {
        printf("!!! inp(str): reading input line %d failed.\n", nline);
        exit(1);
    }
}

/*-------------------------------------------------------------------------
 *
 *	Get the input parameters.
 *
 *------------------------------------------------------------------------*/

static void input(int argc, char **argv, para_t *p)
{
    FILE   *f;
    int     nline=0, nlab=0;
    double  darg[256];
    char   *inpf, buf[1024];

    memset(p, 0, sizeof(para_t));
    if (argc < 2) {
        printf("Usage: %s <input>\n", argv[0]);
        exit(0);
    }
    inpf = argv[1];

    if ((f = fopen(inpf, "rt")) == NULL)
        pstop("!!! inputs: cannot open file: %s\n", inpf);
    while (fgets(buf, 1023, f) != NULL) {
        nline ++;
        if (inp_cutcmt(buf) == 0) continue;
        nlab  ++;

        switch (nlab) {
	case 1:  inp_getSTR(buf, &(p->spotfn), nline);
		 break;

	case 2:  inp_getSTR(buf, &(p->Fstsfn), nline);
		 break;

	case 3:  inp_getDBL(buf, darg, 2, nline);
		 p->bx1 = darg[0];
		 p->by1 = darg[1];
		 break;

	case 4:  inp_getDBL(buf, darg, 2, nline);
		 p->bx2 = darg[0];
		 p->by2 = darg[1];
		 break;

	case 5:  inp_getDBL(buf, &(p->m_density), 1, nline);
		 break;

	case 6:  inp_getDBL(buf, darg, 1, nline);
		 p->spd = darg[0]*darg[0];
		 break;

	case 7:  inp_getINT(buf, &(p->nsp_draw), 1, nline);
		 break;

	case 8:  inp_getINT(buf, &(p->nrec_draw), 1, nline);
		 break;

	case 9:  inp_getINT(buf, &(p->ncls1), 1, nline);
		 break;

	case 10: inp_getINT(buf, &(p->n_intvl), 1, nline);
		 break;

	case 11: inp_getINT(buf, &(p->verb), 1, nline);
		 break;
	}
    }
    fclose(f);
    if (nlab != 11) pstop("!!! input: invalid input file.\n");

    p->Fstart = 0;
}

/*-------------------------------------------------------------------------
 *
 *	Main program
 *
 *------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
    para_t p;

    input(argc, argv, &p);
    readfsts(&p);
    readspot(&p);
    cluster_identify(&p);
    output(&p);
    xcor_output(&p);

    return 0;
}
