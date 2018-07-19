/*-------------------------------------------------------------------------
 *
 *	Plot and generate JPEG image for spots
 *
 *------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "pspot.h"

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
 *	Parse the arguments of the input file.
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
    while (s >= buf && isspace((int)(*s))) {
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
        if (sscanf(buf, "%d", iarg+i) != 1)
            pstop("!!! inp(int): reading input line %d failed.\n", nline);
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
        if (sscanf(buf, "%lf", darg+i) != 1)
            pstop("!!! inp(dbl): reading input line %d failed.\n", nline);
        if ((s = strchr(buf, ' ')) != NULL)
            buf = s;
        while (*buf == ' ' && *buf != '\0')
            buf ++;
    }
}

void inp_getSTR(char *buf, char **pstr, int nline)
{
    *pstr = strdup(buf);
    if (*pstr == NULL)
        pstop("!!! inp(str): reading input line %d failed.\n", nline);
}

/*-------------------------------------------------------------------------
 *
 *	Input parameters
 *
 *------------------------------------------------------------------------*/

void inputs(int argc, char **argv, para_t *p)
{
    FILE   *f;
    int     i=1, nline=0, nlab=0, iarg[256];
    double  darg[256];
    char   *inpf, buf[1024];

    memset(p, 0, sizeof(para_t));
    if (argc < 2) {
	printf("Usage: %s [-v] <input>\n", argv[0]);
	printf("   -v: display the image on the screen.\n");
	exit(0);
    }
    if (strcmp(argv[1], "-v") == 0) {
	p->view = 1;
	i ++;
    }
    inpf = argv[i];

    if ((f = fopen(inpf, "rt")) == NULL)
        pstop("!!! inputs: cannot open file: %s\n", inpf);
    while (fgets(buf, 1023, f) != NULL) {
        nline ++;
        if (inp_cutcmt(buf) == 0) continue;
        nlab  ++;

        switch (nlab) {
	case 1:  inp_getINT(buf, iarg, 1, nline);
		 p->datafmt = iarg[0];
		 break;

        case 2:  inp_getSTR(buf, &(p->inpfn), nline);
                 break;

        case 3:  inp_getSTR(buf, &(p->inpfnF), nline);
                 break;

	case 4:  inp_getSTR(buf, &(p->outmesh), nline);
		 break;

        case 5:  inp_getDBL(buf, darg, 1, nline);
                 p->m_density = darg[0];
                 break;

	case 6:  inp_getINT(buf, iarg, 1, nline);
		 p->imgscale = iarg[0];
		 break;

	case 7:  inp_getINT(buf, iarg, 2, nline);
		 p->imgWidth  = iarg[0];
		 p->imgHeight = iarg[1];
		 break;

	case 8:  inp_getDBL(buf, darg, 2, nline);
		 p->iW = darg[0];
		 p->iH = darg[1];
		 break;

	case 9:  inp_getINT(buf, iarg, 1, nline);
		 p->n_intvl = iarg[0];
		 break;

	case 10: inp_getDBL(buf, darg, 1, nline);
		 p->mesh_size = darg[0];
		 break;

	case 11: inp_getDBL(buf, darg, 1, nline);
		 p->gamma = darg[0];
		 break;

	case 12: inp_getSTR(buf, &(p->fontfn), nline);
		 break;
	}
    }
    fclose(f);
    if (strcmp(p->fontfn, "NULL") == 0) {
	char *s = getenv("PSPOT_FONT");
	if (s != NULL) {
	    free(p->fontfn);
	    if ((p->fontfn = strdup(s)) == NULL)
        	pstop("!!! inputs: not enough memory.\n");
	}
    }
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

    inputs(argc, argv, &p);
    if (p.datafmt == 0) {
	readfsts(&p);
	readspot(&p);
	spot_image(&p);
	if (p.view == 0)
	    output_JPEG(&p, 0, p.inpfn);
    }
    else if (p.datafmt == 1) {
	readfsts(&p);
	readspot(&p);
	mesh_events(&p);
	mesh_output(&p);
	gray_image(&p);
	output_JPEG(&p, 1, p.outmesh);
    }
    else {
	readpixel(&p);
	gray_image(&p);
	output_JPEG(&p, 1, p.inpfn);
    }
    if (p.view == 1)
	display_image(argc, argv, &p);
    return 0;
}
