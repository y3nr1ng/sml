#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <sys/time.h>
#include "pix.h"

/*-------------------------------------------------------------------------
*
*	Stop the code with error messages.
*
*------------------------------------------------------------------------*/

void pstop(char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    fflush(stdout);
    exit(1);
}

/*-------------------------------------------------------------------------
*
*      Timer.
*
*------------------------------------------------------------------------*/

double get_realtime(void) {
    static double t0;
    struct timeval tt;
    double t1, t2, dt=0.0;

    if (gettimeofday(&tt, NULL) != 0) {
        printf("!!! gettimeofday failed.\n");
        exit(1);
    }
    t1 = (double)tt.tv_sec;
    t2 = (double)tt.tv_usec / 1000000.0;
    if (t0 <= 0.0)
        t0 = t1+t2;
    else
        dt = t1+t2-t0;
    return dt;
}

/*-------------------------------------------------------------------------
*
*	Parse the arguments of the input file.
*
*------------------------------------------------------------------------*/

int skipline(FILE *f) {
    int c;

    do {
        c = fgetc(f);
    } while (c != EOF && c != '\n');

    return (c == EOF) ? -1 : 0;
}

int inp_cutcmt(char *buf) {
    char *s;

    if ((s = strrchr(buf, '!')) != NULL) {
        *s = '\0';
        s--;
    }
    s = (s == NULL) ? buf+strlen(buf)-1 : s;
    while ((unsigned long)s >= (unsigned long)buf && isspace((int)(*s))) {
        *s = '\0';
        s--;
    }
    return (s >= buf) ? 1 : 0;
}

void inp_getINT(char *buf, int *iarg, int n, int nline) {
    int i;
    char *s;

    for (i=0; i < n; i++) {
        if (sscanf(buf, "%d", iarg+i) != 1) {
            printf("!!! inp(int): reading input line %d failed.\n", nline);
            exit(1);
        }
        if ((s = strchr(buf, ' ')) != NULL)
            buf = s;
        while (*buf == ' ' && *buf != '\0')
            buf++;
    }
}

void inp_getDBL(char *buf, double *darg, int n, int nline) {
    int i;
    char *s;

    for (i=0; i < n; i++) {
        if (sscanf(buf, "%lf", darg+i) != 1) {
            printf("!!! inp(dbl): reading input line %d failed.\n", nline);
            exit(1);
        }
        if ((s = strchr(buf, ' ')) != NULL)
            buf = s;
        while (*buf == ' ' && *buf != '\0')
            buf++;
    }
}

void inp_getSTR(char *buf, char **pstr, int nline) {
    *pstr = strdup(buf);
    if (*pstr == NULL) {
        printf("!!! inp(str): reading input line %d failed.\n", nline);
        exit(1);
    }
}

/*-------------------------------------------------------------------------
*
*	Read the calibration parameters for find z-position of spots.
*
*------------------------------------------------------------------------*/

static void read_cabf(para_t *p, char *cabf) {
    FILE   *f;
    int i, r, flag[9];
    char   *s, buf[1024];
    double val1, val2;

    for (i=0; i < 9; i++) flag[i]=0;

    if ((f = fopen(cabf, "rt")) == NULL)
        pstop("!!! read_cabf: cannot open file: %s\n", cabf);
    while (fgets(buf, 1023, f) != NULL) {
        s = buf;
        while (*s != '\0' && *s != '=') s++;
        while (*s != '\0' && *s != '-' && !isdigit((int)(*s))) s++;

        if (strncmp(buf, "w0x", 3) == 0 &&
            sscanf(s, "%lf", &val1) == 1) {
            p->cax.w0 = val1;
            if ((s = strstr(buf, "+- ")) != NULL)
                sscanf(s+2, "%lf", &(p->cax.dw0));
            flag[0] = 1;
        }else if (strncmp(buf, "w0y", 3) == 0 &&
                  sscanf(s, "%lf", &val1) == 1) {
            p->cay.w0 = val1;
            if ((s = strstr(buf, "+- ")) != NULL)
                sscanf(s+2, "%lf", &(p->cay.dw0));
            flag[0] = 1;
        }else if (strncmp(buf, "WxA", 3) == 0) {
            r = sscanf(s, "%lf %lf", &val1, &val2);
            if (r >= 1) {
                p->cax.A[0] = val1;
                p->cax.A[1] = (r == 2) ? val2 : 0.0;
                if ((s = strstr(buf, "+- ")) != NULL)
                    sscanf(s+2, "%lf", &(p->cax.dA));
                flag[1] = 1;
            }
        }else if (strncmp(buf, "WxB", 3) == 0) {
            r = sscanf(s, "%lf %lf", &val1, &val2);
            if (r >= 1) {
                p->cax.B[0] = val1;
                p->cax.B[1] = (r == 2) ? val2 : 0.0;
                if ((s = strstr(buf, "+- ")) != NULL)
                    sscanf(s+2, "%lf", &(p->cax.dB));
                flag[2] = 1;
            }
        }else if (strncmp(buf, "Wxc", 3) == 0) {
            r = sscanf(s, "%lf %lf", &val1, &val2);
            if (r >= 1) {
                p->cax.c[0] = val1;
                p->cax.c[1] = (r == 2) ? val2 : 0.0;
                if ((s = strstr(buf, "+- ")) != NULL)
                    sscanf(s+2, "%lf", &(p->cax.dc));
                flag[3] = 1;
            }
        }else if (strncmp(buf, "Wxd", 3) == 0) {
            r = sscanf(s, "%lf %lf", &val1, &val2);
            if (r >= 1) {
                p->cax.d[0] = val1;
                p->cax.d[1] = (r == 2) ? val2 : 0.0;
                if ((s = strstr(buf, "+- ")) != NULL)
                    sscanf(s+2, "%lf", &(p->cax.dd));
                flag[4] = 1;
            }
        }else if (strncmp(buf, "WyA", 3) == 0) {
            r = sscanf(s, "%lf %lf", &val1, &val2);
            if (r >= 1) {
                p->cay.A[0] = val1;
                p->cay.A[1] = (r == 2) ? val2 : 0.0;
                if ((s = strstr(buf, "+- ")) != NULL)
                    sscanf(s+2, "%lf", &(p->cay.dA));
                flag[5] = 1;
            }
        }else if (strncmp(buf, "WyB", 3) == 0) {
            r = sscanf(s, "%lf %lf", &val1, &val2);
            if (r >= 1) {
                p->cay.B[0] = val1;
                p->cay.B[1] = (r == 2) ? val2 : 0.0;
                if ((s = strstr(buf, "+- ")) != NULL)
                    sscanf(s+2, "%lf", &(p->cay.dB));
                flag[6] = 1;
            }
        }else if (strncmp(buf, "Wyc", 3) == 0) {
            r = sscanf(s, "%lf %lf", &val1, &val2);
            if (r >= 1) {
                p->cay.c[0] = val1;
                p->cay.c[1] = (r == 2) ? val2 : 0.0;
                if ((s = strstr(buf, "+- ")) != NULL)
                    sscanf(s+2, "%lf", &(p->cay.dc));
                flag[7] = 1;
            }
        }else if (strncmp(buf, "Wyd", 3) == 0) {
            r = sscanf(s, "%lf %lf", &val1, &val2);
            if (r >= 1) {
                p->cay.d[0] = val1;
                p->cay.d[1] = (r == 2) ? val2 : 0.0;
                if ((s = strstr(buf, "+- ")) != NULL)
                    sscanf(s+2, "%lf", &(p->cay.dd));
                flag[8] = 1;
            }
        }
    }
    fclose(f);

    for (i=0; i < 9; i++) {
        if (flag[i] == 0)
            pstop("!!! read_cabf: reading value error.\n");
    }
}

/*-------------------------------------------------------------------------
*
*	Get parameters from input file.
*
*------------------------------------------------------------------------*/

static void inputs(int argc, char **argv, para_t *p) {
    FILE   *f;
    int nline=0, nlab=0, iarg[256];
    double darg[256];
    char   *inpf, cabf[1024], buf[1024];

    if (argc != 2) {
        printf("Usage: %s <inputf>\n", argv[0]);
        exit(0);
    }
    inpf = argv[1];

    if ((f = fopen(inpf, "rt")) == NULL)
        pstop("!!! inputs: cannot open file: %s\n", inpf);
    while (fgets(buf, 1023, f) != NULL) {
        nline++;
        if (inp_cutcmt(buf) == 0) continue;
        nlab++;

        switch (nlab) {
        case 1:  inp_getINT(buf, iarg, 1, nline);
            p->imgfmt = iarg[0];
            break;

        case 2:  inp_getSTR(buf, &(p->imgfn), nline);
            break;

        case 3:  strcpy(cabf, buf);
            break;

        case 4:  inp_getSTR(buf, &(p->outfn), nline);
            break;

        case 5:  if (strncmp(buf, "NULL", 4) != 0)
                inp_getSTR(buf, &(p->outfnH), nline);
            break;

        case 6:  if (strncmp(buf, "NULL", 4) != 0)
                inp_getSTR(buf, &(p->outfnF), nline);
            break;

        case 7:  if (strncmp(buf, "NULL", 4) != 0)
                inp_getSTR(buf, &(p->outfnS), nline);
            break;

        case 8:  if (strncmp(buf, "NULL", 4) != 0)
                inp_getSTR(buf, &(p->outfnp), nline);
            break;

        case 9:  inp_getINT(buf, iarg, 1, nline);
            if (iarg[0] < 0 || iarg[0] > 1)
                pstop("!!! inp: unknown mode: %d\n", iarg[0]);
            p->mode = iarg[0];

        case 10: inp_getINT(buf, iarg, 2, nline);
            p->frameID1 = iarg[0];
            p->frameID2 = iarg[1];
            break;

        case 11: inp_getINT(buf, iarg, 4, nline);
            p->frame_x1 = iarg[0];
            p->frame_x2 = iarg[1];
            p->frame_y1 = iarg[2];
            p->frame_y2 = iarg[3];
            break;

        case 12: inp_getINT(buf, iarg, 2, nline);
            p->x_find_pixels = iarg[0];
            p->y_find_pixels = iarg[1];
            break;

        case 13: inp_getDBL(buf, darg, 2, nline);
            p->threshold1 = darg[0];
            p->threshold2 = darg[1];
            break;

        case 14: inp_getINT(buf, iarg, 1, nline);
            p->nfsep = iarg[0];
            break;

        case 15: inp_getDBL(buf, darg, 2, nline);
            p->nm_px_x = darg[0];
            p->nm_px_y = darg[1];
            break;

        case 16: inp_getDBL(buf, darg, 1, nline);
            p->i_photon = darg[0];
            break;

        case 17: inp_getINT(buf, iarg, 1, nline);
            p->rmode = iarg[0];
            break;

        case 18: inp_getINT(buf, iarg, 1, nline);
            p->alg = iarg[0];
            break;

        case 19: inp_getDBL(buf, darg, 3, nline);
            p->max_dx  = darg[0];
            p->max_dy  = darg[1];
            p->max_dwx = darg[2];
            p->max_dwy = p->max_dwx;
            break;

        case 20: inp_getDBL(buf, &(p->min_SN), 1, nline);
            break;

        case 21: inp_getDBL(buf, &(p->max_dI_I), 1, nline);
            break;

        case 22: inp_getDBL(buf, darg, 3, nline);
            p->z1 = darg[0];
            p->z2 = darg[1];
            p->dz = darg[2];
            break;

        case 23: inp_getINT(buf, &(p->verb), 1, nline);
            break;
        }
    }
    fclose(f);
    if (nlab != 23) pstop("!!! Invalid input file.\n");

    p->threshold1 = p->threshold1 / p->i_photon;
    p->threshold2 = p->threshold2 / p->i_photon;
    if (p->mode > 0)
        read_cabf(p, cabf);
}

/*-------------------------------------------------------------------------
*
*	Main program.
*
*------------------------------------------------------------------------*/

int main(int argc, char **argv) {
    para_t p;

    memset(&p, 0, sizeof(para_t));
    inputs(argc, argv, &p);
    get_realtime();

    frameIO_init(&p);
    if (p.frameID2 > p.frameID1) {
        spot_dframe(&p);
        qsort(p.sp1, p.n_sp1, sizeof(sp_t *), spot_cmp);
        qsort(p.sp2, p.n_sp2, sizeof(sp_t *), spot_cmp);
        out_spotlist(p.outfnp, p.n_sp1, p.x_find_pixels, p.sp1);
    } else {
        spot_sframe(&p);
        out_spotlist(p.outfnp, p.n_sp1, p.x_find_pixels, p.sp1);
    }

    if (p.rmode == 1)
        spot_output_img(&p);
    else
        spot_fitting(&p);

    return 0;
}
