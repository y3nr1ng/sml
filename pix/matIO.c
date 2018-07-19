#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "pix.h"

#define  miINT8		1
#define  miUINT8	2
#define  miINT16	3
#define  miUINT16	4
#define  miINT32	5
#define  miUINT32	6
#define  miSINGLE	7
#define  miDOUBLE	9
#define  miINT64	12
#define  miUINT64	13
#define  miMATRIX	14
#define  miCOMPRESSED	15
#define  miUTF8		16
#define  miUTF16	17
#define  miUTF32	18

typedef union {			// for byte swap.
    unsigned long  l;
    unsigned int   i[2];
    unsigned short s[4];
    double         d;
    float          f[2];
    char           c[8];
} byte_t;

/*--------------------------------------------------------------------------
 *
 *  The Utility functions for MAT-file reading.
 *
 *-------------------------------------------------------------------------*/

/*
 *  Swap byte ordering for correct endian.
 */
static void
byteswap(int len, char *s)
{
    char c;
    int  i;

    for (i=0; i < len/2; i++) {
	c = s[i];
	s[i] = s[len-1-i];
	s[len-1-i] = c;
    }
}

/*
 *  Read and check the header of MAT-file.
 */
static void
mat_r_header(FILE *f, matfile_t *m)
{
    char   h[128];
    int    i;
    byte_t b;

    if (fread(h, 1, 128, f) != 128)
	pstop("!!! matRead: cannot read header: %s\n", m->fn);
    if (strncmp(h, "MATLAB 5.0 MAT-file", 19) != 0) {
	h[116] = '\0';
	pstop("!!! matRead: unknown header text: %s\n!!!   %s\n", m->fn, h);
    }
    for (i=116; i < 124; i++) {
	if (h[i] != (char)0 && h[i] != ' ')
	    pstop("!!! matRead: header offset not zero: %s\n", m->fn);
    }
    strncpy(m->desc, h, 116);

    b.c[0] = h[124];
    b.c[1] = h[125];
    if (b.s[0] == 0x0100) {
	m->version = b.s[0];
	m->endian  = 0;
    }
    else {
	b.c[0] = h[125];
	b.c[1] = h[124];
	if (b.s[0] == 0x0100) {
	    m->version = b.s[0];
	    m->endian  = 1;
	}
	else
	    pstop("!!! matRead: unknown header version: %s\n", m->fn);
    }
}

/*
 *  Read the data from the file, and perform the byteswap when necessary.
 */
static int
mat_rfile(FILE *f, matfile_t *m, int dsize, int nelem, void *data)
{
    int     i;
    byte_t  b;
    char   *pdata;

    if (fread(data, dsize, nelem, f) != nelem)
        return 1;
    if (m->endian == 0)
	return 0;
    if (dsize <= 1)
	return 0;

    if (dsize > 8)
	pstop("!!! matRead: dsize too large: %s\n", m->fn);
    pdata = (char *)data;
    for (i=0; i < nelem; i++) {
	memcpy(b.c, pdata, dsize);
	byteswap(dsize, b.c);
	memcpy(pdata, b.c, dsize);
	pdata += dsize;
    }
    return 0;
}

/*
 *  Read the data from the buffer, and perform the byteswap when necessary.
 */
static void
mat_rbuf(void *buf, matfile_t *m, int dsize, int nelem, void *data)
{
    int     i;
    byte_t  b;
    char   *pbuf;

    if (m->endian == 0 || dsize <= 1) {
	memcpy(data, buf, dsize*nelem);
	return;
    }
    if (dsize > 8)
	pstop("!!! matRead: dsize too large: %s\n", m->fn);
    pbuf = (char *)buf;
    for (i=0; i < nelem; i++) {
	memcpy(b.c, pbuf, dsize);
	byteswap(dsize, b.c);
	memcpy(data, b.c, dsize);
	pbuf += dsize;
    }
}


/*-------------------------------------------------------------------------
 *
 *  Data processing routings: for miMATRIX.
 *
 *------------------------------------------------------------------------*/

static int
mat_miMATRIX_flags(matfile_t *m, char *data, unsigned int *m_class)
{
    int           mtype;
    unsigned int  mlen, d1, d2;
    unsigned int  flags, class;

    mat_rbuf(data,   m, 4, 1, &mtype);
    mat_rbuf(data+4, m, 4, 1, &mlen);
    if (mtype != miUINT32 && mlen != 8)
	pstop("!!! matRead: invalid miMATRIX flags: %s\n", m->fn);
    mat_rbuf(data+8,  m, 4, 1, &d1);
    mat_rbuf(data+12, m, 4, 1, &d2);

    class = (d1 & 0x00ff);
    flags = (d1 & 0xff00) >> 8;

    if (flags != 0)
	pstop("!!! mat_miMATRIX_flags: un-implemented flag: %s\n", m->fn);
    if (d2 != 0)
	pstop("!!! mat_miMATRIX_flags: un-implemented sparse matrix: %s\n", m->fn);
    *m_class = class;

    return 16;
}

static int
mat_miMATRIX_dimension(matfile_t *m, char *data, int *dim_x, int *dim_y)
{
    int           mtype;
    unsigned int  mlen;

    mat_rbuf(data,   m, 4, 1, &mtype);
    mat_rbuf(data+4, m, 4, 1, &mlen);
    if (mtype != miINT32 && mlen != 8)
	pstop("!!! matRead: invalid miMATRIX dimension: %s\n", m->fn);
    mat_rbuf(data+8,  m, 4, 1, dim_x);
    mat_rbuf(data+12, m, 4, 1, dim_y);

    return 16;
}

static int
mat_miMATRIX_name(matfile_t *m, char *data, char **name)
{
    int           mtype;
    unsigned int  mlen;

    mat_rbuf(data,   m, 4, 1, &mtype);
    mat_rbuf(data+4, m, 4, 1, &mlen);
    if (mtype != miINT8)
	pstop("!!! matRead: invalid miMATRIX name: %s\n", m->fn);
    if ((*name = malloc(mlen+1)) == NULL)
	pstop("!!! matRead: not enough memory for array name: %s\n", m->fn);
    mat_rbuf(data+8,  m, 1, mlen, *name);
    (*name)[mlen] = '\0';
    mlen = (mlen/8 + ((mlen%8 == 0) ? 0 : 1)) * 8;

    return 8+mlen;
}

static void
mat_miMATRIX(matfile_t *m, unsigned int nbyte, void *data, matmx_t *mx)
{
    unsigned int   cnt=0;
    int            i, dtype=0, dsize=0, dlen=0, nelem=0;

    void           *ppp=NULL;
    char           *pcc=NULL;
    short          *pss=NULL;
    int            *pii=NULL;
    long           *pll=NULL;
    float          *pff=NULL;
    double         *pdd=NULL;
    unsigned char  *puc=NULL;
    unsigned short *pus=NULL;
    unsigned int   *pui=NULL;
    unsigned long  *pul=NULL;
    short          *prf=NULL;

    cnt += mat_miMATRIX_flags(m, ((char*)data)+cnt, &(mx->m_class));
    if (cnt >= nbyte) return;
    cnt += mat_miMATRIX_dimension(m, ((char*)data)+cnt, &(mx->dim_x),
				  &(mx->dim_y));
    if (cnt >= nbyte) return;
    cnt += mat_miMATRIX_name(m, ((char*)data)+cnt, &(mx->name));

    while (cnt < nbyte) {
	mat_rbuf(((char*)data)+cnt,   m, 4, 1, &dtype);
	mat_rbuf(((char*)data)+cnt+4, m, 4, 1, &dlen);
	if ((ppp = malloc(dlen)) == NULL)
	    pstop("!!! matRead: cannot allocate array data: %s\n", m->fn);
	switch (dtype) {
	case miINT8:
	    dsize = 1; pcc = (char *)ppp;           break;
	case miUINT8:
	    dsize = 1; puc = (unsigned char *)ppp;  break;
	case miINT16:
	    dsize = 2; pss = (short *)ppp;          break;
	case miUINT16:
	    dsize = 2; pus = (unsigned short *)ppp; break;
	case miINT32:
	    dsize = 4; pii = (int *)ppp;            break;
	case miUINT32:
	    dsize = 4; pui = (unsigned int *)ppp;   break;
	case miSINGLE:
	    dsize = 4; pff = (float *)ppp;          break;
	case miINT64:
	    dsize = 8; pll = (long *)ppp;           break;
	case miUINT64:
	    dsize = 8; pul = (unsigned long *)ppp;  break;
	case miDOUBLE:
	    dsize = 8; pdd = (double *)ppp;         break;
	default:
	    pstop("!!! mat_miMATRIX: un-implemented data type: %s: %d\n",
		         m->fn, dtype);
	}

	nelem    = dlen/dsize;
	mx->type = dtype;
	mat_rbuf(((char *)data)+cnt+8, m, dsize, nelem, ppp);
	cnt += (8+dlen);

	if (dtype != miINT16 && dtype != miUINT16) {
	    if ((prf = malloc(nelem*sizeof(short))) == NULL)
		pstop("!!! matRead: cannot allocate array data: %s\n",
			     m->fn);
	}
	switch (dtype) {
	case miINT8:
	    for (i=0; i<nelem; i++) prf[i] = (short)(pcc[i]); break;
	case miUINT8:
	    for (i=0; i<nelem; i++) prf[i] = (short)(puc[i]); break;
	case miINT16:
	    prf = pss;          break;
	case miUINT16:
	    prf = (short *)pus; break;
	case miINT32:
	    for (i=0; i<nelem; i++) prf[i] = (short)(pii[i]); break;
	case miUINT32:
	    for (i=0; i<nelem; i++) prf[i] = (short)(pui[i]); break;
	case miSINGLE:
	    for (i=0; i<nelem; i++) prf[i] = (short)(pff[i]); break;
	case miINT64:
	    for (i=0; i<nelem; i++) prf[i] = (short)(pll[i]); break;
	case miUINT64:
	    for (i=0; i<nelem; i++) prf[i] = (short)(pul[i]); break;
	case miDOUBLE:
	    for (i=0; i<nelem; i++) prf[i] = (short)(pdd[i]); break;
	}
	mx->data = (void *)prf;
    }
}

/*-------------------------------------------------------------------------
 *
 *  High level routings to read a variable.
 *
 *------------------------------------------------------------------------*/

/*
 *  Process for various data types from a buffer.
 */
static void
mat_data_proc(matfile_t *m, int dtype, unsigned int nbyte, void *data)
{
    switch (dtype) {
    case miMATRIX:
	mat_miMATRIX(m, nbyte, data, &(m->mx));
	break;
    default:
	pstop("!!! matRead: not implemented data type: %s: %d\n", m->fn, dtype);
    }
}

/*
 *  Read and decompress the binary data.
 */
static void
mat_r_decompress(FILE *f, matfile_t *m, unsigned long nbyte)
{
    char          *d, *dd;
    int            cnt=1, r, dtype, dlen;
    unsigned long  nbyte2;

    nbyte2 = nbyte*cnt;
    if ((d  = malloc(nbyte )) == NULL ||
	(dd = malloc(nbyte2)) == NULL)
	pstop("!!! matRead: not enough memory for compressed data: %s\n",
		     m->fn);
    if (mat_rfile(f, m, 1, nbyte, d) != 0)
	pstop("!!! matRead: cannot read compressed data: %s\n", m->fn);
    while (1) {
	r = uncompress((Bytef *)dd, &nbyte2, (const Bytef *)d, nbyte);
	if (r == Z_OK) break;

	if (r != Z_BUF_ERROR)
	    pstop("!!! matRead: data decompress error: %s: %d\n", m->fn, r);
	cnt ++;
	nbyte2 = nbyte*cnt;
	if (cnt > 5)
	    pstop("!!! matRead: cannot decompress data: %s\n", m->fn);
	if ((dd = realloc(dd, nbyte2)) == NULL)
	    pstop("!!! matRead: not enough memory to decompress data: %s\n", m->fn);
    }

    mat_rbuf(dd, m, 4, 1, &dtype);
    mat_rbuf(dd+4, m, 4, 1, &dlen);
    if (dlen+8 != nbyte2)
	pstop("!!! matRead: inconsistant data length for compressed data: %s\n", m->fn);
    mat_data_proc(m, dtype, dlen, dd+8);

    free(d);
    free(dd);
}


/*--------------------------------------------------------------------------
 *
 *  Read the MAT-file (API).
 *
 *--------------------------------------------------------------------------*/

matfile_t *matRead(char *fn)
{
    matfile_t    *m;
    FILE         *f;
    unsigned int  cnt, dtype, nbyte;
    char         *data;

    if ((m = malloc(sizeof(matfile_t))) == NULL)
	pstop("!!! matRead: not enough memory for MAT-file: %s\n", fn);
    memset(m, 0, sizeof(matfile_t));
    if ((f = fopen(fn, "rb")) == NULL)
	pstop("!!! matRead: cannot open file: %s\n", fn);
    m->fn = strdup(fn);
    mat_r_header(f, m);

    cnt = 0;
    while (mat_rfile(f, m, 1, 4, &dtype) == 0) {
	if (mat_rfile(f, m, 1, 4, &nbyte) != 0)
	    pstop("!!! matRead: reading data error: %s: %u\n", fn, cnt);
	switch (dtype) {
	case miCOMPRESSED:
	    mat_r_decompress(f, m, nbyte);
	    break;

	default:
	    if ((data = malloc(nbyte)) == NULL)
		pstop("!!! matRead: not enough memory for %s: dtype=%u\n",
		       m->fn, nbyte);
	    if (fread(data, 1, nbyte, f) != nbyte)
		pstop("!!! matRead: cannot read data for %s: dtype=%u\n",
		       m->fn, nbyte);
	    mat_data_proc(m, dtype, nbyte, data);
	    free(data);
	    break;
	}
	cnt ++;
    }
    fclose(f);

    return m;
}

void matFree(matfile_t *m)
{
    free(m->fn);
    free(m->mx.name);
    free(m->mx.data);
    free(m);
}
