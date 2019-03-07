#include <stdio.h>
#include <stdlib.h>

#ifdef DEBUG
#include <openssl/md5.h>

void spotdump(int sxdim, int sydim, int *sp) {
    int x, y, k, xrng, yrng;

    xrng = sxdim/2;
    yrng = sydim/2;

    k = 0;
    for (y=-yrng; y <= yrng; y++) {
        for (x=-xrng; x <= xrng; x++) {
            printf("(x,y,II):  %03d %03d %d\n", x, y, sp[k]);
            k++;
        }
    }
}

void get_md5sum(const void *buf, size_t n)
{
    unsigned char result[MD5_DIGEST_LENGTH];
    unsigned i;

    MD5((unsigned char*)buf, n, result);
    for (i=0; i < MD5_DIGEST_LENGTH; i++)
        printf("%02x", result[i]);
    printf("\n");
}
#endif
