#pragma once

#include "stdint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t state[4];                                   /* state (ABCD) */
    uint32_t count[2];	        /* number of bits, modulo 2^64 (lsb first) */
    unsigned char buffer[64];                         /* input buffer */
} MD5_CTX;

void MD5Init(MD5_CTX *);
void MD5Update(MD5_CTX *, unsigned char *, unsigned int);
void MD5Final(unsigned char[16], MD5_CTX *);

void md5_raw_input(char output[33], unsigned char *buf, unsigned int len);

#ifdef __cplusplus
}
#endif
