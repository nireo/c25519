#ifndef __SEED_H__
#define __SEED_H__

#include <stdio.h>

static inline int ed25519_create_seed(unsigned char* seed)
{
    FILE* f = fopen("/dev/urandom", "rb");
    if (!f) {
        return -1;
    }

    fread(seed, 1, 32, f);
    fclose(f);

    return 0;
}

#endif
