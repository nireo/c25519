#ifndef C25519_SEED_H
#define C25519_SEED_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* C25519_SEED_H */
