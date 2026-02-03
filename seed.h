#ifndef C25519_SEED_H
#define C25519_SEED_H

#include <stddef.h>
#include <stdio.h>

#if defined(__linux__) && !defined(C25519_NO_GETRANDOM)
#include <errno.h>
#include <sys/random.h>
#include <sys/types.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

static inline int ed25519_create_seed(unsigned char* seed)
{
    size_t offset = 0;
#if defined(__linux__) && !defined(C25519_NO_GETRANDOM)
    while (offset < 32) {
        ssize_t read_len = getrandom(seed + offset, 32U - offset, 0);
        if (read_len > 0) {
            offset += (size_t)read_len;
            continue;
        }
        if (read_len == -1 && errno == EINTR) {
            continue;
        }
        if (read_len == -1 && (errno == ENOSYS || errno == EPERM)) {
            break;
        }
        return -1;
    }
    if (offset == 32) {
        return 0;
    }
#endif
    FILE* f = fopen("/dev/urandom", "rb");
    if (!f) {
        return -1;
    }
    size_t read_len = fread(seed + offset, 1, 32 - offset, f);
    fclose(f);

    return read_len == 32 - offset ? 0 : -1;
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* C25519_SEED_H */
