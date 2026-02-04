# c25519

A portable implementation of Ed25519 in C, somewhat based on the Go standard library implementation.

Files
- `fe.h`: Field element operations for Curve25519 (5x51-bit limbs).
- `fiat_scalar.h`: Fiat-crypto generated scalar field arithmetic (mod l) in Montgomery form.
- `point.h`: Edwards25519 point types and group operations.
- `scalar.h`: Scalar wrapper helpers (encoding/decoding, arithmetic, NAF/radix16 helpers).
- `scalarmult.h`: Scalar multiplication routines and precomputed tables.
- `tables.h`: Lookup table structures and selection helpers for scalar mult.
- `tests/`: C test programs plus a Makefile to build/run them.

## API overview

All public functions are defined in `c25519.h`. Most functions return `0` on success and `-1` on error.

- Types:
  - `ed25519_public_key` (32 bytes)
  - `ed25519_private_key` (64 bytes, seed || public key)
  - `ed25519_seed` (32 bytes)
  - `ed25519_signature` (64 bytes, R || S)
- Key generation:
  - `ed25519_keypair` creates a keypair using the built-in RNG.
  - `ed25519_keypair_from_seed` derives a keypair from a 32-byte seed.
- Signing:
  - `ed25519_sign` signs a message using a private key.
  - `ed25519_sign_from_seed` signs directly from a 32-byte seed.
  - `ed25519_sign_ph` and `ed25519_sign_ph_from_seed` sign a 64-byte prehash.
- Verification:
  - `ed25519_verify` and `ed25519_verify_ph` verify signatures with a public key.

## Usage

```c
#include <stdio.h>
#include <string.h>

#include "c25519.h"

int main(void)
{
    ed25519_public_key pk;
    ed25519_private_key sk;
    if (ed25519_keypair(&pk, &sk) != 0) {
        fprintf(stderr, "keypair failed\n");
        return 1;
    }

    const char* msg = "hello";
    ed25519_signature sig;
    if (ed25519_sign(sig, (const uint8_t*)msg, strlen(msg), &sk) != 0) {
        fprintf(stderr, "sign failed\n");
        return 1;
    }

    if (ed25519_verify(sig, (const uint8_t*)msg, strlen(msg), &pk) != 0) {
        fprintf(stderr, "verify failed\n");
        return 1;
    }

    printf("signature OK\n");
    return 0;
}
```

## Building tests

```
make -C tests run
```

## Disabling fiat-crypto inline assembly

The fiat-crypto scalar code uses a tiny inline-asm value barrier on GCC/Clang. If you need to disable all inline assembly (e.g., for certain toolchains, sanitizers, or restricted build environments), define:

```
-DFIAT_25519_SCALAR_NO_ASM
```

This makes the value barrier a no-op while keeping the same scalar arithmetic implementation.

## Seeded APIs and Linux entropy

You can provide a 32-byte seed directly to avoid relying on the built-in RNG:

- `ed25519_public_key_from_seed`
- `ed25519_keypair_from_seed`
- `ed25519_sign_from_seed`
- `ed25519_sign_ph_from_seed`

On Linux, `ed25519_create_seed` prefers `getrandom(2)` and falls back to `/dev/urandom`. To force the fallback path, define:

```
-DC25519_NO_GETRANDOM
```

## Benchmarks

```
seed generation: 13us (75539 per second)
key generation: 12us (80628 per second)
message signing (short message): 13us (75883 per second)
message verifying (short message): 40us (24901 per second)
scalar addition: 0.01us (183003002 per second)
key exchange: 38us (26130 per second)
```
