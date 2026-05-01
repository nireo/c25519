# c25519

A portable implementation of Ed25519 in C, somewhat based on the Go standard library implementation.

Files
- `fe.h`: Field element operations for Curve25519 (5x51-bit limbs).
- `fiat_scalar.h`: Fiat-crypto generated scalar field arithmetic (mod l) in Montgomery form.
- `point.h`: Edwards25519 point types and group operations.
- `scalar.h`: Scalar wrapper helpers (encoding/decoding, arithmetic, NAF/radix16 helpers).
- `scalarmult.h`: Scalar multiplication routines and precomputed tables.
- `tables.h`: Lookup table structures and selection helpers for scalar mult.
- `tests/`: C test programs used by the top-level make targets.

## API overview

All public functions are defined in `c25519.h`. Most functions return `0` on success and `-1` on error.

### Types

```
ed25519_public_key (32 bytes)
ed25519_private_key (64 bytes, seed || public key)
ed25519_seed (32 bytes)
ed25519_signature (64 bytes, R || S)
```

### Functions

```
ed25519_keypai creates a keypair using the built-in RNG.
ed25519_keypair_from_seed derives a keypair from a 32-byte seed.
ed25519_sign signs a message using a private key.
ed25519_sign_from_seed signs directly from a 32-byte seed.
ed25519_sign_ph and ed25519_sign_ph_from_seed sign a 64-byte prehash.
ed25519_verify and ed25519_verify_ph verify signatures with a public key.
```


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

## Tests 

```
make test
```

## Single-header generation

To generate an amalgamated single-header version of the library:

```
make single-header
```

That writes `c25519_single.h` in the repo root. You can override the output path with `make SINGLE_HEADER=/tmp/c25519_single.h single-header`, or use `python3 scripts/make_single_header.py -` to write to stdout.

## Disabling inline assembly

The fiat-crypto scalar code uses inline assembly. This can be disabled via:

```
-DFIAT_25519_SCALAR_NO_ASM
```

## Seeded APIs and Linux entropy

You can provide a 32-byte seed directly to avoid relying on the built-in RNG:

- `ed25519_public_key_from_seed`
- `ed25519_keypair_from_seed`
- `ed25519_sign_from_seed`
- `ed25519_sign_ph_from_seed`

On Linux, `ed25519_create_seed` creates random bytes from a cryptographically secure random number generator which is implemented using `getrandom()` and on macOS it's implemented using `getentropy()`.

## Benchmarks

Run them with:

```
make benchmark
```

Benchmarks on an Apple M5 MacBook Pro:

```
seed generation: 650ns (1537565 per second)
key generation: 6us (178557 per second)
sha512 (1024 bytes): 1us (713367 per second)
message signing (short message): 6us (162598 per second)
message verifying (short message): 16us (62587 per second)
scalar addition: 2.06ns (485171082 per second)
key exchange: 15us (65194 per second)
```
