# c25519

A portable implementation of Ed25519 in C somewhat based on the Go standard library implementation.

Files
- `fe.h`: Field element operations for Curve25519 (5x51-bit limbs).
- `fiat_scalar.c`: Fiat-crypto generated scalar field arithmetic (mod l) in Montgomery form.
- `point.h`: Edwards25519 point types and group operations.
- `scalar.h`: Scalar wrapper helpers (encoding/decoding, arithmetic, NAF/radix16 helpers).
- `scalarmult.h`: Scalar multiplication routines and precomputed tables.
- `tables.h`: Lookup table structures and selection helpers for scalar mult.
- `tests/`: C test programs plus a Makefile to build/run them.
