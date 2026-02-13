#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "../sha512.h"

static int hex_val(char c) {
  if (c >= '0' && c <= '9') {
    return c - '0';
  }
  if (c >= 'a' && c <= 'f') {
    return 10 + (c - 'a');
  }
  if (c >= 'A' && c <= 'F') {
    return 10 + (c - 'A');
  }
  return -1;
}

static int hex_to_bytes(uint8_t *out, size_t out_len, const char *hex) {
  size_t hex_len = strlen(hex);
  if (hex_len != out_len * 2) {
    return 0;
  }
  for (size_t i = 0; i < out_len; i++) {
    int hi = hex_val(hex[i * 2]);
    int lo = hex_val(hex[i * 2 + 1]);
    if (hi < 0 || lo < 0) {
      return 0;
    }
    out[i] = (uint8_t)((hi << 4) | lo);
  }
  return 1;
}

static int expect_digest(const uint8_t *msg, size_t msg_len,
                         const char *expected_hex) {
  uint8_t expected[SHA512_DIGEST_SIZE];
  if (!hex_to_bytes(expected, sizeof(expected), expected_hex)) {
    fprintf(stderr, "invalid expected hex\n");
    return 0;
  }

  uint8_t digest[SHA512_DIGEST_SIZE];
  sha512_sum(msg, msg_len, digest);
  if (memcmp(digest, expected, sizeof(digest)) != 0) {
    fprintf(stderr, "digest mismatch\n");
    return 0;
  }
  return 1;
}

static int test_sha512_abc(void) {
  const char *msg = "abc";
  const char *expected = "ddaf35a193617abacc417349ae20413112e6fa4e89a97ea2"
                         "0a9eeee64b55d39a2192992a274fc1a836ba3c23a3feebbd"
                         "454d4423643ce80e2a9ac94fa54ca49f";
  return expect_digest((const uint8_t *)msg, strlen(msg), expected);
}

static int test_sha512_empty(void) {
  const char *expected = "cf83e1357eefb8bdf1542850d66d8007d620e4050b5715dc"
                         "83f4a921d36ce9ce47d0d13c5d85f2b0ff8318d2877eec2f"
                         "63b931bd47417a81a538327af927da3e";
  return expect_digest((const uint8_t *)"", 0, expected);
}

static int test_sha512_abcdbc(void) {
  const char *msg = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
  const char *expected = "204a8fc6dda82f0a0ced7beb8e08a41657c16ef468b228a8"
                         "279be331a703c33596fd15c13b1b07f9aa1d3bea57789ca0"
                         "31ad85c7a71dd70354ec631238ca3445";
  return expect_digest((const uint8_t *)msg, strlen(msg), expected);
}

static int test_sha512_long(void) {
  const char *msg =
      "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmnhijklmno"
      "ijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu";
  const char *expected = "8e959b75dae313da8cf4f72814fc143f8f7779c6eb9f7fa1"
                         "7299aeadb6889018501d289e4900f7e4331b99dec4b5433a"
                         "c7d329eeb6dd26545e96e55b874be909";
  return expect_digest((const uint8_t *)msg, strlen(msg), expected);
}

static int test_sha512_million_a(void) {
  const char *expected = "e718483d0ce769644e2e42c7bc15b4638e1f98b13b204428"
                         "5632a803afa973ebde0ff244877ea60a4cb0432ce577c31b"
                         "eb009c5c2c49aa2e4eadb217ad8cc09b";

  uint8_t expected_bytes[SHA512_DIGEST_SIZE];
  if (!hex_to_bytes(expected_bytes, sizeof(expected_bytes), expected)) {
    fprintf(stderr, "invalid expected hex\n");
    return 0;
  }

  sha512_ctx ctx;
  sha512_init(&ctx);

  uint8_t buf[1000];
  memset(buf, 'a', sizeof(buf));
  for (size_t i = 0; i < 1000; i++) {
    sha512_update(&ctx, buf, sizeof(buf));
  }

  uint8_t digest[SHA512_DIGEST_SIZE];
  sha512_final(&ctx, digest);
  if (memcmp(digest, expected_bytes, sizeof(digest)) != 0) {
    fprintf(stderr, "digest mismatch\n");
    return 0;
  }
  return 1;
}

int main(void) {
  struct {
    const char *name;
    int (*fn)(void);
  } tests[] = {
      {"sha512_abc", test_sha512_abc},
      {"sha512_empty", test_sha512_empty},
      {"sha512_abcdbc", test_sha512_abcdbc},
      {"sha512_long", test_sha512_long},
      {"sha512_million_a", test_sha512_million_a},
  };

  int failures = 0;
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
    if (!tests[i].fn()) {
      fprintf(stderr, "[FAIL] %s\n", tests[i].name);
      failures++;
    } else {
      fprintf(stdout, "[OK]   %s\n", tests[i].name);
    }
  }

  if (failures) {
    fprintf(stderr, "%d test(s) failed\n", failures);
    return 1;
  }

  printf("All tests passed.\n");
  return 0;
}
