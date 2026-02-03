CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic
BENCH_CFLAGS ?= $(CFLAGS)

TARGETS = fe_tests point_tests scalar_tests table_tests scalarmult_tests sha512_tests ed25519_tests
BENCH = benchmark
SRC_FE = test_fe.c
SRC_POINT = test_point.c
SRC_SCALAR = test_scalar.c
SRC_TABLE = test_tables.c
SRC_SCALARMULT = test_scalarmult.c
SRC_SHA512 = test_sha512.c
SRC_ED25519 = test_ed25519.c

all: $(TARGETS)

fe_tests: $(SRC_FE) ../fe.h
	$(CC) $(CFLAGS) -I.. $(SRC_FE) -o fe_tests

point_tests: $(SRC_POINT) ../fe.h ../point.h
	$(CC) $(CFLAGS) -I.. $(SRC_POINT) -o point_tests

scalar_tests: $(SRC_SCALAR) ../scalar.h ../fiat_scalar.h
	$(CC) $(CFLAGS) -I.. $(SRC_SCALAR) -o scalar_tests

table_tests: $(SRC_TABLE) ../fe.h ../point.h ../tables.h
	$(CC) $(CFLAGS) -I.. $(SRC_TABLE) -o table_tests

scalarmult_tests: $(SRC_SCALARMULT) ../fe.h ../point.h ../scalar.h ../tables.h ../scalarmult.h
	$(CC) $(CFLAGS) -I.. $(SRC_SCALARMULT) -o scalarmult_tests

sha512_tests: $(SRC_SHA512) ../sha512.h
	$(CC) $(CFLAGS) -I.. $(SRC_SHA512) -o sha512_tests

ed25519_tests: $(SRC_ED25519) ../c25519.h ../sha512.h ../seed.h ../point.h ../scalar.h ../scalarmult.h ../tables.h ../fe.h
	$(CC) $(CFLAGS) -I.. $(SRC_ED25519) -o ed25519_tests

benchmark: benchmark.c ../c25519.h ../sha512.h ../seed.h ../point.h ../scalar.h ../scalarmult.h ../tables.h ../fe.h ../fiat_scalar.h
	$(CC) $(BENCH_CFLAGS) -I.. benchmark.c -o $(BENCH)
	./$(BENCH); status=$$?; rm -f $(BENCH); exit $$status

run: all
	./fe_tests
	./point_tests
	./scalar_tests
	./table_tests
	./scalarmult_tests
	./sha512_tests
	./ed25519_tests
	$(MAKE) clean

clean:
	rm -f $(TARGETS) $(BENCH)
