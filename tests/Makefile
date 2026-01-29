CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic

TARGETS = fe_tests point_tests scalar_tests table_tests
SRC_FE = test_fe.c
SRC_POINT = test_point.c
SRC_SCALAR = test_scalar.c
SRC_TABLE = test_tables.c

all: $(TARGETS)

fe_tests: $(SRC_FE) ../fe.h
	$(CC) $(CFLAGS) -I.. $(SRC_FE) -o fe_tests

point_tests: $(SRC_POINT) ../fe.h ../point.h
	$(CC) $(CFLAGS) -I.. $(SRC_POINT) -o point_tests

scalar_tests: $(SRC_SCALAR) ../scalar.h ../fiat_scalar.c
	$(CC) $(CFLAGS) -I.. $(SRC_SCALAR) -o scalar_tests

table_tests: $(SRC_TABLE) ../fe.h ../point.h ../tables.h
	$(CC) $(CFLAGS) -I.. $(SRC_TABLE) -o table_tests

run: all
	./fe_tests
	./point_tests
	./scalar_tests
	./table_tests

clean:
	rm -f $(TARGETS)
