CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic

TARGETS = fe_tests point_tests
SRC_FE = test_fe.c
SRC_POINT = test_point.c

all: $(TARGETS)

fe_tests: $(SRC_FE) ../fe.h
	$(CC) $(CFLAGS) -I.. $(SRC_FE) -o fe_tests

point_tests: $(SRC_POINT) ../fe.h ../point.h
	$(CC) $(CFLAGS) -I.. $(SRC_POINT) -o point_tests

run: all
	./fe_tests
	./point_tests

clean:
	rm -f $(TARGETS)
