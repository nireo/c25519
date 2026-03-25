PYTHON ?= python3
SINGLE_HEADER ?= c25519_single.h
SINGLE_HEADER_SCRIPT ?= scripts/make_single_header.py

.PHONY: all test tests benchmark single-header single_header clean

all: test

test:
	$(MAKE) -C tests run

tests: test

benchmark:
	$(MAKE) -C tests benchmark

single-header:
	$(PYTHON) $(SINGLE_HEADER_SCRIPT) "$(SINGLE_HEADER)"

single_header: single-header

clean:
	$(MAKE) -C tests clean
	rm -f "$(SINGLE_HEADER)"
