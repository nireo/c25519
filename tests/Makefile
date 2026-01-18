CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic

TARGET = fe_tests
SRC = test_fe.c

all: $(TARGET)

$(TARGET): $(SRC) ../fe.h
	$(CC) $(CFLAGS) -I.. $(SRC) -o $(TARGET)

run: all
	./$(TARGET)

clean:
	rm -f $(TARGET)
