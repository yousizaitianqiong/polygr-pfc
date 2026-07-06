CC ?= gcc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -fopenmp
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?= -lfftw3_omp -lfftw3 -lm

BIN := polygr
SRC := src/polygr.c

.PHONY: all clean

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) $(BIN) $(BIN).exe
