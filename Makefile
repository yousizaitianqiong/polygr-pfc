CC ?= gcc
MPICC ?= mpicc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -fopenmp
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?= -lfftw3_omp -lfftw3 -lm
MPI_LDLIBS ?= -lfftw3_mpi -lfftw3_omp -lfftw3 -lm

BIN := polygr
MPI_BIN := polygr_mpi
COMMON_SRC := src/xyz.c
SRC := src/polygr.c $(COMMON_SRC)
MPI_SRC := src/polygr_mpi.c $(COMMON_SRC)

.PHONY: all mpi clean

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SRC) $(LDFLAGS) $(LDLIBS) -o $@

mpi: $(MPI_BIN)

$(MPI_BIN): $(MPI_SRC)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(MPI_SRC) $(LDFLAGS) $(MPI_LDLIBS) -o $@

clean:
	$(RM) $(BIN) $(BIN).exe $(MPI_BIN)
