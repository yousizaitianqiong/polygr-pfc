CC ?= gcc
MPICC ?= mpicc
NVCC ?= nvcc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -fopenmp
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?= -lfftw3_omp -lfftw3 -lm
MPI_LDLIBS ?= -lfftw3_mpi -lfftw3_omp -lfftw3 -lm
CUDA_CFLAGS ?= -O3
CUDA_LDLIBS ?= -lcufft -lcudart -lm
FIGURE_LDLIBS ?= -lm
ifeq ($(OS),Windows_NT)
FIGURE_LDLIBS += -lgdi32
endif

BIN := polygr
MPI_BIN := polygr_mpi
CUDA_BIN := polygr_cuda
FIGURE_BIN := polygr_figure
SAMPLE_BIN := polygr_figure_samples
GBMAP_BIN := polygr_gbmap
ORIENT_GBMAP_BIN := polygr_orient_gbmap
COMMON_SRC := src/xyz.c src/field_image.c
SRC := src/polygr.c $(COMMON_SRC)
MPI_SRC := src/polygr_mpi.c $(COMMON_SRC)
FIGURE_SRC := src/figure.c src/field_image.c
SAMPLE_SRC := src/figure_samples.c
GBMAP_SRC := src/gbmap.c src/field_image.c
ORIENT_GBMAP_SRC := src/orient_gbmap.c src/field_image.c

.PHONY: all mpi cuda clean

all: $(BIN) $(FIGURE_BIN) $(SAMPLE_BIN) $(GBMAP_BIN) $(ORIENT_GBMAP_BIN)

$(BIN): $(SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SRC) $(LDFLAGS) $(LDLIBS) -o $@

mpi: $(MPI_BIN)

$(MPI_BIN): $(MPI_SRC)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(MPI_SRC) $(LDFLAGS) $(MPI_LDLIBS) -o $@

cuda: $(CUDA_BIN)

$(CUDA_BIN): src/polygr_cuda.cu
	$(NVCC) $(CUDA_CFLAGS) $(CPPFLAGS) src/polygr_cuda.cu $(LDFLAGS) $(CUDA_LDLIBS) -o $@

$(FIGURE_BIN): $(FIGURE_SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(FIGURE_SRC) $(LDFLAGS) $(FIGURE_LDLIBS) -o $@

$(SAMPLE_BIN): $(SAMPLE_SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SAMPLE_SRC) $(LDFLAGS) -lm -o $@

$(GBMAP_BIN): $(GBMAP_SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(GBMAP_SRC) $(LDFLAGS) -lm -o $@

$(ORIENT_GBMAP_BIN): $(ORIENT_GBMAP_SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(ORIENT_GBMAP_SRC) $(LDFLAGS) -lm -o $@

clean:
	$(RM) $(BIN) $(BIN).exe $(MPI_BIN) $(CUDA_BIN) $(CUDA_BIN).exe $(FIGURE_BIN) $(FIGURE_BIN).exe $(SAMPLE_BIN) $(SAMPLE_BIN).exe $(GBMAP_BIN) $(GBMAP_BIN).exe $(ORIENT_GBMAP_BIN) $(ORIENT_GBMAP_BIN).exe
