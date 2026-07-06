# Polycrystalline graphene PFC generator

This fork provides a pure C workflow from phase-field calculation to final
graphene XYZ coordinates. The Java coordinator, Python, and OVITO conversion
steps are no longer needed.

Two parallel executables are provided:

| Program | Parallel model | Target |
| --- | --- | --- |
| `polygr` / `polygr.exe` | OpenMP + FFTW threads | Windows or one Linux host |
| `polygr_mpi` | MPI + OpenMP + FFTW-MPI threads | Linux clusters and SLURM |

Both programs use the same C coordinate module (`src/xyz.c`) for periodic
local-minimum detection, Delaunay triangulation, and direct XYZ output.
The original `src/pfc.c` and Java JAR files remain as reference material.

## Build

### Windows single-host version

Install the MSYS2 UCRT64 packages:

```powershell
C:\msys64\usr\bin\bash.exe -lc "pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-fftw mingw-w64-ucrt-x86_64-make"
```

Build from PowerShell:

```powershell
.\build-windows.ps1
```

The script creates `polygr.exe` and copies the required runtime DLLs beside it.

### Linux single-host version

Install GCC, FFTW with its OpenMP library, and GNU Make, then run:

```bash
make
```

### Linux hybrid MPI version

Install an MPI implementation, FFTW built with MPI and thread support, GCC,
and GNU Make. Then run:

```bash
make mpi
```

The MPI link libraries can be overridden when a cluster uses custom paths:

```bash
make mpi CPPFLAGS="-I/opt/fftw/include" \
    LDFLAGS="-L/opt/fftw/lib" \
    MPI_LDLIBS="-lfftw3_mpi -lfftw3_omp -lfftw3 -lm"
```

## Run

Each positional run name refers to `<name>.in`. Multiple names execute
sequentially in one process or MPI job, so the original two PFC stages can be
run with one command.

### Quick Windows test

```powershell
cd example
..\polygr.exe quick1 quick2 --xyz quick.xyz --threads 8
```

### Quick MPI test

```bash
cd example
mpirun -np 2 ../polygr_mpi quick1 quick2 \
    --xyz quick-mpi.xyz \
    --threads 4
```

This uses two MPI processes and four OpenMP/FFTW threads per process.

### Full input templates

Before using `step1.in` and `step2.in`, replace `WW`, `HH`, `NN`, `RR`, and
`xabc` in `step1.in`, and replace `WW` and `HH` in `step2.in`.

Single-host command:

```bash
../polygr step1 step2 --xyz graphene.xyz --threads 8
```

Hybrid MPI command:

```bash
mpirun -np 4 ../polygr_mpi step1 step2 \
    --xyz graphene.xyz \
    --threads 12
```

The provided `slurm/polygr.slm` is a four-node example using one MPI rank and
48 OpenMP threads per node.

## Options

```text
--xyz FILE                 XYZ output path (default: <last-run>.xyz)
--no-xyz                   Skip coordinate extraction
--threads N                OpenMP and FFTW thread count per process
--pfc-lattice VALUE        PFC lattice constant (default: 7.3)
--angstrom-lattice VALUE   physical lattice constant in angstrom (default: 2.46)
```

## MPI data flow

FFTW-MPI distributes the density field as row slabs. OpenMP handles local
point-wise loops inside each MPI process. After the final relaxation:

1. each process packs its local real-space rows;
2. `MPI_Gatherv` restores the global row order on rank 0;
3. rank 0 runs the shared pure C coordinate extractor;
4. the final XYZ file is written directly.

Rank 0 must have enough memory for one complete final density field during
coordinate extraction.

## OVITO / Extended XYZ format

```text
atom_count
Lattice="..." Origin="0 0 -10" Properties="id:I:1:species:S:1:pos:R:3:radius:R:1:color:R:3" pbc="T T F"
atom_id C x y z radius red green blue
```

The output follows the Extended XYZ convention understood directly by OVITO
Basic 3.15.5. Atom IDs start at 1, the species is carbon (`C`), `z` is 0, and
coordinates are in angstrom. The file stores the actual calculation cell,
periodic boundaries in x/y, and a non-periodic 20-angstrom z dimension.
Per-particle radius and color properties make the structure immediately
visible when opened in OVITO without manual display configuration.

## Verification

The Windows quick test produces 52 carbon atoms. The C coordinate extractor
and the original `coordinator.jar` produce the same atom count, with a maximum
nearest-coordinate difference of approximately `6.3e-9` angstrom.
