# Polycrystalline graphene PFC generator

This repository now includes a pure C end-to-end implementation:

- phase-field crystal relaxation in C;
- shared-memory parallelism through OpenMP and FFTW threads;
- periodic local-minimum detection and Delaunay triangulation in C;
- direct XYZ output after the last calculation stage;
- no Java, MPI, OVITO, or Python runtime dependency.

The original MPI source (`src/pfc.c`) and Java JAR files are retained only as
references. The new program is `src/polygr.c`.

## Build

### Windows

Install MSYS2 UCRT64 packages:

```powershell
C:\msys64\usr\bin\bash.exe -lc "pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-fftw mingw-w64-ucrt-x86_64-make"
```

Then build from PowerShell:

```powershell
.\build-windows.ps1
```

The script creates `polygr.exe` and copies its required runtime DLLs beside it.

### Linux

Install GCC, FFTW (including its OpenMP library), and GNU Make, then run:

```bash
make
```

## Run

The input file format is compatible with the original PFC program. One or more
run names may be supplied. Each name refers to `<name>.in`.

Run the included small two-stage test:

```powershell
cd example
..\polygr.exe quick1 quick2 --xyz quick.xyz --threads 4
```

For the original two-stage workflow, first replace `WW`, `HH`, `NN`, `RR`, and
`xabc` in `example/step1.in`, and replace `WW` and `HH` in `example/step2.in`.
Then run:

```powershell
..\polygr.exe step1 step2 --xyz graphene.xyz --threads 8
```

Options:

```text
--xyz FILE                 XYZ output path (default: <last-run>.xyz)
--no-xyz                   Skip coordinate extraction
--threads N                OpenMP and FFTW thread count
--pfc-lattice VALUE        PFC lattice constant (default: 7.3)
--angstrom-lattice VALUE   physical lattice constant in angstrom (default: 2.46)
```

## XYZ format

The output is written directly after the final relaxation:

```text
atom_count

atom_id atom_type x y z
```

Atom IDs start at 1, atom type is 1, and `z` is 0. Coordinates are in angstrom.

## Verification

For the included quick test, the C coordinate extractor and the original
`coordinator.jar` both produce 52 atoms. The maximum nearest-coordinate
difference is approximately `6.3e-9` angstrom.
