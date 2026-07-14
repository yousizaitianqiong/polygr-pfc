# Polycrystalline graphene PFC generator

This fork provides a pure C workflow from phase-field calculation to graphene
XYZ coordinates, grayscale PNG diagnostics, and per-frame density-field data
that can be post-processed. Scientific results should be archived from the
numeric `.dat` and Extended XYZ outputs; PNG figures are only presentation
products.

Two parallel executables are provided:

| Program | Parallel model | Target |
| --- | --- | --- |
| `polygr` / `polygr.exe` | OpenMP + FFTW threads | Windows or one Linux host |
| `polygr_mpi` | MPI + OpenMP + FFTW-MPI threads | Linux clusters and SLURM |

Both programs use the same C coordinate module (`src/xyz.c`) for periodic
local-minimum detection, Delaunay triangulation, and direct XYZ output.
They also share a small built-in PNG writer (`src/field_image.c`) for density
field images.
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

### CUDA acceleration notes

CUDA acceleration is not part of the default build. The optional single-GPU
`polygr_cuda` executable uses cuFFT and CUDA kernels for the PFC solver loop.
See `docs/cuda-optimization.md` for the target design, hot paths, validation
checks, and current backend limitations.

On Windows, `nvidia-smi` only confirms that the NVIDIA driver can run CUDA
programs. Building CUDA code also requires the CUDA Toolkit, including
`nvcc.exe` and cuFFT development libraries. After installing the Toolkit, build
the CUDA probe/backend entry with:

```powershell
.\build-cuda-windows.ps1
```

The current CUDA backend is intended for solver validation and benchmarking.
Run it with `--no-xyz --no-png`; final XYZ/PNG postprocessing remains in the
CPU executable.

## Run

Each positional run name refers to `<name>.in`. Multiple names execute
sequentially in one process or MPI job, so the original two PFC stages can be
run with one command.

### Quick Windows test

```powershell
cd example
..\polygr.exe quick1 quick2 --xyz quick.xyz --png quick.png --threads 8
```

### Larger OVITO validation test

The `medium1.in` and `medium2.in` inputs use a `144 x 144` grid and 18 initial
grains. They produce a denser, more useful structure for checking the final
XYZ file in OVITO Basic 3.15.5 while still running comfortably on a Windows
workstation.

```powershell
cd example
..\polygr.exe medium1 medium2 --xyz medium-ovito.xyz --png medium-ovito.png --threads 8
```

To generate an OVITO-readable coordinate sequence for every written PFC
snapshot, add `--xyz-frames`:

```powershell
cd example
..\polygr.exe medium1 medium2 --xyz medium-ovito.xyz --png medium-ovito.png `
    --xyz-frames medium-frame --threads 8
```

This keeps the original density fields, for example `medium1-t-0.dat` and
`medium2-t-200.dat`, and also writes Extended XYZ frames such as
`medium-frame-medium1-t-0.xyz` and `medium-frame-medium2-t-200.xyz`.
Open the XYZ sequence in OVITO for visual inspection and defect analysis. The
XYZ files include `coordination`, `ring0`, `ring1`, `ring2`, and `gb`
properties, so atoms next to non-six-member rings or abnormal coordination are
marked as grain-boundary candidates. Use the PNG only as a quick density-field
diagnostic.

### old_pfc-compatible two-stage workflow

This fork follows the `petenez/old_pfc` workflow: first grow randomly placed
honeycomb crystallites close to melting with conserved dynamics, then reload
that state and relax it further with nonconserved dynamics and box-size
optimization. The paper connection is the one-mode PFC generation path used for
large polycrystalline graphene samples, plus conversion of relaxed density-field
minima into carbon coordinates and nonhexagonal rings.

The original Java `plotter.jar` only linearly maps density-field minima and
maxima to grayscale images. The original Java `coordinator.jar` detects local
density minima, finds nearest-neighbor triplets of those minima, places carbon
atoms at the triplet centers, then writes carbon coordinates and nonhexagonal
ring membership. This fork keeps equivalent C implementations for grayscale PNG
output and XYZ coordinate extraction.

### Mentor-reference CVD composite figure test

The `paper_cvd1.in` and `paper_cvd2.in` inputs are tuned for the composite
figure workflow. They keep the `288 x 192` grid but use 128 initial nuclei, a
smaller initial radius, and a 50-step write interval so panel `(a)` can be built
from a smoother PFC growth sequence instead of widely spaced snapshots. The
target visual style is the mentor-provided CVD graphene reference figure: a
top-row growth sequence ending in a relaxed nc-graphene orientation map, plus
panels for island-spacing scaling, a stress/orientation sheet schematic, and a
grain-boundary network map.

```powershell
cd example
..\polygr.exe paper_cvd1 paper_cvd2 --xyz paper-cvd.xyz --png paper-cvd.png --png-scale 4 --threads 8
```

The composite figure is generated by the pure C postprocessor
`polygr_figure`. It reads four growth snapshots plus one relaxed snapshot,
then writes a single PNG with panels `(a)`, `(b)`, `(c)`, and `(d)`: grayscale
density growth images, a final orientation map, a CSV-driven island-spacing panel,
a 3D-style orientation sheet, and a grain-boundary map. This is a visual
summary for the generated PFC sample and should be treated as a CVD composite
figure, not as a literal reproduction of PRB 94, 035414 Figure 1. The PRB
Figure 1 compares equilibrium density fields from four PFC variants
(`PFC1`, `APFC`, `PFC3`, and `XPFC`), while this fork follows the
`old_pfc`-compatible one-mode growth-and-relaxation path.

```powershell
..\polygr_figure.exe --width 288 --height 192 --output paper-cvd-figure.png `
    --style pfc `
    --panel-b-csv ..\batch-panel-b\panel-b-stats.csv `
    paper_cvd1-t-0.dat paper_cvd1-t-50.dat paper_cvd1-t-150.dat `
    paper_cvd1-t-300.dat paper_cvd2-t-300.dat
```

Use the following accepted command sequence when regenerating the main composite
figure for visual checks. The first four panel `(a)` frames are fixed to
`panel_a_seedgrow1` at `t=0,5,20,60`, and the fifth relaxed nc-graphene frame is
fixed to `panel_a_seedgrow2-t-300.dat`. Prefer this input order when checking
the composite figure so the intermediate `panel_a_seedgrow2` frames are not
mistaken for the accepted first four growth frames.

```powershell
cd D:\09yuan\polygr-pfc-fork\example
..\polygr.exe panel_a_seedgrow1 panel_a_seedgrow2 --threads 8
..\polygr_figure.exe --style pfc --width 432 --height 288 `
    --output panel-a-growth-pfc-dat-check.png `
    panel_a_seedgrow1-t-0.dat `
    panel_a_seedgrow1-t-5.dat `
    panel_a_seedgrow1-t-20.dat `
    panel_a_seedgrow1-t-60.dat `
    panel_a_seedgrow2-t-300.dat
```

Those snapshots correspond to random nucleation, small island growth, island
coalescence, a near-continuous polycrystalline film, and the final relaxed
nc-graphene orientation map.

Panel data sources:

1. panel `(a)` defaults to `--style pfc`, which directly renders the real PFC
   `.dat` snapshots supplied on the command line. The first four frames are
   true density-field growth frames displayed with a density-threshold and
   local-peak postprocessor so the images resemble OVITO-style atom/island
   views rather than raw continuous grayscale fields. The fifth frame is a
   relaxed local orientation map computed from the final `.dat`; panel `(a)`
   uses a smaller orientation window than panels `(c)` and `(d)` so the top-row
   relaxed map retains more local grain detail. Its colors encode the local
   lattice orientation angle `theta(x,y)` in the range `-pi/12` to `pi/12`;
   different colors mean different in-plane crystalline orientations, not
   different chemical species or arbitrary labels. Use `--style schematic` only
   for a reference-style fixed-Voronoi explanatory graphic.
2. panel `(b)` reads batch statistics from CSV via `--panel-b-csv`. The minimum
   required columns are `seed,nucleus_count,a_s_nm,grain_size_nm`; extra
   provenance columns are ignored by the figure renderer.
3. panel `(c)` is an orientation field with a deterministic height proxy. It is
   not a true 3D buckling calculation or a true stress calculation.
4. panel `(d)` is extracted from local orientation discontinuities in the final
   PFC density field; the red lines are not painted on top of a PNG.

To generate panel `(b)` CSV data on Windows, run a batch of small PFC simulations
over several initial nucleus counts and random seeds:

```powershell
cd D:\09yuan\polygr-pfc-fork
.\scripts\run-panel-b-batch.ps1 `
    -NucleusCounts 36,72,128,180 `
    -Seeds 101,202,303 `
    -OutputDir batch-panel-b `
    -Csv batch-panel-b\panel-b-stats.csv `
    -Threads 4
```

The batch runner writes one row per run. `a_s_nm` is the initial island spacing
estimated from the final cell area and `nucleus_count`. `grain_size_nm` is a
grain-boundary spacing estimate from the final Extended XYZ: the script counts
`gb=1` atoms and converts that boundary-atom density into an areal spacing. The
CSV preserves the final `.xyz` and `.dat` paths so each point can be traced back
to its PFC output.

For a larger, high-nucleation real PFC example that gives a more reference-like
top row and an OVITO-checkable grain-boundary network, use `gb_reference1.in`
and `gb_reference2.in`:

```powershell
cd D:\09yuan\polygr-pfc-fork\example
..\polygr.exe gb_reference1 gb_reference2 `
    --xyz gb-reference-final.xyz `
    --png gb-reference-final.png `
    --xyz-frames gb-reference-frame `
    --threads 4

cd D:\09yuan\polygr-pfc-fork
.\polygr_figure.exe --width 432 --height 288 --style pfc `
    --panel-b-csv batch-panel-b-wide\panel-b-stats-wide.csv `
    --output example\composite-reference-inputs-pfc.png `
    example\gb_reference1-t-0.dat `
    example\gb_reference1-t-60.dat `
    example\gb_reference1-t-180.dat `
    example\gb_reference1-t-360.dat `
    example\gb_reference2-t-300.dat

.\polygr_gbmap.exe example\gb-reference-final.xyz `
    --width 432 --height 288 `
    --output example\gb-reference-ovito-gb-network.png
```

For fast visual testing without running the PFC solver, `polygr_figure_samples`
generates a synthetic CVD growth series in the same `.dat` format. It uses
random nucleation, radial growth, and Voronoi grain assignment so panel `(a)`
shows sparse nuclei, isolated islands, island coalescence, near-continuous film
formation, and a final relaxed nc-graphene orientation map.

```powershell
..\polygr_figure_samples.exe --width 288 --height 192 --grains 36 --prefix discrete
..\polygr_figure.exe --width 288 --height 192 --output discrete-figure.png `
    discrete-t-0.dat discrete-t-1.dat discrete-t-2.dat `
    discrete-t-3.dat discrete-t-4.dat
```

For a denser reference-style visual stress test, increase the synthetic grain
count. This produces a more complex orientation map, denser boundary network,
and a stronger 3D surface in panel `(c)`.

```powershell
..\polygr_figure_samples.exe --width 288 --height 192 --grains 144 --prefix reference-rich
..\polygr_figure.exe --width 288 --height 192 --output reference-rich-figure.png `
    reference-rich-t-0.dat reference-rich-t-1.dat reference-rich-t-2.dat `
    reference-rich-t-3.dat reference-rich-t-4.dat
```

Use `--plot-variant scatter` when the comparison panel should represent a
more dispersed built-in fallback trend instead of the default plateau-to-growth
fallback trend. For production figures, prefer `--panel-b-csv`:

```powershell
..\polygr_figure.exe --width 288 --height 192 --plot-variant scatter `
    --output reference-rich-scatter-figure.png `
    reference-rich-t-0.dat reference-rich-t-1.dat reference-rich-t-2.dat `
    reference-rich-t-3.dat reference-rich-t-4.dat
```

### Quick MPI test

```bash
cd example
mpirun -np 2 ../polygr_mpi quick1 quick2 \
    --xyz quick-mpi.xyz \
    --png quick-mpi.png \
    --threads 4
```

This uses two MPI processes and four OpenMP/FFTW threads per process.

### Full input templates

Before using `step1.in` and `step2.in`, replace `WW`, `HH`, `NN`, `RR`, and
`xabc` in `step1.in`, and replace `WW` and `HH` in `step2.in`.

Single-host command:

```bash
../polygr step1 step2 --xyz graphene.xyz --png graphene.png --threads 8
```

Hybrid MPI command:

```bash
mpirun -np 4 ../polygr_mpi step1 step2 \
    --xyz graphene.xyz \
    --png graphene.png \
    --threads 12
```

The provided `slurm/polygr.slm` is a four-node example using one MPI rank and
48 OpenMP threads per node.

## Options

```text
--xyz FILE                 Final XYZ output path (default: <last-run>.xyz)
--no-xyz                   Skip coordinate extraction
--xyz-frames PREFIX        Write an Extended XYZ file for every saved .dat frame
--png FILE                 grayscale density PNG path (default: <last-run>.png)
--no-png                   Skip density PNG output
--png-scale N              bilinear PNG scale factor (default: 4)
--threads N                OpenMP and FFTW thread count per process
--pfc-lattice VALUE        PFC lattice constant (default: 7.3)
--angstrom-lattice VALUE   physical lattice constant in angstrom (default: 2.46)
```

`polygr_figure` additionally accepts:

```text
--panel-b-csv FILE         CSV for panel (b), requiring a_s_nm and grain_size_nm
--style pfc|schematic      pfc is real .dat post-processing; schematic is fixed-Voronoi explanatory rendering
--plot-variant NAME        built-in fallback panel (b): plateau or scatter
```

## MPI data flow

FFTW-MPI distributes the density field as row slabs. OpenMP handles local
point-wise loops inside each MPI process. After the final relaxation:

1. each process packs its local real-space rows;
2. `MPI_Gatherv` restores the global row order on rank 0;
3. rank 0 runs the shared pure C coordinate extractor and PNG writer;
4. the final XYZ and grayscale PNG files are written directly.

Rank 0 must have enough memory for one complete final density field during
coordinate extraction and PNG output.

## Grayscale density PNG

The final density field is written as an 8-bit grayscale PNG. Pixel intensity is
linearly normalized from the final field minimum to maximum, matching the
default scalar-field path used by the original Java `plotter.jar`. The image is
flipped vertically so it has the same orientation as the Java plotter output.
The default `--png-scale 4` uses bilinear interpolation for a smoother figure
export; use `--png-scale 1` for a one-pixel-per-grid-point diagnostic image.

## OVITO / Extended XYZ workflow

The paper workflow is data-first: relax PFC density fields, extract atomic
coordinates from the density extrema, then analyze grain-boundary defects such
as 5|7 and 5|8|7 rings. For supervisor-facing results, deliver:

1. the input files and random seeds;
2. all saved `*-t-*.dat` density frames;
3. matching per-frame Extended XYZ files from `--xyz-frames`;
4. OVITO screenshots or exported images generated from those XYZ files;
5. a short note describing which frame, coloring, and defect criterion were
   used.

In OVITO, open the final or per-frame `.xyz` file. The file already declares
periodic boundaries in `x` and `y` and colors `gb=1` atoms red. For a cleaner
view, add an `Expression selection` modifier with `gb == 1`, or color by the
`coordination`, `ring0`, `ring1`, and `ring2` particle properties. For
graphene, grain boundaries should be identified from non-hexagonal ring
statistics or coordination/topological analysis, not manually painted onto a
PNG. In this exporter, each carbon atom stores the sizes of its three adjacent
PFC rings; any adjacent ring size other than 6 or any coordination other than 3
sets `gb=1`.

To make a clean panel `(d)` style grain-boundary network from the same OVITO
XYZ data, use `polygr_gbmap`. It projects only `gb=1` atoms to a white
background, connects nearby boundary atoms with thin red lines, and optionally
draws the `N_x`/`N_y` grid:

```powershell
..\polygr_gbmap.exe medium-ovito-final.xyz `
    --output medium-ovito-gb-map.png --width 512 --height 512 --grid 4 4
```

This PNG is the figure-panel export. The corresponding `.xyz` remains the
post-processing data source that can be reopened in OVITO.

For a closer match to the reference panel `(d)`, use the final PFC density
field directly and extract the grain-boundary centerlines from local
orientation jumps:

```powershell
..\polygr_orient_gbmap.exe gb_large2-t-300.dat `
    --width 432 --height 288 `
    --output gb-large-orient-gb-map-thin-t010.png `
    --image-width 768 --image-height 512 --grid 4 4 `
    --threshold 0.10 --radius 9
```

`polygr_gbmap` shows the atomistic defect network from `gb=1` atoms in the
OVITO-readable XYZ file. `polygr_orient_gbmap` shows the mesoscopic continuous
grain-boundary centerline network from `theta(x,y)` discontinuities in the PFC
density field, which is the better choice for panel `(d)`.

## Extended XYZ format

```text
atom_count
Lattice="..." Origin="0 0 -10" Properties="id:I:1:species:S:1:pos:R:3:coordination:I:1:ring0:I:1:ring1:I:1:ring2:I:1:gb:I:1:Radius:R:1:Color:R:3" pbc="T T F"
atom_id C x y z coordination ring0 ring1 ring2 gb Radius red green blue
```

The output follows the Extended XYZ convention understood directly by OVITO
Basic 3.15.5. Atom IDs start at 1, the species is carbon (`C`), `z` is 0, and
coordinates are in angstrom. The file stores the actual calculation cell,
periodic boundaries in x/y, and a non-periodic 20-angstrom z dimension.
Per-particle radius and color properties make the structure immediately
visible when opened in OVITO without manual display configuration. Normal
three-coordinated atoms adjacent only to six-member rings are gray; candidate
grain-boundary atoms are red.

## Verification

The Windows quick test produces 52 carbon atoms. The C coordinate extractor
and the original `coordinator.jar` produce the same atom count, with a maximum
nearest-coordinate difference of approximately `6.3e-9` angstrom.

On Windows, run the full smoke test and PNG figure generation check with:

```powershell
.\tests\run-tests.ps1
```

The script runs the quick two-stage solver, validates the Extended XYZ header,
validates generated PNG metadata, creates synthetic CVD growth snapshots, and
writes default and scatter composite PNG figures under `test-output`.
