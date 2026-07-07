# CUDA optimization plan

This project is currently optimized for a single host with OpenMP plus FFTW
threads, and for Linux clusters with MPI plus FFTW-MPI. CUDA can help most when
the full PFC time-step loop stays resident on the GPU. Moving only one pointwise
loop to CUDA is not recommended because every step would pay host/device copy
costs while the FFTs still run on the CPU.

## Best target

The separate single-GPU executable `polygr_cuda` mirrors the solver portion of
`polygr` and uses:

- `cufftPlan2d` for the real-to-complex and complex-to-real transforms.
- Device arrays for `p`, `q`, `linear`, and `nonlinear`.
- CUDA kernels for pointwise real-space and frequency-space operations.
- Host copies for input initialization, `.dat` snapshots, and progress output.

The existing `polygr` and `polygr_mpi` paths remain unchanged. CUDA is a
single-node accelerator path; it does not replace the MPI path unless a future
multi-GPU design is added.

Current limitation: `polygr_cuda` is for solver validation and benchmarking.
Pass `--no-xyz --no-png`; final XYZ and PNG postprocessing still use the CPU
executable.

## Hot paths to move together

The CPU implementation is in `src/polygr.c`.

1. `update_operators`
   - Computes one real coefficient pair per half-spectrum point.
   - Good CUDA kernel candidate.

2. `solver_step`
   - Computes `gamma * q^2 + delta * q^3` in real space.
   - Runs `q_forward`.
   - Applies the semi-implicit update to `p_frequency` and copies it to
     `q_frequency`.
   - Runs `q_inverse`.
   - This is the main CUDA target and should be kept entirely on device.

3. `calculate_properties`
   - Needed for progress output and box optimization.
   - Can start with CPU fallback by copying the field back at print intervals,
     but this limits performance when `print_interval` or `optimize_interval`
     is small.
   - A later CUDA version should use device reductions for free energy and
     density.

4. `optimize_box`
   - Calls `calculate_properties` five times, so it can dominate short runs.
   - Prefer keeping `calculate_properties` on GPU before optimizing this path.

## cuFFT normalization

FFTW and cuFFT both leave transforms unnormalized. The current CPU code scales
after forward transforms:

```c
scale_complex(..., 1.0 / (w * h));
```

The CUDA path must preserve the same convention. Either scale after R2C
transforms or fold the same factor into the frequency-space kernels. Keep this
choice explicit, because it is the easiest place to introduce result drift.

## Data layout

The CPU path allocates real arrays with FFTW's padded R2C layout:

```text
real_count    = h * (2 * (w / 2 + 1))
complex_count = h * (w / 2 + 1)
stride        = 2 * (w / 2 + 1)
```

Use the same logical layout with cuFFT:

- Real field: `double *`, length `h * stride`.
- Complex field: reinterpret or store as `cufftDoubleComplex *`, length
  `h * (w / 2 + 1)`.

Keeping this layout avoids changing state-file, PNG, and XYZ code.

## Build shape

Recommended optional target:

```make
NVCC ?= nvcc
CUDA_CFLAGS ?= -O3
CUDA_LDLIBS ?= -lcufft -lcudart -lm

cuda: polygr_cuda
```

Do not make CUDA part of the default `make all`; many Windows and cluster
systems can build the current CPU/MPI tools without an NVIDIA toolchain.

On Windows, use:

```powershell
.\build-cuda-windows.ps1
```

If this reports that `nvcc.exe` was not found, the machine may still have a
CUDA-capable GPU and driver. Install the CUDA Toolkit or add its `bin`
directory to `PATH` before building CUDA code.

## Expected wins

CUDA should help when the grid is large enough for FFT cost to dominate, such as
the `paper*.in` and production `step*.in` workloads. The `quick*.in` workload is
too small to be a meaningful GPU benchmark.

Use these checks when validating a CUDA backend:

1. Run the CPU quick test and the CUDA quick test from a clean `example`
   directory.
2. Compare final `.out` free-energy and density values with a small tolerance.
3. Compare final XYZ atom count.
4. Benchmark larger runs with `--no-xyz --no-png` first to isolate solver time.
5. Repeat with normal output enabled to measure transfer and postprocessing
   overhead.

## Suggested implementation order

1. Create `src/polygr_cuda.cu` by copying only the single-host solver path from
   `src/polygr.c`.
2. Replace FFTW allocation and plans with CUDA allocation and cuFFT plans.
3. Implement CUDA kernels for `update_operators`, `scale_complex`,
   real-space nonlinear update, frequency-space semi-implicit update, and final
   inverse-transform scaling if needed.
4. Keep file I/O, argument parsing, XYZ, and PNG output on the host.
5. Add a `make cuda` target and document required CUDA Toolkit version.
6. Add a quick regression command that compares CPU and CUDA outputs.

This staged approach avoids a half-CPU, half-GPU design that is usually slower
than the current OpenMP/FFTW solver.
