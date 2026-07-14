$ErrorActionPreference = "Stop"

$msys = "C:\msys64\ucrt64"
$gcc = Join-Path $msys "bin\gcc.exe"
if (-not (Test-Path $gcc)) {
    throw "MSYS2 UCRT64 GCC was not found at $gcc"
}

$env:PATH = "$(Join-Path $msys 'bin');$env:PATH"
& $gcc -O3 -std=c11 -Wall -Wextra -fopenmp `
    -I"$msys\include" "src\polygr.c" "src\xyz.c" "src\field_image.c" `
    -L"$msys\lib" -lfftw3_omp -lfftw3 -lm `
    -o "polygr.exe"
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

& $gcc -O3 -std=c11 -Wall -Wextra `
    -I"$msys\include" "src\figure.c" "src\field_image.c" `
    -L"$msys\lib" -lm -lgdi32 `
    -o "polygr_figure.exe"
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

& $gcc -O3 -std=c11 -Wall -Wextra `
    -I"$msys\include" "src\figure_samples.c" `
    -L"$msys\lib" -lm `
    -o "polygr_figure_samples.exe"
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

& $gcc -O3 -std=c11 -Wall -Wextra `
    -I"$msys\include" "src\gbmap.c" "src\field_image.c" `
    -L"$msys\lib" -lm `
    -o "polygr_gbmap.exe"
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

& $gcc -O3 -std=c11 -Wall -Wextra `
    -I"$msys\include" "src\orient_gbmap.c" "src\field_image.c" `
    -L"$msys\lib" -lm `
    -o "polygr_orient_gbmap.exe"
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

$runtimeDlls = @(
    "libfftw3-3.dll",
    "libfftw3_omp-3.dll",
    "libgomp-1.dll",
    "libgcc_s_seh-1.dll",
    "libwinpthread-1.dll"
)
foreach ($dll in $runtimeDlls) {
    Copy-Item (Join-Path $msys "bin\$dll") -Destination $dll -Force
}

Write-Host "Built $PWD\polygr.exe, $PWD\polygr_figure.exe, $PWD\polygr_figure_samples.exe, $PWD\polygr_gbmap.exe, and $PWD\polygr_orient_gbmap.exe"
