$ErrorActionPreference = "Stop"

$cl = Get-Command cl.exe -ErrorAction SilentlyContinue
if (-not $cl -and -not $env:POLYGR_CUDA_VSDEV) {
    $vsDevCmd = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat"
    if (Test-Path $vsDevCmd) {
        $env:POLYGR_CUDA_VSDEV = "1"
        & cmd /c "call `"$vsDevCmd`" -arch=x64 && powershell -NoProfile -ExecutionPolicy Bypass -File `"$PSCommandPath`""
        exit $LASTEXITCODE
    }
}

$cl = Get-Command cl.exe -ErrorAction SilentlyContinue
if (-not $cl) {
    throw "MSVC cl.exe was not found. Install Visual Studio Build Tools with the C++ toolchain."
}

$nvcc = Get-Command nvcc.exe -ErrorAction SilentlyContinue
if (-not $nvcc) {
    $candidateRoots = @(
        "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA",
        "D:\Program Files\NVIDIA GPU Computing Toolkit\CUDA"
    )
    foreach ($root in $candidateRoots) {
        if (Test-Path $root) {
            $candidate = Get-ChildItem $root -Recurse -Filter nvcc.exe -ErrorAction SilentlyContinue |
                Select-Object -First 1
            if ($candidate) {
                $nvcc = Get-Command $candidate.FullName
                break
            }
        }
    }
}

if (-not $nvcc) {
    throw "CUDA Toolkit nvcc.exe was not found. The NVIDIA driver can run CUDA programs, but building polygr_cuda.exe requires the CUDA Toolkit."
}

& $nvcc.Source -O3 `
    -Xcompiler "/utf-8" `
    "src\polygr_cuda.cu" `
    -lcufft -lcudart `
    -o "polygr_cuda.exe"
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

Write-Host "Built $PWD\polygr_cuda.exe"
