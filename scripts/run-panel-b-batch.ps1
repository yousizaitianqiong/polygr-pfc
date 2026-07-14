param(
    [int[]]$NucleusCounts = @(36, 72, 128, 180),
    [int[]]$Seeds = @(101, 202, 303),
    [string]$OutputDir = "batch-panel-b",
    [string]$Csv = "batch-panel-b\panel-b-stats.csv",
    [int]$Width = 288,
    [int]$Height = 192,
    [int]$GrowthSteps = 120,
    [int]$RelaxSteps = 80,
    [int]$Threads = 4,
    [double]$Dx = 0.7,
    [double]$Dy = 0.7,
    [double]$PfcLattice = 7.2552,
    [double]$AngstromLattice = 2.46
)

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Polygr = Join-Path $RepoRoot "polygr.exe"
if (-not (Test-Path -LiteralPath $Polygr -PathType Leaf)) {
    throw "polygr.exe was not found. Run build-windows.ps1 first."
}

$OutputPath = if ([System.IO.Path]::IsPathRooted($OutputDir)) {
    $OutputDir
} else {
    Join-Path $RepoRoot $OutputDir
}
$CsvPath = if ([System.IO.Path]::IsPathRooted($Csv)) {
    $Csv
} else {
    Join-Path $RepoRoot $Csv
}

New-Item -ItemType Directory -Path $OutputPath -Force | Out-Null
New-Item -ItemType Directory -Path (Split-Path -Parent $CsvPath) -Force | Out-Null

function Get-XyzStats {
    param(
        [Parameter(Mandatory = $true)][string]$Path,
        [Parameter(Mandatory = $true)][int]$NucleusCount,
        [Parameter(Mandatory = $true)][double]$FallbackAreaA2,
        [Parameter(Mandatory = $true)][double]$AngstromLattice
    )

    $reader = [System.IO.File]::OpenText($Path)
    try {
        [void]$reader.ReadLine()
        $header = $reader.ReadLine()
        $areaA2 = $FallbackAreaA2
        if ($header -match 'Lattice="([^"]+)"') {
            $parts = $Matches[1] -split '\s+' | ForEach-Object { [double]$_ }
            if ($parts.Count -ge 5) {
                $areaA2 = [math]::Abs($parts[0] * $parts[4] - $parts[1] * $parts[3])
            }
        }

        $atomCount = 0
        $gbCount = 0
        while (($line = $reader.ReadLine()) -ne $null) {
            if ([string]::IsNullOrWhiteSpace($line)) {
                continue
            }
            $fields = $line -split '\s+'
            if ($fields.Count -lt 10) {
                continue
            }
            $atomCount++
            if ([int]$fields[9] -eq 1) {
                $gbCount++
            }
        }
    }
    finally {
        $reader.Close()
    }

    $areaNm2 = $areaA2 / 100.0
    $asNm = [math]::Sqrt($areaNm2 / [double]$NucleusCount)
    if ($gbCount -gt 0) {
        $boundaryLengthNm = $gbCount * $AngstromLattice / 10.0
        $grainSizeNm = $areaNm2 / $boundaryLengthNm
    } else {
        $grainSizeNm = $asNm
    }

    [pscustomobject]@{
        AtomCount = $atomCount
        GbCount = $gbCount
        AreaNm2 = $areaNm2
        ASNm = $asNm
        GrainSizeNm = $grainSizeNm
    }
}

$rows = New-Object System.Collections.Generic.List[object]
$previousLocation = Get-Location
Set-Location $OutputPath
try {
    foreach ($nucleusCount in $NucleusCounts) {
        foreach ($seed in $Seeds) {
            $prefix = "n${nucleusCount}-s${seed}"
            $growth = "$prefix-growth"
            $relax = "$prefix-relax"
            $xyz = "$prefix-final.xyz"
            $png = "$prefix-final.png"

            @"
S $seed
O $GrowthSteps $GrowthSteps
A $Width $Height
I 1 $Dx $Dy $PfcLattice 0.25 -0.18 $nucleusCount 0.58
M -0.18 1.0 0.0 1.0 1
R $GrowthSteps $Dx $Dy 0.2 0
"@ | Set-Content -LiteralPath "$growth.in" -Encoding ASCII

            @"
S $seed
O $RelaxSteps $RelaxSteps
A $Width $Height
I 2 $growth-t-$GrowthSteps.dat -0.072 0.1
M -0.15 1.0 0.87 1.0 0
R $RelaxSteps $Dx $Dy 1.0 100
"@ | Set-Content -LiteralPath "$relax.in" -Encoding ASCII

            Write-Host "RUN nucleus_count=$nucleusCount seed=$seed"
            & $Polygr $growth $relax --xyz $xyz --png $png --png-scale 2 --threads $Threads --pfc-lattice $PfcLattice --angstrom-lattice $AngstromLattice
            if ($LASTEXITCODE -ne 0) {
                throw "polygr failed for nucleus_count=$nucleusCount seed=$seed"
            }

            $fallbackAreaA2 = $Width * $Dx * $Height * $Dy * ($AngstromLattice / $PfcLattice) * ($AngstromLattice / $PfcLattice)
            $stats = Get-XyzStats -Path (Join-Path $OutputPath $xyz) -NucleusCount $nucleusCount -FallbackAreaA2 $fallbackAreaA2 -AngstromLattice $AngstromLattice
            $rows.Add([pscustomobject]@{
                seed = $seed
                nucleus_count = $nucleusCount
                a_s_nm = [math]::Round($stats.ASNm, 6)
                grain_size_nm = [math]::Round($stats.GrainSizeNm, 6)
                atom_count = $stats.AtomCount
                gb_atom_count = $stats.GbCount
                area_nm2 = [math]::Round($stats.AreaNm2, 6)
                final_xyz = $xyz
                final_dat = "$relax-t-$RelaxSteps.dat"
            }) | Out-Null
        }
    }
}
finally {
    Set-Location $previousLocation
}

$rows | Export-Csv -LiteralPath $CsvPath -NoTypeInformation
Write-Host "Wrote $CsvPath"
