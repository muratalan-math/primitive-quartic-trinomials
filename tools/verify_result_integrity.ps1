param(
    [string]$ResultDir = ""
)

$ErrorActionPreference = "Stop"

function Write-Line {
    param([string]$Text = "")
    Add-Content -LiteralPath $script:ReportPath -Value $Text -Encoding UTF8
}

function Test-GzipFile {
    param([Parameter(Mandatory=$true)][string]$Path)
    $fs = $null
    $gz = $null
    try {
        $fs = [System.IO.File]::OpenRead($Path)
        $gz = New-Object System.IO.Compression.GZipStream(
            $fs,
            [System.IO.Compression.CompressionMode]::Decompress
        )
        $buffer = New-Object byte[] 1048576
        while (($n = $gz.Read($buffer, 0, $buffer.Length)) -gt 0) { }
        return $true
    }
    catch {
        return $false
    }
    finally {
        if ($gz -ne $null) { $gz.Dispose() }
        if ($fs -ne $null) { $fs.Dispose() }
    }
}

if ([string]::IsNullOrWhiteSpace($ResultDir)) {
    $ResultDir = $PSScriptRoot
}
$ResultDir = [System.IO.Path]::GetFullPath($ResultDir.Trim('"'))

$script:ReportPath = Join-Path $ResultDir "RUN_ALL_INTEGRITY_REPORT.txt"
$AllManifest = Join-Path $ResultDir "MANIFEST_SHA256_ALL.txt"
$ZenodoManifest = Join-Path $ResultDir "MANIFEST_SHA256_ZENODO.txt"
$ZenodoList = Join-Path $ResultDir "ZENODO_FILE_LIST.txt"

Remove-Item -LiteralPath $script:ReportPath,$AllManifest,$ZenodoManifest,$ZenodoList -Force -ErrorAction SilentlyContinue

$started = Get-Date
Write-Line "RUN_ALL integrity report"
Write-Line "Generated: $($started.ToString('yyyy-MM-dd HH:mm:ss zzz'))"
Write-Line "Result directory: $ResultDir"
Write-Line ""

$expectedStarts = @()
for ([Int64]$q = 100000000; $q -le 25950000000; $q += 50000000) {
    $expectedStarts += $q
}

$partGz = @(Get-ChildItem -LiteralPath $ResultDir -File -Filter "part_*.csv.gz")
$partStats = @(Get-ChildItem -LiteralPath $ResultDir -File -Filter "part_*.stats")
$partDone = @(Get-ChildItem -LiteralPath $ResultDir -File -Filter "part_*.done")

$ignoreNames = @(
    "RUN_ALL_INTEGRITY_REPORT.txt",
    "MANIFEST_SHA256_ALL.txt",
    "MANIFEST_SHA256_ZENODO.txt",
    "ZENODO_FILE_LIST.txt",
    "verify_result_integrity.ps1",
    "VERIFY_RESULT.bat",
    "README.txt"
)

$allFilesBefore = @(Get-ChildItem -LiteralPath $ResultDir -File |
    Where-Object { $_.Name -notin $ignoreNames })

$missing = New-Object System.Collections.Generic.List[string]
$unexpectedParts = New-Object System.Collections.Generic.List[string]

foreach ($q in $expectedStarts) {
    foreach ($suffix in @(".csv.gz",".stats",".done")) {
        $name = "part_${q}${suffix}"
        if (-not (Test-Path -LiteralPath (Join-Path $ResultDir $name))) {
            $missing.Add($name)
        }
    }
}

$expectedBaseNames = @{}
foreach ($q in $expectedStarts) {
    $expectedBaseNames["part_${q}.csv.gz"] = $true
    $expectedBaseNames["part_${q}.stats"] = $true
    $expectedBaseNames["part_${q}.done"] = $true
}

foreach ($f in @($partGz + $partStats + $partDone)) {
    if (-not $expectedBaseNames.ContainsKey($f.Name)) {
        $unexpectedParts.Add($f.Name)
    }
}

foreach ($n in @("pp.csv.gz","pp.stats","pp.done","summary.txt","failures.csv")) {
    if (-not (Test-Path -LiteralPath (Join-Path $ResultDir $n))) {
        $missing.Add($n)
    }
}

$failuresPath = Join-Path $ResultDir "failures.csv"
$failuresBytes = if (Test-Path -LiteralPath $failuresPath) {
    (Get-Item -LiteralPath $failuresPath).Length
} else { -1 }

Write-Line "STRUCTURE"
Write-Line "---------"
Write-Line "Expected prime-range parts : $($expectedStarts.Count)"
Write-Line "Found part_*.csv.gz        : $($partGz.Count)"
Write-Line "Found part_*.stats         : $($partStats.Count)"
Write-Line "Found part_*.done          : $($partDone.Count)"
Write-Line "Expected raw file count    : 1559"
Write-Line "Observed raw file count    : $($allFilesBefore.Count)"
Write-Line "failures.csv bytes         : $failuresBytes"
Write-Line "Missing expected files     : $($missing.Count)"
Write-Line "Unexpected part files      : $($unexpectedParts.Count)"

if ($missing.Count -gt 0) {
    Write-Line ""
    Write-Line "Missing:"
    foreach ($x in $missing) { Write-Line "  $x" }
}
if ($unexpectedParts.Count -gt 0) {
    Write-Line ""
    Write-Line "Unexpected part files:"
    foreach ($x in $unexpectedParts) { Write-Line "  $x" }
}

$summaryPath = Join-Path $ResultDir "summary.txt"
$summaryText = ""
if (Test-Path -LiteralPath $summaryPath) {
    $summaryText = Get-Content -LiteralPath $summaryPath -Raw
}

$summaryChecks = [ordered]@{
    "primes processed = 1127956365" = ($summaryText -match '1127956365')
    "prime powers = 13992" = ($summaryText -match '13992')
    "total fields = 1127970357" = ($summaryText -match '1127970357')
    "C4F = 1127970357" = ($summaryText -match 'C4F.*1127970357')
    "C4 = 0" = ($summaryText -match 'C4\s+.*:\s*0')
    "DS = 0" = ($summaryText -match 'DS\s+.*:\s*0')
    "UNRESOLVED = 0" = ($summaryText -match '(UNRESOLVED|other).*:\s*0')
    "failures lines = 0" = ($summaryText -match 'failures\.csv.*:\s*0')
}

Write-Line ""
Write-Line "SUMMARY CHECKS"
Write-Line "--------------"
foreach ($k in $summaryChecks.Keys) {
    $status = if ($summaryChecks[$k]) {"PASS"} else {"CHECK"}
    Write-Line ("{0,-38} : {1}" -f $k, $status)
}

$gzipFiles = @($partGz)
$ppGz = Join-Path $ResultDir "pp.csv.gz"
if (Test-Path -LiteralPath $ppGz) {
    $gzipFiles += Get-Item -LiteralPath $ppGz
}
$gzipFiles = @($gzipFiles | Sort-Object Name)

Write-Line ""
Write-Line "GZIP INTEGRITY"
Write-Line "--------------"
Write-Line "Files to test: $($gzipFiles.Count)"

$badGzip = New-Object System.Collections.Generic.List[string]
$i = 0
foreach ($f in $gzipFiles) {
    $i++
    Write-Progress -Activity "Testing gzip integrity" -Status "$i / $($gzipFiles.Count): $($f.Name)" -PercentComplete (($i * 100.0) / $gzipFiles.Count)
    if (-not (Test-GzipFile -Path $f.FullName)) {
        $badGzip.Add($f.Name)
    }
}
Write-Progress -Activity "Testing gzip integrity" -Completed
Write-Line "Corrupt/unreadable gzip files: $($badGzip.Count)"
foreach ($x in $badGzip) { Write-Line "  $x" }

$rawFiles = @($allFilesBefore | Sort-Object Name)
$zenodoFiles = @($rawFiles | Where-Object {
    $_.Name -eq "summary.txt" -or
    $_.Name -eq "failures.csv" -or
    $_.Name -eq "pp.csv.gz" -or
    $_.Name -eq "pp.stats" -or
    $_.Name -like "part_*.csv.gz" -or
    $_.Name -like "part_*.stats"
} | Sort-Object Name)

Write-Line ""
Write-Line "SHA-256 MANIFESTS"
Write-Line "----------------"
Write-Line "Hashing all raw files       : $($rawFiles.Count)"
Write-Line "Hashing Zenodo-stage files  : $($zenodoFiles.Count)"
Write-Line "Note: .done checkpoint files are excluded from the Zenodo-stage manifest."

$i = 0
foreach ($f in $rawFiles) {
    $i++
    Write-Progress -Activity "SHA-256: all raw files" -Status "$i / $($rawFiles.Count): $($f.Name)" -PercentComplete (($i * 100.0) / $rawFiles.Count)
    $hash = (Get-FileHash -LiteralPath $f.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    "{0}  {1}  {2}" -f $hash, $f.Length, $f.Name | Add-Content -LiteralPath $AllManifest -Encoding UTF8
}
Write-Progress -Activity "SHA-256: all raw files" -Completed

$i = 0
foreach ($f in $zenodoFiles) {
    $i++
    Write-Progress -Activity "SHA-256: Zenodo-stage files" -Status "$i / $($zenodoFiles.Count): $($f.Name)" -PercentComplete (($i * 100.0) / $zenodoFiles.Count)
    $hash = (Get-FileHash -LiteralPath $f.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    "{0}  {1}  {2}" -f $hash, $f.Length, $f.Name | Add-Content -LiteralPath $ZenodoManifest -Encoding UTF8
    $f.Name | Add-Content -LiteralPath $ZenodoList -Encoding UTF8
}
Write-Progress -Activity "SHA-256: Zenodo-stage files" -Completed

$structureOK = (
    $expectedStarts.Count -eq 518 -and
    $partGz.Count -eq 518 -and
    $partStats.Count -eq 518 -and
    $partDone.Count -eq 518 -and
    $allFilesBefore.Count -eq 1559 -and
    $missing.Count -eq 0 -and
    $unexpectedParts.Count -eq 0 -and
    $failuresBytes -eq 0
)
$summaryOK = -not ($summaryChecks.Values -contains $false)
$gzipOK = ($gzipFiles.Count -eq 519 -and $badGzip.Count -eq 0)

$finished = Get-Date
Write-Line ""
Write-Line "FINAL STATUS"
Write-Line "------------"
Write-Line "Structure                    : $(if ($structureOK) {'PASS'} else {'CHECK'})"
Write-Line "summary.txt expected counts  : $(if ($summaryOK) {'PASS'} else {'CHECK'})"
Write-Line "gzip integrity (519 files)   : $(if ($gzipOK) {'PASS'} else {'CHECK'})"
Write-Line "Manifest ALL entries         : $($rawFiles.Count)"
Write-Line "Manifest ZENODO entries      : $($zenodoFiles.Count)"
Write-Line "Started                      : $($started.ToString('yyyy-MM-dd HH:mm:ss'))"
Write-Line "Finished                     : $($finished.ToString('yyyy-MM-dd HH:mm:ss'))"
Write-Line "Elapsed                      : $([string]($finished - $started))"
Write-Line ""

if ($structureOK -and $summaryOK -and $gzipOK) {
    Write-Line "OVERALL: PASS"
} else {
    Write-Line "OVERALL: CHECK REQUIRED"
}

Write-Host ""
Write-Host "Finished."
Write-Host "Report: $script:ReportPath"
Write-Host "All-files manifest: $AllManifest"
Write-Host "Zenodo-stage manifest: $ZenodoManifest"
Write-Host "Zenodo file list: $ZenodoList"
Write-Host ""

if ($structureOK -and $summaryOK -and $gzipOK) {
    Write-Host "OVERALL: PASS" -ForegroundColor Green
} else {
    Write-Host "OVERALL: CHECK REQUIRED" -ForegroundColor Yellow
}
