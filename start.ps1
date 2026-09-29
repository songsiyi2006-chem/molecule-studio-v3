param(
    [ValidateRange(1, 65535)][int]$Port = 8003,
    [switch]$NoBrowser,
    [switch]$Train,
    [switch]$Test,
    [switch]$PrepareData
)

$ErrorActionPreference = 'Stop'
$projectRoot = $PSScriptRoot
$localPython = Join-Path $projectRoot '.venv\Scripts\python.exe'
if (Test-Path -LiteralPath $localPython -PathType Leaf) {
    $projectPython = $localPython
} else {
    $pythonCommand = Get-Command python -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
    if (-not $pythonCommand) {
        throw 'Python was not found. Install Python 3.12, create .venv, and install requirements.txt.'
    }
    $projectPython = $pythonCommand.Source
}

$selectedActions = @($Train, $Test, $PrepareData) | Where-Object { $_.IsPresent }
if ($selectedActions.Count -gt 1) {
    throw 'Choose only one of -Train, -Test, or -PrepareData.'
}

# A venv may reuse an existing Conda runtime. Resolve its actual base prefix
# instead of storing a machine-specific installation path in the project.
$prefixOutput = & $projectPython -c 'import sys; print(sys.base_prefix)'
if ($LASTEXITCODE -ne 0 -or -not $prefixOutput) {
    throw 'The Python interpreter could not start. Recreate the project environment.'
}
$runtimePrefix = ($prefixOutput | Select-Object -Last 1).Trim()
$runtimePaths = @(
    @(
        (Join-Path $runtimePrefix 'Library\bin'),
        (Join-Path $runtimePrefix 'Scripts'),
        $runtimePrefix
    ) | Where-Object { Test-Path -LiteralPath $_ -PathType Container }
)
$env:PATH = (($runtimePaths + @($env:PATH)) -join [System.IO.Path]::PathSeparator)
$env:OMP_NUM_THREADS = '2'
$env:MKL_NUM_THREADS = '2'
$env:OPENBLAS_NUM_THREADS = '2'
$env:PYTHONUTF8 = '1'

Push-Location -LiteralPath $projectRoot
try {
    if ($Train) {
        Write-Host 'Training all three physicochemical models using the project environment...'
        & $projectPython 'scripts\train_models.py'
    } elseif ($Test) {
        Write-Host 'Checking the supplied project and model artifacts...'
        & $projectPython -m unittest discover -s tests -v
    } elseif ($PrepareData) {
        Write-Host 'Rebuilding the database from verified local raw snapshots (offline)...'
        & $projectPython 'scripts\prepare_data.py' --offline
        if ($LASTEXITCODE -eq 0) {
            & $projectPython 'scripts\build_catalog.py'
        }
    } else {
        Write-Host ('Molecule Studio V3: http://127.0.0.1:{0}' -f $Port)
        if ($NoBrowser) {
            & $projectPython 'run.py' --port $Port
        } else {
            & $projectPython 'run.py' --port $Port --open
        }
    }
    $projectExitCode = $LASTEXITCODE
} finally {
    Pop-Location
}
exit $projectExitCode
