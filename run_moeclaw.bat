@echo off
setlocal EnableDelayedExpansion

REM ============================================================================
REM  MOE-CLAW: build the browser from a CSV and serve it.
REM
REM    1. Activates the conda "moeclaw" environment
REM    2. Fetches RDKit-JS runtime the first time (if not already present)
REM    3. Builds moe_browser.html from your CSV
REM    4. Starts a local server and opens the page
REM       (the server is required so Chrome loads the WASM for substructure
REM        search -- opening the .html by double-click does NOT work)
REM
REM  Usage:
REM    run_moeclaw.bat                    (uses moe_index.csv in this folder)
REM    run_moeclaw.bat path\to\index.csv
REM ============================================================================

REM ---- CONFIG ----------------------------------------------------------------
set "APPDIR=%~dp0"
set "SCRIPT=MOEClaw_buildbrowser.py"
set "CONDA_ENV=rdkit"
set "PORT=8000"
REM ----------------------------------------------------------------------------

if "%~1"=="" (
    set "INPUT_CSV=%APPDIR%moe_index.csv"
) else (
    set "INPUT_CSV=%~1"
)
set "OUTPUT_HTML=%APPDIR%moe_browser.html"

echo.
echo ==== MOE-CLAW ====
echo   Folder    : %APPDIR%
echo   Input CSV : %INPUT_CSV%
echo   Conda env : %CONDA_ENV%
echo.

if not exist "%APPDIR%%SCRIPT%" (
    echo ERROR: %SCRIPT% not found in this folder.
    goto :end
)
if not exist "%INPUT_CSV%" (
    echo ERROR: input CSV not found: %INPUT_CSV%
    echo Pass a path as the first argument, or put moe_index.csv here.
    goto :end
)

REM ---- Activate conda --------------------------------------------------------
set "CONDA_FOUND="
for %%D in (
    "%USERPROFILE%\anaconda3"
    "%USERPROFILE%\miniconda3"
    "%USERPROFILE%\Anaconda3"
    "%USERPROFILE%\Miniconda3"
    "%ProgramData%\anaconda3"
    "%ProgramData%\miniconda3"
) do (
    if exist "%%~D\Scripts\activate.bat" (
        set "CONDA_FOUND=%%~D"
        goto :got_conda
    )
)
:got_conda
if not defined CONDA_FOUND (
    echo Could not auto-locate Anaconda/Miniconda; trying conda from PATH...
    call conda activate %CONDA_ENV%
) else (
    echo Using conda at: !CONDA_FOUND!
    call "!CONDA_FOUND!\Scripts\activate.bat" %CONDA_ENV%
)
if errorlevel 1 (
    echo ERROR: failed to activate conda env "%CONDA_ENV%".
    echo Create it first:  conda env create -f environment.yml
    goto :end
)

REM ---- Fetch RDKit-JS runtime if missing -------------------------------------
@REM  if not exist "%APPDIR%RDKit_minimal.wasm" (
@REM      echo.
@REM      echo First run: fetching RDKit-JS runtime ^(~7 MB^)...
@REM      python "%APPDIR%fetch_rdkit_js.py" --dir "%APPDIR%"
@REM      if errorlevel 1 (
@REM          echo ERROR: could not fetch RDKit-JS. Check your internet connection.
@REM          goto :end
@REM      )
@REM  )

REM ---- Build -----------------------------------------------------------------
echo.
echo Building HTML from CSV...
python "%APPDIR%%SCRIPT%" "%INPUT_CSV%" -o "%OUTPUT_HTML%"
if errorlevel 1 (
    echo ERROR: the Python generator failed. See messages above.
    goto :end
)
echo Build complete: %OUTPUT_HTML%

REM ---- Serve + open ----------------------------------------------------------
echo.
echo Starting local server on http://localhost:%PORT%/
echo   Leave this window open while you use the browser.
echo   Press Ctrl+C here to stop the server when done.
echo.
cd /d "%APPDIR%"
start "" "http://localhost:%PORT%/moe_browser.html"
python -m http.server %PORT%

:end
echo.
pause
endlocal
