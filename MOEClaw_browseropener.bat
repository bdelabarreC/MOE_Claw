@echo off
setlocal EnableDelayedExpansion

REM ============================================================================
REM  MOE-CLAW browser builder + server launcher
REM
REM  What it does:
REM    1. Activates the conda "rdkit" environment
REM    2. Runs the Python generator to turn moe_index.csv -> moe_browser.html
REM    3. Starts a local web server and opens the page in your browser
REM       (the server is required so Chrome will load the RDKit WASM for
REM        substructure search -- file:// double-click does NOT work)
REM
REM  Usage:
REM    run_moeclaw_browser.bat                 (uses the defaults below)
REM    run_moeclaw_browser.bat path\to\index.csv
REM ============================================================================

REM ---- CONFIG: adjust these three paths to match your machine -----------------
REM Folder that holds MOEClaw_buildbrowser_v7.py AND RDKit_minimal.js/.wasm.
REM %~dp0 means "the folder this .bat file lives in" -- keep the .bat in that
REM same folder and you won't need to edit this line.
set "APPDIR=%~dp0"

REM Name of the Python generator script (in APPDIR).
set "SCRIPT=MOEClaw_buildbrowser.py"

REM Conda environment name that has rdkit + pandas installed.
set "CONDA_ENV=rdkit"

REM Port for the local server.
set "PORT=8000"
REM ----------------------------------------------------------------------------

REM ---- Input CSV: first argument, or default to moe_index.csv in APPDIR -------
if "%~1"=="" (
    set "INPUT_CSV=%~dp0moe_files"
) else (
    set "INPUT_CSV=%~1"
)

REM Output HTML goes into APPDIR so it sits beside RDKit_minimal.js/.wasm.
set "OUTPUT_HTML=%APPDIR%MOEClaw_browser.html"

echo.
echo ==== MOE-CLAW browser builder ====
echo   App folder : %APPDIR%
echo   Input CSV  : %INPUT_CSV%
echo   Output     : %OUTPUT_HTML%
echo   Conda env  : %CONDA_ENV%
echo.

REM ---- Sanity checks ---------------------------------------------------------
if not exist "%INPUT_CSV%" (
    echo ERROR: input CSV not found:
    echo   %INPUT_CSV%
    echo Pass the path as an argument, or place moe_index.csv next to this .bat.
    goto :end
)
if not exist "%APPDIR%%SCRIPT%" (
    echo ERROR: generator script not found:
    echo   %APPDIR%%SCRIPT%
    goto :end
)

REM ---- Activate conda --------------------------------------------------------
REM "conda activate" only works in a shell that has been initialized. Inside a
REM .bat we call the conda hook first. We try a few common install locations.
set "CONDA_FOUND="
for %%D in (
    "%USERPROFILE%\anaconda3"
    "%USERPROFILE%\miniconda3"
    "%USERPROFILE%\Anaconda3"
    "%USERPROFILE%\Miniconda3"
    "%LOCALAPPDATA%\Continuum\anaconda3"
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
    echo WARNING: could not auto-locate your Anaconda/Miniconda install.
    echo Trying "conda" straight from PATH instead...
    call conda activate %CONDA_ENV%
) else (
    echo Using conda at: !CONDA_FOUND!
    call "!CONDA_FOUND!\Scripts\activate.bat" %CONDA_ENV%
)

if errorlevel 1 (
    echo ERROR: failed to activate conda environment "%CONDA_ENV%".
    echo Check the environment name, or activate it manually and run the
    echo python command shown below.
    goto :end
)



REM cd into APPDIR so the server's document root contains the HTML + JS + WASM.
cd /d "%APPDIR%"

REM Open the page in the default browser after a short delay so the server
REM has a moment to start. "start" returns immediately; the timeout gives slack.
start "" "http://localhost:%PORT%/MOEClaw_browser.html"

REM This call blocks and keeps the server running until you Ctrl+C.
python -m http.server %PORT%

:end
echo.
pause
endlocal