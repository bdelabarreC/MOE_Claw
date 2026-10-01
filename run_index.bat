@echo off
setlocal

REM ============================================================================
REM  MOE-CLAW indexer -- runs moe_claw.svl under MOE (headless) to turn a
REM  folder of .moe files into moe_index.csv.
REM
REM  Shows progress on screen AND writes a log.
REM ============================================================================

REM --- Adjust this to your MOE install path ---
set MOE_BIN=D:\Program Files\moe2024.0601\bin
REM --- Adjust to location of moe_claw script
set SCRIPT=%~dp0moe_claw.svl
REM --- Adjust to where you want log output deposited (helps with trouble shooting)
set LOGDIR=%~dp0logs

REM --- Directory to index: use first argument, or a default as shown below ---
if "%~1"=="" (
    set TARGET_DIR=%~dp0moe_files
) else (
    set TARGET_DIR=%~1
)

REM --- Derive CSV path from the target dir ---
set OUT_CSV=%TARGET_DIR%/moe_index.csv

if not exist "%LOGDIR%" mkdir "%LOGDIR%"
set LOGFILE=%LOGDIR%\moe_index_%date:~-4%%date:~4,2%%date:~7,2%.log

REM --- Tell the user what's happening (these go to the SCREEN, no redirect) ---
echo.
echo ============================================
echo   MOE-CLAW indexer
echo ============================================
echo   Indexing : %TARGET_DIR%
echo   Output   : %OUT_CSV%
echo   Log      : %LOGFILE%
echo.

REM --- Sanity checks ----------------------------------------------------------
if not exist "%MOE_BIN%\moebatch.exe" (
    echo ERROR: moebatch.exe not found at:
    echo   %MOE_BIN%
    goto :end
)
if not exist "%SCRIPT%" (
    echo ERROR: SVL script not found at:
    echo   %SCRIPT%
    goto :end
)

echo Starting MOE. Reading .moe files can take SEVERAL MINUTES for a large
echo folder -- each file must be fully loaded before it can be indexed.
echo Progress appears below. Do not press Ctrl+C unless you want to abort.
echo.
echo [%date% %time%] Indexing %TARGET_DIR% >> "%LOGFILE%"

REM --- Run MOE: show output live on screen AND append it to the log -----------
REM (Batch has no native 'tee', so we pipe through PowerShell's Tee-Object.)
"%MOE_BIN%\moebatch.exe" -licwait -run "%SCRIPT%" -dir "%TARGET_DIR%" -csv "%OUT_CSV%" 2>&1 | powershell -NoProfile -Command "$input | Tee-Object -FilePath '%LOGFILE%' -Append"

echo [%date% %time%] Exit code %ERRORLEVEL% >> "%LOGFILE%"
echo.
echo ============================================
echo   Finished. CSV written to:
echo   %OUT_CSV%
echo ============================================

:end
echo.
pause
endlocal
