@echo off
setlocal
cd /d "%~dp0"
echo.
echo RUN_ALL result integrity check
echo ==============================
echo.
echo This can take some time because every .csv.gz file is fully tested
echo and every raw output file is hashed with SHA-256.
echo.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0verify_result_integrity.ps1"
echo.
pause
