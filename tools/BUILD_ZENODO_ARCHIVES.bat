@echo off
setlocal
cd /d "%~dp0"
echo.
echo Building Zenodo archives from this completed result folder...
echo.
python "%~dp0build_zenodo_archives.py"
echo.
pause
