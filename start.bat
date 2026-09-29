@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0start.ps1" %*
if errorlevel 1 (
    echo.
    echo Molecule Studio V3 stopped with an error. Read the message above.
    pause
    exit /b 1
)
