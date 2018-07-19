@echo off

if "%BIOIMG_ROOT%" == "" (
    set BIOIMG_ROOT=c:\bioimg
)
set PSPOT_FONT=%BIOIMG_ROOT%\share\pspot\font72.glf
set CYGWIN=nodosfilewarning
PATH=%BIOIMG_ROOT%\bin;%PATH%

if "%1" == "" (
    echo "Usage: spotview <input_file>"
    goto END
) 
%BIOIMG_ROOT%\bin\pspot.exe -v %1

:END
