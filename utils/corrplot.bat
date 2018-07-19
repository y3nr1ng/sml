@echo off
if "%BIOIMG_ROOT%" == "" (
    set BIOIMG_ROOT=c:\bioimg
)
set XAPPLRESDIR=%BIOIMG_ROOT%\share\gnuplot\4.7\app-defaults
set GNUPLOT_PS_DIR=%BIOIMG_ROOT%\share\gnuplot\4.7\PostScript
set CYGWIN=nodosfilewarning
set PERL5LIB=%BIOIMG_ROOT%\lib\perl5\5.34

perl %BIOIMG_ROOT%\bin\evcorr xcor.dat
