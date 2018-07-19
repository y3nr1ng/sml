@echo off
if "%BIOIMG_ROOT%" == "" (
    set BIOIMG_ROOT=c:\bioimg
)
if "%1" == "" (
    set inpf=input
) else (
    set inpf=%1
)
set XAPPLRESDIR=%BIOIMG_ROOT%\share\gnuplot\4.7\app-defaults
set GNUPLOT_PS_DIR=%BIOIMG_ROOT%\share\gnuplot\4.7\PostScript
set CYGWIN=nodosfilewarning
set PERL5LIB=%BIOIMG_ROOT%\lib\perl5\5.34
set PSPOT_FONT=%BIOIMG_ROOT%\share\pspot\font72.glf

perl %BIOIMG_ROOT%\bin\geninp %inpf%
pix inp.pix
pclst inp.cl > log.cl
perl %BIOIMG_ROOT%\bin\evcorr xcor.dat
pspot inp2.pspot
pspot inp3.pspot
pspot -v inp1.pspot
