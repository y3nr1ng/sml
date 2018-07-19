#!/bin/sh

if [ -z "$BIOIMG_ROOT" ]; then
    echo "The environment variable is not set: BIOIMG_ROOT"
    exit 0
fi

if [ $# = 0 ]; then
    inpf=input
else
    inpf=$1
fi
export PSPOT_FONT=$BIOIMG_ROOT/share/pspot/font72.glf

$BIOIMG_ROOT/bin/geninp $inpf
$BIOIMG_ROOT/bin/pix    inp.pix
$BIOIMG_ROOT/bin/pclst  inp.cl  > log.cl
$BIOIMG_ROOT/bin/evcorr xcor.dat
$BIOIMG_ROOT/bin/pspot  inp2.pspot
$BIOIMG_ROOT/bin/pspot  inp3.pspot
$BIOIMG_ROOT/bin/pspot -v inp1.pspot
