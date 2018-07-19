#!/bin/sh

if [ -z "$BIOIMG_ROOT" ]; then
    echo "The environment variable is not set: BIOIMG_ROOT"
    exit 0
fi

if [ $# = 0 ]; then
    inpf=inp1.pspot
else
    inpf=$1
fi
export PSPOT_FONT=$BIOIMG_ROOT/share/pspot/font72.glf

$BIOIMG_ROOT/bin/pspot -v $inpf
