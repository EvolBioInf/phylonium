#!/bin/sh -f

ARG1="$1"
PHYLONIUM=${ARG1:-./src/phylonium}
SIMF=./test/simf

# Test if the tools exists, and can be executed
$PHYLONIUM --version > /dev/null || exit 1
$SIMF -h > /dev/null || exit 1

SEED=${RANDOM_SEED:-0}
SEED2=0
SEED3=0
if test $SEED -ne 0; then
    SEED=$((SEED + 1))
    SEED2=$((SEED + 2))
    SEED3=$((SEED + 3))
fi

$SIMF -s $SEED -l 100000 -p simple
$PHYLONIUM simple0.fasta simple1.fasta > /dev/null || exit 1

rm simple0.fasta simple1.fasta
