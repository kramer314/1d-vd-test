#!/bin/bash

export OMP_NUM_THREADS=1

# $1 is the input file

./build/main.x $1
