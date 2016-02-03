#!/bin/bash

mkdir -P ./deps
git clone https://github.com/kramer314/fortran-lib.git ./deps/fortran-lib
cd ./deps/fortran-lib
scons
