Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
See the LICENSE.txt file at the top-level of this distribution.

Description
===========
Test of single-dimensional flux-based method for extracting momentum
distributions from wavepackets in time-dependent quantum mechanical
calculations. Both a semi-classical method and a method incorporating
first-order quantum corrections are implemented. Tests are done using the
free-particle propagation of a Gaussian wavepacket and the extracted momentum
spectra are compared with the analytically known results.

Example usage
=============
The following will run a simulation with the parameters given in the default
input input file `run.inp`:

    ./run.sh input/run.inp

After execution, the output directory in the input file (by default `output/`)
will have five updated files present:

* run.inp: Copy of input file
* run.log: Runtime logging
* vd_px.dat: Extracted momentum spectrum
* vd_resid.dat: Expected momentum spectrum, residuals of momentum spectrum, and
  cumulative residuals of momentum spectrum
* vd_resid_analysis.txt: Basic statistical analysis of results

Required dependencies
=====================
* Fortran compiler (most recently tested on GNU Fortran 5.3)
* SCons (most recently tested on SCons 2.4)

* Fortran generics library, available at:
  https://github.com/kramer314/fortran-lib
  Note that this library is still in beta status, so the API is technically
  still in flux. Always use the most recent version.

Recommended Dependencies
* OpenMP library (most recently tested on GCC OpenMP 5.3)

Basic setup
===========
The following will download and build this program and its library
dependencies:

    git clone https://github.com/kramer314/1d-vd-test.git
    cd ./1d-vd-test/
    ./build-deps.sh
    scons

Manual setup
============
Distribution of Fortran libraries is pretty nasty until submodules from F2015
are supported. What the `./build-deps.sh` script does is clone the Fortran
generics library project that this program uses into the directory `deps/`
and then builds it into the folder `deps/fortran-lib/build/`. You can have
this generics library project at some other location; you just need to build
the library as normal and create a symlink as follows:

    ln -s /path/to/fortran-lib-project ./deps/fortran-lib

This is a link to the project directory itself, not the directory where the
compiled library is actually located. Right now, the build script adds the
path `./deps/fortran-lib/build` to the compiler path. This may change.

Notes on `convergence-tests` directory
======================================
This is a collection of scripts that systematically batch together and run a
large number of simulations in parallel to test convergence sensitivity to grid
sizes. These are fairly hacked together, and use notation that isn't
documented anywhere but private notes.
