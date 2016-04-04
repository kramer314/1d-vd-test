# Copyright (c) 2015 Alex Kramer <kramer.alex.kramer@gmail.com>
# See the LICENSE.txt file at the top-level directory of this distribution.

import os

source_dir = "#src/"
build_dir = "#build/"

libs = ["flib"]
lib_path= "#deps/fortran-lib/build"

env = DefaultEnvironment(ENV = os.environ, TOOLS = ['default', "gfortran"])

IEEE_flags = "-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans "
openmp_flags = "-fopenmp "
debug_flags = "-Og -g3 -Wall -Wextra -Wconversion -Wunused-parameter " + \
    "-pedantic -std=f2008 -fcheck=all -fbacktrace "
prod_flags = "-O3 -march=native "

flags = debug_flags + openmp_flags + IEEE_flags

env.Replace(F90FLAGS = flags)
env.Replace(LINKFLAGS = flags)
env.Replace(FORTRANMODDIRPREFIX = "-J ")
env.Replace(FORTRANMODDIR = build_dir)
env.Replace(F90PATH = lib_path)

Export("env")
Export("libs")
Export("lib_path")

SConscript(source_dir+"SConscript", variant_dir=build_dir, duplicate=1)

# For whatever reason, we can't use duplicate=0 and have *.mod files in the
# build directory. But, if we duplicate the source tree into the build
# directory SCons doesn't automatically clean the source files, so we have to
# manually define the entire build directory as a cleaning target.
Clean(".", build_dir)
