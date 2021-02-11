#ifndef DEFS_HPP_
#define DEFS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file defs.hpp.in
//  \brief Template file for defs.hpp.  When the configure.py script is run, a new
//  defs.hpp file will be created (overwriting the last) from this template.  This new
//  file contains Athena++ specific cpp macros and definitions set by configure.

//----------------------------------------------------------------------------------------
// macros which define physics and algorithms

// problem generator
#define PROBLEM_GENERATOR "globAccDisk"

// coordinate system
#define COORDINATE_SYSTEM "spherical_polar"

// enable shearing box? default=0 (false)
#define SHEARING_BOX 0

// non-barotropic equation of state (i.e. P not simply a func of rho)? default=1 (true)
#define NON_BAROTROPIC_EOS 1

// Riemann solver
#define RIEMANN_SOLVER "roe"

// include magnetic fields? default=0 (false)
#define MAGNETIC_FIELDS_ENABLED 1

// include self gravity? default=0 (false)
#define SELF_GRAVITY_ENABLED 0

// make use of FFT? default=0 (false)
#define FFT_ENABLED 0

#define NO_FFT

// include radiative transfer? default=0 (false)
#define RADIATION_ENABLED 0

// enable special or general relativity? default=0 (false)
#define RELATIVISTIC_DYNAMICS 0

// enable general relativity? default=0 (false)
#define GENERAL_RELATIVITY 0

// enable GR frame transformations? default=0 (false)
#define FRAME_TRANSFORMATIONS 0

// floating point single precision option
#define SINGLE_PRECISION_ENABLED 0

// MPI parallelization (MPI_PARALLEL or NOT_MPI_PARALLEL)
#define NOT_MPI_PARALLEL

// openMP parallelization (OPENMP_PARALLEL or NOT_OPENMP_PARALLEL)
#define NOT_OPENMP_PARALLEL

// HDF5 output (HDF5OUTPUT or NO_HDF5OUTPUT)
#define HDF5OUTPUT

// HDF5 floating point ouptut precision option
#define H5_DOUBLE_PRECISION_ENABLED 0

// compiler options
#define COMPILED_WITH "g++"
#define COMPILER_COMMAND "mpicxx"
#define COMPILED_WITH_OPTIONS " -I/usr/local/Cellar/hdf5/1.10.4/include -O3 -std=c++11  -L/usr/local/Cellar/hdf5/1.10.4/lib  -lhdf5"

//----------------------------------------------------------------------------------------
// macros associated with numerical algorithm (rarely modified)

#define NHYDRO 5
#define NFIELD 3
#define NWAVE 7
#define NGHOST 2
#define MAX_NSTAGE 5     // maximum number of stages per cycle for time-integrator
#define MAX_NREGISTER 3  // maximum number of (u, b) register pairs for time-integrator

//----------------------------------------------------------------------------------------
// general purpose macros (never modified)

#define PI 3.1415926535897932
#define SQRT2 1.4142135623730951
#define ONE_OVER_SQRT2 0.7071067811865475
#define ONE_3RD 0.3333333333333333
#define TWO_3RD 0.6666666666666667
#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+36
#define SQR(x) ( (x)*(x) )
#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 )

#endif // DEFS_HPP_
