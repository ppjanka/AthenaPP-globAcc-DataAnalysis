#ifndef ASPECT_RATIOS_HPP
#define ASPECT_RATIOS_HPP

/*

Read from Athena++ hdf5 and calculate the minimum and maximum aspect ratios of the mesh cells.

By: Patryk Pjanka, unless stated otherwise.

*/

#import<iostream>
#include "hdf5.h"

#include "athena_hdf5_reader/athena.hpp"         // Real
#include "athena_hdf5_reader/athena_arrays.hpp"  // AthenaArray

// declarations
void HDF5ReadRealArray(const char *filename, const char *dataset_name, int rank_file,
    const int *start_file, const int *count_file, int rank_mem, const int *start_mem,
    const int *count_mem, AthenaArray<Real> &array, bool collective, bool noop);

int* read_metadata (char* filename, char* dataset_name);

AthenaArray<double>* read_xnf (char* filename, char* dataset_name, int meshblock_no=-1);
AthenaArray<double>* read_quantity (char* filename, char* dataset_name, int meshblock_no=-1, int quantity=-1);

#endif