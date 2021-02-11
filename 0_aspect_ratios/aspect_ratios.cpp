/*

A program to read from Athena++ hdf5 and calculate the minimum and maximum aspect ratios of the mesh cells.

By: Patryk Pjanka, unless stated otherwise.

WARNING: Only spherical_polar coordinates implemented!

*/

#define VERBOSE

#include <fstream>
#include <vector>
#include <utility> // pair

#include "aspect_ratios.hpp"

using namespace std;

const double adiab_idx = 1.1;

int main (int argc, char** argv) {

  char* filename = argv[1];
  int mb,i,j,k;

  fstream outfile;
  if (argc > 2) {
    outfile.open(argv[2], ios::out);
    outfile << "x1f, x2f, x3f, dx1, dx2, dx3, dx1/dx2, dx2/dx3, dx3/dx1, dt1, dt2, dt3, dts1, dts2, dts3, dts, dtms, dt, vel1, vel2, vel3, bcc1, bcc2, bcc3, valfven, rho" << endl;
    outfile << scientific;
  }

  // read the shape of the data
  int* data_shape = read_metadata(filename, "prim");
  for (int i = 0; i <= data_shape[0]+1; i++)
    cout << data_shape[i] << "  ";
  cout << endl;

  // read the coordinate values
  AthenaArray<double>* x1f = read_xnf (filename, "x1f");
  AthenaArray<double>* x2f = read_xnf (filename, "x2f");
  AthenaArray<double>* x3f = read_xnf (filename, "x3f");
  AthenaArray<double>* rho = read_quantity (filename, "prim", -1, 0);
  AthenaArray<double>* press = read_quantity (filename, "prim", -1, 1);
  AthenaArray<double>* vel1 = read_quantity (filename, "prim", -1, 2);
  AthenaArray<double>* vel2 = read_quantity (filename, "prim", -1, 3);
  AthenaArray<double>* vel3 = read_quantity (filename, "prim", -1, 4);
  AthenaArray<double>* bcc1 = read_quantity (filename, "B", -1, 0);
  AthenaArray<double>* bcc2 = read_quantity (filename, "B", -1, 1);
  AthenaArray<double>* bcc3 = read_quantity (filename, "B", -1, 2);

  // now go through each meshblock and check the aspect ratios
  vector< pair<double,double> > aspect_ratios; // [1/2, 2/3, 3/1]
  for (i = 0; i < 3; i++)
    aspect_ratios.push_back(pair<double,double>(1.e6,0.));
  double x1, x2, x3;
  double dx1, dx2, dx3;
  double dt1, dt2, dt3, dts1, dts2, dts3, dts, dt, dtms;
  double one_two, two_three, three_one;
  double csound;
  double bcc, valfven;
  for (mb = 0; mb < data_shape[2]; mb++) { // meshblock
    k = 0; // unlikely to depend on phi
    //for (k = 0; k < data_shape[3]-1; k++) {  // phi
      for (j = 0; j < data_shape[4]-1; j++) { // theta
        for (i = 0; i < data_shape[5]-1; i++) { // r
          // calculate roughly mean mid-cell dimensions
          x1 = 0.5 * ((*x1f)(mb,i+1) + (*x1f)(mb,i));
          x2 = 0.5 * ((*x2f)(mb,j+1) + (*x2f)(mb,j));
          x3 = 0.5 * ((*x3f)(mb,k+1) + (*x3f)(mb,k));
          // calculate cell size
          dx1 = (*x1f)(mb,i+1) - (*x1f)(mb,i);
          dx2 = x1 * ((*x2f)(mb,j+1) - (*x2f)(mb,j));
          dx3 = x1 * sin(x2) * ((*x3f)(mb,k+1) - (*x3f)(mb,k));
          // calculate ratios
          one_two = dx1/dx2;
          two_three = dx2/dx3;
          three_one = dx3/dx1;
          // calculate fluid speed time steps
          dt1 = fabs(dx1 / (*vel1)(mb, k,j,i));
          dt2 = fabs(dx2 / (*vel2)(mb, k,j,i));
          dt3 = fabs(dx3 / (*vel3)(mb, k,j,i));
          // calculate sound speed time step
          csound = sqrt(adiab_idx * (*press)(mb, k,j,i)/(*rho)(mb, k,j,i));
          dts1 = dx1 / csound;
          dts2 = dx2 / csound;
          dts3 = dx3 / csound;
          dts = min(min(dx1, dx2), dx3) / csound;
          // calculate minimal magnetosonic time step
          bcc = /*sqrt(4.*M_PI) */ sqrt((*bcc1)(mb, k,j,i)*(*bcc1)(mb, k,j,i) + (*bcc2)(mb, k,j,i)*(*bcc2)(mb, k,j,i) + (*bcc3)(mb, k,j,i)*(*bcc3)(mb, k,j,i));
          valfven = bcc / sqrt(/*4.*M_PI */ (*rho)(mb, k,j,i));
          dtms = min(min(dx1, dx2), dx3) / sqrt(valfven*valfven + csound*csound);
          // final time step
          dt = min(min(min(min(dt1, dt2), dt3), dts), dtms);
          // report to file if requested
          if (argc > 2)
            outfile << x1 << ", " << x2 << ", " << x3 << ", "
               << dx1 << ", " << dx2 << ", " << dx3 << ", "
               << one_two << ", "  << two_three << ", "  << three_one << ", " 
               << dt1 << ", " << dt2 << ", " << dt3 << ", " 
               << dts1 << ", " << dts2 << ", " << dts3 << ", " 
               << dts << ", " << dtms << ", " << dt << ", " 
               << fabs((*vel1)(mb, k,j,i)) << ", " << fabs((*vel2)(mb, k,j,i)) << ", " << fabs((*vel3)(mb, k,j,i)) << ", " 
               << fabs((*bcc1)(mb, k,j,i)) << ", " << fabs((*bcc2)(mb, k,j,i)) << ", " << fabs((*bcc3)(mb, k,j,i)) << ", "
               << valfven << ", "
               << fabs((*rho)(mb, k,j,i)) << endl;
          // update minmax if needed
          if (one_two < aspect_ratios[0].first)
            aspect_ratios[0].first = one_two;
          else if (one_two > aspect_ratios[0].second)
            aspect_ratios[0].second = one_two;
          if (two_three < aspect_ratios[1].first)
            aspect_ratios[1].first = two_three;
          else if (two_three > aspect_ratios[1].second)
            aspect_ratios[1].second = two_three;
          if (three_one < aspect_ratios[2].first)
            aspect_ratios[2].first = three_one;
          else if (three_one > aspect_ratios[2].second)
            aspect_ratios[2].second = three_one;
        }
      }
    //}
  }

  // print the results
  for (i = 0; i < 3; i++)
    cout << "dx" << i+1 << " / dx" << (i+1)%3+1 << " ratio: min " << aspect_ratios[i].first << ", max " << aspect_ratios[i].second << endl;

  if (argc > 2)
    outfile.close();

  return 0;
}

/*
The following is based on the example file supplied by the HDF5 team:
https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html
Returns for dataset_name = prim/cons: <File rank> <number of quantities if more than one> <# of meshblocks> <meshblock size: 3 values> <total no of points>
Returns for dataset_name = xnf: <File rank> <# of meshblocks> <meshblock size> <total no of points>
*/
// check the file details
int* read_metadata (char* filename, char* dataset_name) {

  hid_t file, dataset; /* handles */
  hid_t       datatype, dataspace;
  hid_t       memspace;
  H5T_class_t t_class;                 /* data type class */
  H5T_order_t order;                 /* data order */
  size_t      size;                  /*
              * size of the data element
              * stored in file
              */
  hsize_t     dimsm[3];              /* memory space dimensions */
  hssize_t    npoints;               /* number of points in the dataset */
  herr_t      status;

  hsize_t      count[2];              /* size of the hyperslab in the file */
  hsize_t      offset[2];             /* hyperslab offset in the file */
  hsize_t      count_out[3];          /* size of the hyperslab in memory */
  hsize_t      offset_out[3];         /* hyperslab offset in memory */
  int          i, j, k, status_n, rank;

  /*
   * Open the file and the dataset.
   */
  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);

  /*
   * Get datatype and dataspace handles and then query
   * dataset class, order, size, rank and dimensions.
   */
  datatype  = H5Dget_type(dataset);     /* datatype handle */
  t_class     = H5Tget_class(datatype);
  order     = H5Tget_order(datatype);
  size  = H5Tget_size(datatype);
  #ifdef VERBOSE
  if (t_class == H5T_INTEGER) printf("Data set has INTEGER type \n");
  if (order == H5T_ORDER_LE) printf("Little endian order \n");
  printf(" Data size is %d \n", (int)size);
  #endif


  dataspace = H5Dget_space(dataset);    /* dataspace handle */
  rank      = H5Sget_simple_extent_ndims(dataspace);
  hsize_t     dims_out[rank];           /* dataset dimensions */
  status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  npoints   = H5Sget_simple_extent_npoints(dataspace);

  // print out info

  #ifdef VERBOSE
  cout << " rank " << rank << ", ";
  for (i = 0; i < rank-1; i++)
    cout << (unsigned long)(dims_out[i]) << " x ";
  cout << (unsigned long)(dims_out[rank-1]) << ", no of points " << (unsigned long)(npoints) << endl;
  #endif

  // output the details needed to initialize the AthenaArray
  int* result = new int [rank+2]; // rank, dimensions, npoints
  result[0] = rank;
  for (i = 0; i < rank; i++) result[i+1] = (unsigned long)(dims_out[i]);
  result[rank+1] = (unsigned long)(npoints);
  return result;
}

// read the coordinate grid
AthenaArray<double>* read_xnf (char* filename, char* dataset_name, int meshblock_no) {

  #ifdef VERBOSE
  cout << "Reading " << dataset_name << " from " << filename << ", meshblock " << meshblock_no << "... " << endl;
  #endif

  int i,j,k;

  // check the file details
  int* file_shape = read_metadata(filename, dataset_name);

  // select from file
  int rank_file = file_shape[0];
  int rank_mem = rank_file;
  int start_file [rank_file];
  int count_file [rank_file];
  int start_mem [rank_mem];
  int count_mem [rank_mem];
  AthenaArray<Real>* array = new AthenaArray<Real>;

  // meshblocks to read
  if (meshblock_no < 0) {
    start_file[0] = start_mem[0] = 0;
    count_file[0] = count_mem[0] = file_shape[1];
    array->NewAthenaArray(file_shape[1], file_shape[2]);
  } else {
    start_file[0] = meshblock_no;
    start_mem[0] = 0;
    count_file[0] = count_mem[0] = 1;
    array->NewAthenaArray(1, file_shape[2]);
  }

  // read the whole meshblock
  start_file[1] = start_mem[1] = 0;
  count_file[1] = count_mem[1] = file_shape[2];

  HDF5ReadRealArray(filename, dataset_name, rank_file, start_file, count_file, rank_mem, start_mem, count_mem, *array, true, false);

  #ifdef VERBOSE
  cout << "done." << endl;
  #endif

  return array;
}

// read a physical quantity
AthenaArray<double>* read_quantity (char* filename, char* dataset_name, int meshblock_no, int quantity) {

  #ifdef VERBOSE
  cout << "Reading " << dataset_name << " from " << filename << ", meshblock " << meshblock_no << "... " << endl;
  #endif

  int i;

  // check the file details
  int* file_shape = read_metadata(filename, dataset_name);

  // select from file
  int rank_file = file_shape[0];
  int rank_mem = rank_file;
  int start_file [rank_file];
  int count_file [rank_file];
  int start_mem [rank_mem];
  int count_mem [rank_mem];
  int array_shape [rank_file]; // stores no of physical quantities and meshblocks read
  AthenaArray<Real>* array = new AthenaArray<Real>;

  // quantities to read
  if (quantity < 0) {
    start_file[0] = start_mem[0] = 0;
    count_file[0] = count_mem[0] = file_shape[1];
    array_shape[0] = file_shape[1];
  } else {
    start_file[0] = quantity;
    start_mem[0] = 0;
    count_file[0] = count_mem[0] = 1;
    array_shape[0] = 1;
  }

  // meshblocks to read
  if (meshblock_no < 0) {
    start_file[1] = start_mem[1] = 0;
    count_file[1] = count_mem[1] = file_shape[2];
    array_shape[1] = file_shape[2];
  } else {
    start_file[1] = meshblock_no;
    start_mem[1] = 0;
    count_file[1] = count_mem[1] = 1;
    array_shape[1] = 1;
  }

  // read the whole meshblock
  for (i = 0; i < rank_file-2; i++) {
    start_file[i+2] = start_mem[i+2] = 0;
    count_file[i+2] = count_mem[i+2] = file_shape[i+3];
    array_shape[i+2] = file_shape[i+3];
  }

  delete file_shape;

  array->NewAthenaArray(array_shape[0], array_shape[1], array_shape[2], array_shape[3], array_shape[4]);
  HDF5ReadRealArray(filename, dataset_name, rank_file, start_file, count_file, rank_mem, start_mem, count_mem, *array, true, false);

  #ifdef VERBOSE
  cout << "done." << endl;
  #endif

  return array;
}