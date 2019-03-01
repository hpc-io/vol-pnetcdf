#ifndef ADIOS_VOL_H
#define ADIOS_VOL_H

//#define H5VL_FRIEND // junmin added
#define H5F_FRIEND /*suppress error about including H5Fpkg   */

#include "H5Fpkg.h"
#include "hdf5.h"

#include "mpi.h"
#include "reader.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ADIOS 511
#define VOLNAME "ADIOS-VOL"

#define ADIOS_TIMESTEP "ADIOS_TIMESTEP"

extern hid_t GLOBAL_FILE_ID;

typedef struct H5VL_complex_float_t {
  float real;
  float img;
} H5VLT_complex_float_t;

typedef struct H5VL_complex_double_t {
  double real;
  double img;
} H5VLT_complex_double_t;

typedef struct H5VL_adios_dim_t {
  hsize_t ndim;
  hsize_t *dims;
} H5VL_adios_dim_t;

typedef struct H5VL_adios_var_t {
  char *name;
  unsigned int ntimestep;

  H5VL_adios_dim_t *dimInfo;

  unsigned int curr_timestep;
  unsigned int is_global_dim;

  H5_bp_file_t *fileReader;

  // hid_t assigned_gid;
} H5VL_adios_var_t;

typedef struct H5VL_adios_attr_t {
  char *name;
  // H5VL_adios_dim_t *dimInfo;
  int ts; // only needed forthe ADIOS_TIMESTEP
  H5_bp_file_t *fileReader;
} H5VL_adios_attr_t;

typedef struct H5VL_adios_t {
  H5_bp_file_t *bpFileReader;
} H5VL_adios_t;

typedef struct H5VL_link_adios_varnames_t {
  hid_t fid;
  int rank;
  int size;
  char **varNames;
} H5VL_adios_varnames_t;

typedef struct {
  hsize_t posCompact;  // e.g. 0,1,2,3
  hsize_t posInSource; // e.g. 3,5,8.9
} H5_posMap;

void GetSelOrder(hid_t space_id, H5S_sel_type space_type, H5_posMap **result);

void assignToMemSpace(H5_posMap *sourceSelOrderInC,
                      H5_posMap *targetSelOrderInC, hsize_t total,
                      size_t dataTypeSize, char *adiosData, char *buf);

herr_t H5Pset_fapl_adios(hid_t *acc_tpl, hid_t *fapl_id);
herr_t H5P_unset_adios();

#endif
