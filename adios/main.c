#include "H5Epublic.h"
#include "adios_vol.h"
#include "hdf5.h"
#include "mpi.h"

hid_t GLOBAL_FILE_ID = 0;

void browseVar(hid_t file_id, const char *varName);
void browseVarGroup(hid_t file_id, const char *varName);

herr_t iter_callback_func(hid_t loc_id, const char *name,
                          const H5L_info_t *info, void *op_data) {
  herr_t status, return_val = 0;

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);
  /*
  if (strlen(name) != 18) {
    if (rank == 0)
      printf("\n skipping .. %s \n", name);
    return return_val;
  }
  */
  // printf("\n=====> handled in callback ...name=%s len=%lu\n", name,
  // strlen(name));
  printf("\n=====> handled in callback ");

  H5VL_adios_varnames_t *data = (H5VL_adios_varnames_t *)op_data;
  browseVarGroup(data->fid, name);
  // browseVar(data->fid, name);
  return return_val;
}

void printBBox(hsize_t *count, hsize_t *offset, int ndim) {
  int i;
  printf("\t bbox size:: ");
  for (i = 0; i < ndim; i++) {
    printf("[%llu + %llu]", offset[i], count[i]);
  }
  printf("\n");
}

void readData(hid_t dataType, hid_t memspace, hid_t dspace, hsize_t size,
              hid_t dset) {
  H5T_class_t typeClass = H5Tget_class(dataType);
  size_t s = H5Tget_size(dataType);

  hsize_t memSpaceSize = 1;
  int n;

  if (memspace == H5S_ALL) {
    memSpaceSize = size;
  } else {
    const int ndims = H5Sget_simple_extent_ndims(memspace);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(memspace, dims, NULL);

    for (n = 0; n < ndims; n++) {
      memSpaceSize *= dims[n];
    }
  }

  printf("\tmem space base size = %llu, selection size=%llu\n", memSpaceSize,
         size);
  if (typeClass == H5T_FLOAT) {
    double data_out[memSpaceSize];
    for (n = 0; n < memSpaceSize; n++) {
      data_out[n] = (double)0;
    }
    herr_t status =
        H5Dread(dset, dataType, memspace, dspace, H5P_DEFAULT, data_out);

    if (size > 1) {
      printf("\t... first two values: %.3f, %.3f\n", data_out[0], data_out[1]);
      printf("\t... last two values: %.3f, %.3f\n", data_out[size - 2],
             data_out[size - 1]);
    } else {
      printf("\t... value=: %.3f\n", data_out[0]);
    }
  } else if (typeClass == H5T_BITFIELD) {
    printf(" ......... print char?? ");
  } else if (typeClass == H5T_INTEGER) {
    if (s == 1) {
      // read int8_t
      int8_t data_out[memSpaceSize];
      herr_t status =
          H5Dread(dset, dataType, memspace, dspace, H5P_DEFAULT, data_out);
      if (size > 1) {
        // printf("... first two values: %" PRIi8, data_out[0]);
        // printf(",%" PRIi8 "\n", data_out[1]);
        printf("... first two values: %" PRIi8 ", %" PRIi8 "\n", data_out[0],
               data_out[1]);
        printf("... last two values: %d, %d\n", data_out[size - 2],
               data_out[size - 1]);
      } else {
        printf("... value=: %d\n", data_out[0]);
      }

    } else if (s == 4) { // similarly for unsigned
      // read int
      int data_out[memSpaceSize];
      herr_t status =
          H5Dread(dset, dataType, memspace, dspace, H5P_DEFAULT, data_out);
      if (size > 1) {
        printf("... first two values: %d , %d\n", data_out[0], data_out[1]);
        printf("... last two values: %d, %d\n", data_out[size - 2],
               data_out[size - 1]);
      } else {
        printf("... value=: %d\n", data_out[0]);
      }
    } else if (s == 8) { // similarly for unsigned
      int64_t data_out[memSpaceSize];
      herr_t status =
          H5Dread(dset, dataType, memspace, dspace, H5P_DEFAULT, data_out);
      if (size > 1) {
        printf("... first two values: %lld , %lld\n", data_out[0], data_out[1]);
        printf("... last two values: %lld, %lld\n", data_out[size - 2],
               data_out[size - 1]);
      } else {
        printf("... value=: %lld\n", data_out[0]);
      }
    }
  } else if (typeClass == H5T_COMPOUND) {
    if (s == 8) {

    } else if (s == 16) {
      struct H5VL_complex_double_t data_out[memSpaceSize];
      herr_t status =
          H5Dread(dset, dataType, memspace, dspace, H5P_DEFAULT, data_out);
      if (size > 1) {
        printf("... first two values: (%.3f+i%.3f) , (%.3f+i%.3f)\n",
               data_out[0].real, data_out[0].img, data_out[1].real,
               data_out[1].img);
        printf("... last two values:  (%.3f+i%.3f), (%.3f+i%.3f)\n",
               data_out[size - 2].real, data_out[size - 2].img,
               data_out[size - 1].real, data_out[size - 1].img);
      } else {
        printf("... value=: %.3f+i%.3f\n", data_out[0].real, data_out[0].img);
      }
    }
  } else {
    printf(" ......... print type??  %d  size=%ld\n", typeClass, s);
  }
}

void testHyperslabAsBB(int ndims, hsize_t *dims, hid_t dspace, hid_t dset,
                       hid_t dataType) {
  // read box
  hsize_t count[ndims], offset[ndims];
  int i = 0;
  for (i = 0; i < ndims; i++) {
    count[i] = dims[i] / 2;
    offset[i] = dims[i] / 3;
    if ((count[i] == 0) && dims[i] > 0) {
      count[i] = 1;
    }
    // printf(" \t bbox: %dth dim: s=%llu, c=%llu ", i, offset[i], count[i]);
  }

  hsize_t size = 1;
  for (i = 0; i < ndims; i++) {
    // hsize_t size=count[0] * count[1];
    size *= count[i];
  }

  if (size == 0) {
    printf("Nothing to retrieve.\n");
    return;
  }

  printBBox(count, offset, ndims);

  hid_t memspace = H5Screate_simple(ndims, count, NULL);
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  readData(dataType, memspace, dspace, size, dset);
  H5Sclose(memspace);

  // H5Tclose(nativeType);
}

void testSLABALL(int ndims, hsize_t *dims, hid_t dspace, hid_t dset,
                 hid_t dataType) {
  // create hyperslab
  hsize_t count[ndims], offset[ndims], stride[ndims], block[ndims];
  int i = 0;
  for (i = 0; i < ndims; i++) {
    if (dims[i] > 1) {
      offset[i] = 1;
    } else {
      offset[i] = 0;
    }
    block[i] = dims[i] / 3;
    if (block[i] == 0) {
      block[i] = 1;
    }

    if (dims[i] >= 4) {
      stride[i] = block[i] + 1;
      count[i] = 2;
    } else {
      stride[i] = 1;
      count[i] = 1;
    }
  }

  hsize_t size = 1;
  hsize_t mcount[ndims];
  for (i = 0; i < ndims; i++) {
    size *= count[i] * block[i];
    mcount[i] = count[i] * block[i];
  }

  hid_t memspace = H5Screate_simple(ndims, mcount, NULL);

  readData(dataType, memspace, dspace, size, dset);
  H5Sclose(memspace);

  // H5Tclose(nativeType);
}

void testSLABSLAB(int ndims, hsize_t *dims, hid_t dspace, hid_t dset,
                  hid_t dataType) {
  // create hyperslab
  hsize_t count[ndims], offset[ndims], stride[ndims], block[ndims];
  int i = 0;
  for (i = 0; i < ndims; i++) {
    if (dims[i] > 1) {
      offset[i] = 1;
    } else {
      offset[i] = 0;
    }
    block[i] = dims[i] / 3;
    if (block[i] == 0) {
      block[i] = 1;
    }

    if (dims[i] >= 4) {
      stride[i] = block[i] + 1;
      count[i] = 2;
    } else {
      stride[i] = 1;
      count[i] = 1;
    }
  }

  hid_t memspace = H5Screate_simple(ndims, dims, NULL);
  hsize_t moffset[ndims], mcount[ndims];

  hsize_t size = 1;
  for (i = 0; i < ndims; i++) {
    // hsize_t size=count[0] * count[1];
    size *= count[i] * block[i];
    moffset[i] = 0;
    mcount[i] = count[i] * block[i];
  }

  if (size == 0) {
    printf("Nothing to retrieve.\n");
    return;
  }

  printf("\n\t::Using hyperslab ::\n");
  for (i = 0; i < ndims; i++) {
    printf("\t%d th dim: total=%llu,\t ", i, dims[i]);
    printf("\t   offset:%llu, count %llu, stride=[%llu], block=[%llu]",
           offset[i], count[i], stride[i], block[i]);
    printf("\t      moffset %llu, mcount %llu \n", moffset[i], mcount[i]);
  }
  printf("\n");

#ifndef NEVER
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, moffset, NULL, mcount, NULL);
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, stride, count, block);

  readData(dataType, memspace, dspace, size, dset);
  H5Sclose(memspace);
#else
  printf(" do nothing \n");
#endif
  // H5Tclose(nativeType);
}

void testParallelism(int ndims, hsize_t *dims, hid_t dspace, hid_t dset,
                     hid_t dataType) {
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpisize;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpisize);

  // create hyperslab
  hsize_t count[ndims], offset[ndims], stride[ndims], block[ndims];
  int i = 0;
  for (i = 0; i < ndims; i++) {
    int step = dims[i] / mpisize;
    offset[i] = rank * step;

    block[i] = 1;
    stride[i] = 1;

    count[i] = 1;
    if (step > 2)
      count[i] = step - 2;
    else
      return;
  }

  printf("  CHECKing parallelism rank=%d\n", rank);
  if (rank == 0) {
    for (i = 0; i < ndims; i++) {
      printf(":: %llu ", dims[i]);
    }
    printf("\n");
  }
  printBBox(count, offset, ndims);

  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, stride, count, block);

  hsize_t size = 1;
  hsize_t mcount[ndims];
  for (i = 0; i < ndims; i++) {
    size *= count[i] * block[i];
    mcount[i] = count[i] * block[i];
  }

  hid_t memspace = H5Screate_simple(ndims, mcount, NULL);
  readData(dataType, memspace, dspace, size, dset);
  H5Sclose(memspace);
}

void testALLSLAB(int ndims, hsize_t *dims, hid_t dspace, hid_t dset,
                 hid_t dataType) {
  int i = 0;
  hsize_t largeDims[ndims];

  for (i = 0; i < ndims; i++) {
    largeDims[i] = 2 * dims[i];
  }

  hid_t memspace = H5Screate_simple(ndims, largeDims, NULL);
  hsize_t moffset[ndims], mcount[ndims];

  hsize_t size = 1;
  for (i = 0; i < ndims; i++) {
    // hsize_t size=count[0] * count[1];
    size *= dims[i];
    moffset[i] = dims[i] / 2;
    mcount[i] = dims[i];
  }

  if (size == 0) {
    printf("Nothing to retrieve.\n");
    return;
  }

  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, moffset, NULL, mcount, NULL);

  readData(dataType, memspace, H5S_ALL, size, dset);
  H5Sclose(memspace);
}

/*
void testHyperslab(hsize_t* count, hsize_t* offset, int ndim, hid_t dspace,
hid_t dset, hid_t dataType)
{
  hsize_t size = 1;
  int i=0;
  for (int i=0; i<ndim; i++) {
    //hsize_t size=count[0] * count[1];
    size *= count[i];
  }

  if (size == 0) {
    return;
  }

  printBBox(count, offset, ndim);

  hid_t memspace = H5Screate_simple(ndim, count, NULL);
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  readData(dataType, memspace, dspace, size, dset);
  H5Sclose(memspace);

  //H5Tclose(nativeType);
}
*/

void readAllFloats(hid_t dset, hid_t dataType, const char *varName) {
  size_t s = H5Tget_size(dataType);
  printf("\tvar %s  is -- scalar -- type is h5 float ,... size=%lu\n", varName,
         s);

  if (s == 4) {
    float value;
    H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    printf("value=%.3f\n", value);
  } else if (s == 8) {
    double value;
    H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    printf("value=%.3f\n", value);
  }
}

void readAllInts(hid_t dset, hid_t dataType, const char *varName) {
  size_t s = H5Tget_size(dataType);
  H5T_sign_t checkSign = H5Tget_sign(dataType);

  if (checkSign == H5T_SGN_NONE) {
    printf("\tvar %s  is -- scalar -- type is <unsigned> int ,... size=%lu\n",
           varName, s);
    if (s == 4) {
      uint value;
      H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      printf("value=%u\n", value);
    } else if (s == 1) {
      uint8_t value;
      H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      printf("value=%hhu\n", value);
    } else if (s == 2) {
      uint16_t value;
      H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      printf("value=%hu\n", value);
    } else if (s == 8) {
      int64_t value;
      H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      printf("value=%llu\n", value);
    }
  } else {
    printf("\tvar %s  is -- scalar -- type is int ,... size=%lu\n", varName, s);
    if (s == 4) {
      int value;
      H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      printf("value=%d\n", value);
    } else if (s == 1) {
      int8_t value;
      H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      printf("value=%hhd\n", value);
    } else if (s == 2) {
      int16_t value;
      H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      printf("value=%hd\n", value);
    } else if (s == 8) {
      int64_t value;
      H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
      printf("value=%lld\n", value);
    }
  }
}

void browseVarGroup(hid_t file_id, const char *varName) {

  hid_t varGroup_id = H5Gopen(file_id, varName, H5P_DEFAULT);
  if (varGroup_id < 0)
    return;

  //
  // attrs not yet implemented
  //
  // hid_t timestepAttrId = H5Aopen(varGroup_id, "timestep", H5P_DEFAULT);
  hid_t timestepAttrId = H5Aopen(varGroup_id, ADIOS_TIMESTEP, H5P_DEFAULT);
  unsigned int varTimeSteps;
  H5Aread(timestepAttrId, H5T_NATIVE_UINT, &varTimeSteps);
  H5Aclose(timestepAttrId);
  // status = H5Literate (varGroup_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
  // iter_callback_func, &results);
  int i = 0;
  int MAX_NAME_LENGTH = 100;
  for (i = 0; i < varTimeSteps; i++) {
    char tmpName[MAX_NAME_LENGTH];
    // depreciated H5Gget_objname_by_idx(varGroup_id, i, tmpName,
    // MAX_NAME_LENGTH);
    H5Lget_name_by_idx(varGroup_id, varName, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE,
                       i, tmpName, MAX_NAME_LENGTH, H5P_DEFAULT);
    hid_t dset = H5Dopen(varGroup_id, tmpName, H5P_DEFAULT);

    hid_t dataType = H5Dget_type(dset);

    if (dataType == 0) {
      printf("Unable to handle this type in h5\n");
      continue;
    }

    hid_t dspace = H5Dget_space(dset);

    // If your data space is simple (i.e. not null or scalar), then you can get
    // the number of dimensions using H5Sget_simple_extent_ndims:
    if (dspace != H5S_NULL) {
      if (!H5Sis_simple(dspace)) {
        return;
      }
      const int ndims = H5Sget_simple_extent_ndims(dspace);
      if (ndims == 0) {
        // printf("var %s  is -- scalar -- ... assuming .. integer ..\n",
        // varName);  uint32_t value;  H5Dread(dset, H5T_NATIVE_UINT, H5S_ALL,
        // H5S_ALL, H5P_DEFAULT, &value);

        H5T_class_t typeClass = H5Tget_class(dataType);
        if (typeClass == H5T_INTEGER) {
          readAllInts(dset, dataType, varName);
        } else if (typeClass == H5T_FLOAT) {
          size_t s = H5Tget_size(dataType);
          readAllFloats(dset, dataType, varName);
        } else if (typeClass == H5T_STRING) {
          size_t s = H5Tget_size(dataType);
          printf("\tvar %s  is -- scalar -- type is h5 string ,... size=%lu\n",
                 varName, s);

          char value[s + 1];
          H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
          printf("value=%s\n", value);

        } else if (typeClass == H5T_COMPOUND) {
          size_t s = H5Tget_size(dataType);
          printf("compound type: size=%lu\n", s);
          if (s == 8) {
            // float value[2];
            struct H5VL_complex_float_t value;
            H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
            printf("value=%.3f+i%.3f\n", value.real, value.img);
          } else if (s == 16) {
            // double value[2];
            struct H5VL_complex_double_t value;
            H5Dread(dset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
            printf("value=%.3f+i%.3f\n", value.real, value.img);
          }
        } else {
          printf(
              "\tscalar var: %s. Sorry the assocaited type is not handled \n",
              varName);
        }
      } else { // ndims > 0;
        int i = 0;
        printf("\tvar %s  is -- regular -- \n", varName);
        hsize_t dims[ndims];
        H5Sget_simple_extent_dims(dspace, dims, NULL);
#ifndef SERIAL
        // testHyperslab(count, offset, ndims, dspace, dset, dataType);
        // testHyperslabAsBB(ndims, dims, dspace, dset, dataType);

        printf("\n.... (1) test SLAB to SLAB .. \n");
        testSLABSLAB(ndims, dims, dspace, dset, dataType);

        printf("\n .... (2) test ALL to ALL .. \n");
        hsize_t size = 1;
        for (i = 0; i < ndims; i++) {
          size *= dims[i];
        }
        readData(dataType, H5S_ALL, H5S_ALL, size, dset);

        printf("\n .... (3) test ALL to SLAB .. \n");
        testALLSLAB(ndims, dims, dspace, dset, dataType);

// printf("\n .... (4) test SLAB to ALL .. \n");
// testSLABALL(ndims, dims, dspace, dset, dataType);

#else
        if (strlen(varName) >= 16) {
          testParallelism(ndims, dims, dspace, dset, dataType);
        }
#endif
      }
    }

    // close the dspace too?
    // H5Tclose(dataType); // dataType is from adios_vol.c which returns
    // immutable types, so dont delete
    H5Sclose(dspace);
    // H5Gclose(varGroup_id);
  }
  H5Gclose(varGroup_id);
}

void readMe(const char *file_name) {
  herr_t status;
  hid_t file_id = 0, acc_tpl = 0;
  ssize_t len;
  hid_t fapl_vol = H5Pcreate(H5P_FILE_ACCESS);

#ifdef SERIAL_ACCESS
  acc_tpl = H5Pcreate(H5P_FILE_ACCESS);

  if (H5Pset_fapl_adios(&acc_tpl, &fapl_vol) < 0)
    return;
#else
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
  assert(acc_tpl != FAIL);
  /* set Parallel access with communicator */
  int ret = H5Pset_fapl_mpio(acc_tpl, comm, info);
  if (H5Pset_fapl_adios(&acc_tpl, &fapl_vol) < 0)
    return;
  assert(ret != FAIL);
#endif

  //
  // file_id = H5Fcreate(file_name, H5F_ACC_RDONLY, H5P_DEFAULT, acc_tpl);
  file_id = H5Fopen(file_name, H5F_ACC_RDONLY, acc_tpl);
  if (file_id < 0)
    return;

  int namesize = strlen(file_name);
  char name[namesize + 1], fullpath[500];

  // just a test for H5S
  // hsize_t fdim[] = {10};
  // hid_t fid = H5Screate_simple(1, fdim, fdim);

  len = H5VLget_driver_name(file_id, name, 25);
  printf("FILE VOL name = %s \t ", name);

  len = H5Fget_name(file_id, name, 50);
  printf("Input file name=%s \t", name);

  hsize_t fsize;
  H5Fget_filesize(file_id, &fsize);
  printf("File size => %lld\t", fsize);

// testing whether a parameter exists
#ifdef NEVER
  herr_t check = H5Lexists(file_id, "T", H5P_DEFAULT);
  printf("var: T exists?  %d\n", check);
// expect to fail when checking T1:
// check = H5Lexists(file_id, "T1", H5P_DEFAULT); printf("var: T1 exists ?
// %d\n", check);
#endif

  // get # of vars in file
  H5G_info_t ginfo;
  status = H5Gget_info(file_id, &ginfo);
  ssize_t varcount = ginfo.nlinks;
  printf("Variable  count=%ld\n", varcount);

// now there are two wayas to browse the datasets
// A: ask H5G to return by idx
//    H5Gget_objname_by_idx(loc_id, idx, name, size) is depreciated
//    and replaced by a more complicated H5Lget_name_by_idx() which needs to
//    call H5Lget_name_by_idx() to get the name @ idx and size This is Not as
//    easy as H5Literate() So will not try here
//
// B: use H5Literate to get varNames[varcount]
//    will be explored below
//

// hid_t*  varid_list = (hid_t *)HDcalloc((size_t)varcount, sizeof(hid_t));
// ssize_t idcount =H5Fget_obj_ids(file_id, H5F_OBJ_DATASET, varcount,
// varid_list);  assert(varcount == idcount);  HDfree(oid_list);

//
// note: adios_vol returns a list of var names via h5literate
//
//#define ITERATE_AND_BROWSE_EACH_VAR

#ifdef ITERATE_AND_BROWSE_EACH_VAR
  char **varNames = SAFE_CALLOC(varcount, sizeof(char *));
  printf("iterating and browsing each variable .. \n");
  status = H5Literate(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
                      iter_callback_func, (void *)varNames);

  int i = 0;
  for (i = 0; i < varcount; i++) {
    browseVar(file_id, varNames[i]);
  }
  free(varNames);
#else
  printf("Using call backs to handle each var \n");

  H5VL_adios_varnames_t results;
  results.fid = file_id;
  results.varNames = SAFE_CALLOC(varcount, sizeof(char *));

  status = H5Literate(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
                      iter_callback_func, &results);

  int i = 0;
  for (i = 0; i < varcount; i++) {
    SAFE_FREE(results.varNames[i]);
  }
  SAFE_FREE(results.varNames);
#endif
#undef ITERATE_AND_BROWSE_EACH_VAR

  if (file_id)
    H5Fclose(file_id);
  if (acc_tpl)
    H5Pclose(acc_tpl);
  H5Pclose(fapl_vol);

  // H5Pclose(under_fapl);

  // H5VLunregister(H5VL_ADIOS_g);
  H5P_unset_adios();
  // MPI_Finalize();
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  if (argc == 1) {
    // default
    readMe("tiny.bp");
  } else {
    readMe(argv[1]);
  }

  MPI_Finalize();
  return 0;
}
