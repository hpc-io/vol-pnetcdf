//#define H5VL_FRIEND // junmin added

#include "H5Dprivate.h" // H5D_t
#include "H5Iprivate.h" // H5I_object_verify
#include "hdf5.h"
#define H5S_FRIEND  // suppress error for H5Spkg
#include "H5Spkg.h" // H5S_hyper_dim
#define H5O_FRIEND  // suppress error for H5Opkg
#include "H5Opkg.h"
//#include "mpi.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "H5Epublic.h"
#include "adios_vol.h"
#include "error.h"

static hid_t H5VL_ADIOS_g = -1;
static MPI_Comm H5VL_MPI_COMM = MPI_COMM_SELF;
static MPI_Info H5VL_MPI_INFO = MPI_INFO_NULL;
static int HDF5_MAX_STR_SIZE_IN_ARRAY = 30;

static int h5i_search_func(void *obj, hid_t id, void *key) {
  if (key == obj) {
    printf("id = %llu\n", id);
    return 1;
  } else
    return 0;
}

//
// either "ts_number" or ".../ts_number"
// is a valid dataset name
//
static int isDatasetName(const char *name) {
  int len = strlen(name);
  while (true) {
    char c = name[len - 1];
    if (c >= '0' && c <= '9') {
      len--;
      if (len == 0) {
        return 1;
      }
    } else if ((len < strlen(name)) && (c == '/')) {
      return 1;
    } else {
      break;
    }
  }
  return 0;
}

//
// returns numOfDigits of timestep
//
static int getDatasetTimeStep(const char *name, int *ts) {
  int power = -1;
  int len = strlen(name);
  while (true) {
    char c = name[len - 1];
    if (c >= '0' && c <= '9') {
      if (*ts < 0)
        *ts = 0; // initialize

      power++;
      int digit = c - '0';
      *ts += digit * pow(10, power);
      len--;
      if (len == 0) {
        return (power + 1);
      }
    } else if ((len < strlen(name)) && (c == '/')) {
      return (power + 1);
    } else {
      break;
    }
  }

  return power;
}

static void *assignDim(struct adios_index_comp_struct_v1 *varInAdios /* in */,
                       H5VL_adios_var_t *var /* out */) {
  int idx = adios_get_ch_idx(varInAdios, var->curr_timestep);

  // ADIOS_VOL_LOG_ERR("out of bound timestep: %d\n", val->curr_timestep);

  REQUIRE_SUCC_MSG((idx >= 0), NULL, "out of bound timestep: %d\n",
                   var->curr_timestep);

  int ndim = adios_get_var_ndim(var->fileReader, idx, varInAdios,
                                &(var->is_global_dim));
  if (ndim < 0) {
    ADIOS_VOL_LOG_ERR("Not a file or file object");
    return NULL;
  }

  var->dimInfo->ndim = ndim;

  // var->is_global_dim = 1;
  if (var->dimInfo->ndim > 0) {
    var->dimInfo->dims =
        (hsize_t *)SAFE_CALLOC(var->dimInfo->ndim, sizeof(hsize_t));

    int i;
    for (i = 0; i < var->dimInfo->ndim; i++) {
      // dims order: local/global/offset for each dimension
      uint64_t temp = varInAdios->characteristics[idx].dims.charDims[i * 3 + 1];
      if (temp > 0) {
        var->dimInfo->dims[i] = temp;
      } else {
        var->dimInfo->dims[i] =
            varInAdios->characteristics[idx].dims.charDims[i * 3];
      }
      // printf(" [%dth: %llu] ", i, var->dimInfo->dims[i]);
    }

    if (ADIOS_YES == var->fileReader->adios_host_language_fortran) {
      // file was in fortran order
      for (i = 0; i < ndim / 2; i++) {
        uint64_t temp = var->dimInfo->dims[i];
        var->dimInfo->dims[i] = var->dimInfo->dims[ndim - i - 1];
        var->dimInfo->dims[ndim - i - 1] = temp;
      }
    }
  }
  return var;
}

static void *createVar(H5_bp_file_t *fileReader, const char *pathWithTimeStep) {
  int ts = -1;
  int digits = getDatasetTimeStep(pathWithTimeStep, &ts);

  int pathLen = strlen(pathWithTimeStep);
  char varName[pathLen];
  memset(varName, '\0', sizeof(varName));
  if ('/' == pathWithTimeStep[0]) {
    strncpy(varName, pathWithTimeStep + 1, pathLen - 1 - digits - 1);
  } else {
    strncpy(varName, pathWithTimeStep, pathLen - digits - 1);
  }

  struct adios_index_comp_struct_v1 *varInAdios =
      adios_get_var_byname(fileReader, varName);

  REQUIRE_NOT_NULL_ERR(varInAdios, NULL);

  H5VL_adios_var_t *var =
      (H5VL_adios_var_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_var_t));

  var->name = (char *)SAFE_CALLOC(strlen(varName) + 1, sizeof(char));
  sprintf(var->name, "%s", varName);
  var->fileReader = fileReader;

  var->dimInfo = (H5VL_adios_dim_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_dim_t));

  var->ntimestep = adios_get_var_nsteps(varInAdios);

  var->curr_timestep = ts;
  if (ts >= var->ntimestep) {
    SAFE_FREE(var);
    ADIOS_VOL_LOG_ERR(" timestep given overflowed.");
    return NULL;
  }

  assignDim(varInAdios, var);

  return var;
}

//
// static void adios_vol_not_supported(const char* msg) {printf("%s is not yet
// supported in ADIOS VOL\n", msg);}

// static hid_t toHDF5type(enum ADIOS_DATATYPES adiosType) {
static hid_t toHDF5type(struct adios_index_comp_struct_v1 *adiosVar) {
  enum ADIOS_DATATYPES adiosType = adiosVar->type;
  //
  // this is here to make sure H5Tget_class(type) returns H5T_STRING
  // actual size here is not determined here
  //

  int size_of_type =
      adios_get_type_size(adiosType, adiosVar->characteristics[0].value);

  switch (adiosType) {
  case adios_string: {
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, size_of_type);
    return atype;
  }
  case adios_byte:
    return H5T_NATIVE_CHAR;

  case adios_unsigned_byte:
    return H5T_NATIVE_UCHAR;

  case adios_string_array: {
    int maxStrSize = HDF5_MAX_STR_SIZE_IN_ARRAY;
    hid_t atype = H5Tcopy(H5T_C_S1);
    herr_t ret = H5Tset_size(atype, maxStrSize);
    ret = H5Tset_strpad(atype, H5T_STR_NULLTERM);
    return atype;
    // return H5T_STRING;
  }
  case adios_short:
    return H5T_NATIVE_SHORT;

  case adios_unsigned_short:
    return H5T_NATIVE_USHORT;

  case adios_integer:
    return H5T_NATIVE_INT;

  case adios_unsigned_integer:
    return H5T_NATIVE_UINT;

  case adios_long:
    return H5T_NATIVE_LONG;

  case adios_unsigned_long:
    return H5T_NATIVE_ULONG;

  case adios_real:
    return H5T_NATIVE_FLOAT;

  case adios_double:
    return H5T_NATIVE_DOUBLE;

  case adios_long_double:
    return H5T_NATIVE_LDOUBLE;

  case adios_complex: {
    size_t h5float_size = H5Tget_size(H5T_NATIVE_FLOAT);
    hid_t m_DefH5TypeComplexFloat = H5Tcreate(H5T_COMPOUND, 2 * h5float_size);

    H5Tinsert(m_DefH5TypeComplexFloat, "freal", 0, H5T_NATIVE_FLOAT);
    H5Tinsert(m_DefH5TypeComplexFloat, "fimg", h5float_size, H5T_NATIVE_FLOAT);

    return m_DefH5TypeComplexFloat;
  }
  case adios_double_complex: {
    size_t h5double_size = H5Tget_size(H5T_NATIVE_DOUBLE);
    hid_t m_DefH5TypeComplexDouble = H5Tcreate(H5T_COMPOUND, 2 * h5double_size);

    H5Tinsert(m_DefH5TypeComplexDouble, "dreal", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert(m_DefH5TypeComplexDouble, "dimg", h5double_size,
              H5T_NATIVE_DOUBLE);

    return m_DefH5TypeComplexDouble;
  }
  // return mH5T_COMPLEX_DOUBLE;
  default:
    SHOW_ERROR_MSG("unknown adios type %d found.", adiosType);
    return 0;
  }
}

static herr_t H5VL_adios_init(hid_t vipl_id) {
  // printf("------- ADIOS INIT, called at H5VLinitialize\n");
  return 0;
}

static herr_t H5VL_adios_term(hid_t vtpl_id) {
  // printf("------- ADIOS TERMINATE called at H5VLunregister\n");
  return 0;
}

int call_adios_get_num_vars(H5_bp_file_t *f) {
  int num = 0;

  struct adios_index_comp_struct_v1 *v = f->vars_root;
  while (v) {
    num++;
    v = v->next;
  }
  return num;
}

int call_adios_get_num_attributes(H5_bp_file_t *f) {
  int num = 0;

  struct adios_index_comp_struct_v1 *v = f->attrs_root;
  while (v) {
    num++;
    v = v->next;
  }
  return num;
}

int call_is_valid_bpfile(const char *name, H5_bp_file_t *bpFile) {
  REQUIRE_NOT_NULL(bpFile);
  // printf("  check adios validity here.\n");

  int isBpFile = 1;
  int err, flag = 0;

  MPI_File fh = 0;

#ifdef NEVER
  err = MPI_File_open(MPI_COMM_SELF, (char *)name, MPI_MODE_RDONLY,
                      (MPI_Info)MPI_INFO_NULL, &fh);
#else
  err = MPI_File_open(H5VL_MPI_COMM, (char *)name, MPI_MODE_RDONLY,
                      H5VL_MPI_INFO, &fh);

#endif
  // H5Adios_list_MPI_Error(err);
  REQUIRE_MPI_SUCC(err);

  REQUIRE_SUCC((0 == call_init(bpFile, fh)), -1);

  REQUIRE_SUCC((0 == call_seek_back(bpFile, ADIOS_BP_MINIFOOTER_SIZE)), -1);

  int offset = ADIOS_BP_MINIFOOTER_SIZE - 4;

  REQUIRE_SUCC((0 == call_adios_get_version(bpFile, offset)), -1);

  // printf("version=%d\n", bpFile->bp_version);

  if (bpFile->bp_version < 3) {
    printf("bp version is not up to date. \n");
    return 0;
  }

  bpFile->file_name = SAFE_MALLOC(strlen(name) + 1);

  REQUIRE_NOT_NULL(bpFile->file_name);

  sprintf(bpFile->file_name, "%s", name);
  bpFile->file_name[strlen(name)] = '\0';

  return isBpFile;
}

static void cleanUp(H5VL_adios_t *file) {
  call_close(file->bpFileReader);
  SAFE_FREE(file->bpFileReader);
  // free(file->bpFileReader);
  free(file);
}

static void *H5VL_adios_file_open(const char *name, unsigned flags,
                                  hid_t fapl_id, hid_t dxpl_id, void **req) {
  H5VL_adios_t *file;

  file = (H5VL_adios_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_t));
  // file->bpFileReader = (H5_bp_file_t*) SAFE_MALLOC(sizeof (H5_bp_file_t));
  file->bpFileReader = (H5_bp_file_t *)SAFE_CALLOC(1, sizeof(H5_bp_file_t));

  REQUIRE_SUCC_ACTION((1 == call_is_valid_bpfile(name, file->bpFileReader)),
                      cleanUp(file), NULL);

  // test read metadata
  REQUIRE_SUCC_ACTION(
      (0 == call_seek(file->bpFileReader, file->bpFileReader->pg_index)),
      cleanUp(file), NULL);

  return (void *)file;
}

static herr_t H5VL_adios_file_get(void *file, H5VL_file_get_t get_type,
                                  hid_t dxpl_id, void **req,
                                  va_list arguments) {
  herr_t ret_value = SUCCEED;

  H5VL_adios_t *f = (H5VL_adios_t *)file;

  switch (get_type) {
  /* H5Fget_access_plist */
  case H5VL_FILE_GET_FAPL:
    ADIOS_VOL_NOT_SUPPORTED_ERR("H5VL_FILE_GET_FAPL\n");
    break;
  case H5VL_FILE_GET_FCPL:
    ADIOS_VOL_NOT_SUPPORTED_ERR("H5VL_FILE_GET_FCPL\n");
    break;
  case H5VL_FILE_GET_OBJ_COUNT: {
    int type = va_arg(arguments, unsigned int);
    ssize_t *ret = va_arg(arguments, ssize_t *);

    if (type == H5F_OBJ_DATASET) {
      ssize_t result = call_adios_get_num_vars(f->bpFileReader);
      *ret = (ssize_t)result;
    } else if (type == H5F_OBJ_ATTR) {
      ssize_t result = call_adios_get_num_attributes(f->bpFileReader);
      *ret = (ssize_t)result;
    } else {
      ADIOS_VOL_NOT_SUPPORTED_ERR("count types other than dataset\n");
    }
    break;
  }
  case H5VL_FILE_GET_OBJ_IDS: {
    unsigned type = va_arg(arguments, unsigned int);
    size_t max_objs = va_arg(arguments, size_t);
    hid_t *oid_list = va_arg(arguments, hid_t *);
    ssize_t *ret = va_arg(arguments, ssize_t *);
    size_t obj_count = 0; /* Number of opened objects */

    if (type == H5F_OBJ_DATASET) {
      printf("todo: return var ids \n");
      // temp set to obj_count
      int i;
      for (i = 0; i < obj_count; i++) {
        oid_list[i] = i;
      }
      *ret = (ssize_t)max_objs;
    } else {
      ADIOS_VOL_NOT_SUPPORTED_ERR("get ids of types  other than dataset\n");
    }
    break;
  }
  case H5VL_FILE_GET_INTENT:
    ADIOS_VOL_NOT_SUPPORTED_ERR("H5VL_FILE_GET_INTENT\n");
    break;
  case H5VL_FILE_GET_NAME: {
    H5I_type_t type = va_arg(arguments, H5I_type_t);
    size_t size = va_arg(arguments, size_t);
    char *name = va_arg(arguments, char *);
    ssize_t *ret = va_arg(arguments, ssize_t *);

    size_t len = strlen(f->bpFileReader->file_name);
    if (name) {
      sprintf(name, "%s", f->bpFileReader->file_name);
      name[len] = '\0';
    }

    // Set the return value for the API call
    *ret = (ssize_t)len;

    break;
  }
  case H5VL_OBJECT_GET_FILE: {
    break;
  }
  default:
    ADIOS_VOL_NOT_SUPPORTED_ERR("H5VL_FILE_GET flag");
  }

  return ret_value;
}

static herr_t H5VL_adios_file_optional(void *file, hid_t dxpl_id,
                                       void /*H5_ATTR_UNUSED*/ **req,
                                       va_list arguments) {
  herr_t ret_value = SUCCEED; /* Return value */
  H5VL_adios_t *f = (H5VL_adios_t *)file;
  if (NULL == f) {
    return FAIL;
  }

  H5VL_file_optional_t optional_type = va_arg(arguments, H5VL_file_optional_t);

  switch (optional_type) {
  /* H5Fget_filesize */
  case H5VL_FILE_GET_SIZE: {
    hsize_t *ret = va_arg(arguments, hsize_t *);
    *ret = f->bpFileReader->file_size;
    break;
  }
  default:
    ret_value = FAIL;
    // ADIOS_VOL_NOT_SUPPORTED_ERR("Only file option supported is
    // get_file_size\n");
    ADIOS_VOL_NOT_SUPPORTED_ERR(
        "Only file option supported is get_file_size\n");
  }
  return ret_value;
}

static herr_t H5VL_adios_file_close(void *file, hid_t dxpl_id, void **req) {
  H5VL_adios_t *f = (H5VL_adios_t *)file;

  if (f->bpFileReader != NULL) {
    // call_seek(f->bpFileReader, f->bpFileReader->pg_index);
    call_close(f->bpFileReader);
    MPI_File_close(&(f->bpFileReader->mpi_handle));
    SAFE_FREE(f->bpFileReader->file_name);
    SAFE_FREE(f->bpFileReader);
  }
  free(f);

  // printf("\n------- ADIOS Close a file at H5FClose\n");
  return 1;
}

static void *H5VL_adios_group_open(void *file, H5VL_loc_params_t loc_params,
                                   const char *varName, hid_t gapl_id,
                                   hid_t dxpl_id, void H5_ATTR_UNUSED **req) {
  H5_bp_file_t *bpFile = NULL;
  if (H5I_GROUP == loc_params.obj_type) {
    H5VL_adios_var_t *var = (H5VL_adios_var_t *)file;
    bpFile = var->fileReader;
  } else if (H5I_FILE == loc_params.obj_type) {
    bpFile = ((H5VL_adios_t *)file)->bpFileReader;
  }

  if ((varName != NULL) && (1 == strlen(varName)) &&
      ('/' == varName[0])) { // h5ls request
    H5VL_adios_var_t *dummy =
        (H5VL_adios_var_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_var_t));
    dummy->name = NULL;
    dummy->fileReader = bpFile;
    return dummy;
  }

  H5VL_adios_var_t *var =
      (H5VL_adios_var_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_var_t));

  var->name = (char *)SAFE_CALLOC(strlen(varName) + 1, sizeof(char));
  sprintf(var->name, "%s", varName);
  var->fileReader = bpFile;

  var->dimInfo = (H5VL_adios_dim_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_dim_t));

  struct adios_index_comp_struct_v1 *varInAdios =
      adios_get_var_byname(var->fileReader, varName);

  var->ntimestep = adios_get_var_nsteps(varInAdios);

  return var;
  ;
}

static herr_t H5VL_adios_group_close(void *grp, hid_t H5_ATTR_UNUSED dxpl_id,
                                     void H5_ATTR_UNUSED **req) {
  herr_t ret_value = SUCCEED; /* Return value */

  H5VL_adios_var_t *v = (H5VL_adios_var_t *)grp;
  SAFE_FREE(v->name);
  if (v->dimInfo)
    SAFE_FREE(v->dimInfo->dims);
  SAFE_FREE(v->dimInfo);
  free(v);
  return ret_value;
}

static herr_t H5VL_adios_group_get(void *obj, H5VL_group_get_t get_type,
                                   hid_t dxpl_id, void H5_ATTR_UNUSED **req,
                                   va_list arguments) {
  herr_t ret_value = SUCCEED; /* Return value */

  switch (get_type) {
  case H5VL_GROUP_GET_INFO: {
    H5VL_adios_t *f = (H5VL_adios_t *)obj;
    H5G_info_t *grp_info = va_arg(arguments, H5G_info_t *);
    ssize_t result;
    if (f->bpFileReader != NULL) {
      result = call_adios_get_num_vars(f->bpFileReader);
    } else {
      H5VL_adios_var_t *var = (H5VL_adios_var_t *)obj;
      result = call_adios_get_num_vars(var->fileReader);
    }

    // H5VL_loc_params_t loc_params = va_arg (arguments, H5VL_loc_params_t);
    // H5G_info_t  *grp_info = va_arg (arguments, H5G_info_t *);

    H5VL_loc_params_t loc_params = va_arg(arguments, H5VL_loc_params_t);
    if (loc_params.type == H5VL_OBJECT_BY_SELF) { /* H5Gget_info */
      grp_info->nlinks = result;
    }

    break;
  }
  case H5VL_GROUP_GET_GCPL: {
    hid_t *new_gcpl_id = va_arg(arguments, hid_t *);
    *new_gcpl_id = H5Pcreate(H5P_GROUP_CREATE);

    break;
  }
  default:
    ADIOS_VOL_NOT_SUPPORTED_ERR("adios_group_get() only supports get info");
    ret_value = FAIL;
  }
  return ret_value;
}

static herr_t H5VL_adios_link_get(void *obj, H5VL_loc_params_t loc_params,
                                  H5VL_link_get_t get_type, hid_t dxpl_id,
                                  void H5_ATTR_UNUSED **req,
                                  va_list arguments) {
  herr_t ret_value = FAIL; /* Return value */
  switch (get_type) {
  /* H5Lget_name_by_idx */
  case H5VL_LINK_GET_NAME: {
    char *name = va_arg(arguments, char *);
    size_t size = va_arg(arguments, size_t);
    ssize_t *ret = va_arg(arguments, ssize_t *);

    hsize_t idx = loc_params.loc_data.loc_by_idx.n;
    sprintf(name, "%llu", idx);
    *ret = strlen(name);
  }
    ret_value = SUCCEED; /* Return value */
    break;
  case H5VL_LINK_GET_INFO: {
    H5L_info_t *linfo = va_arg(arguments, H5L_info_t *);
#ifdef NEVER
    /* Get the link information */
    H5G_loc_t loc;
    if (loc_params.type == H5VL_OBJECT_BY_NAME) { /* H5Lget_info */
      if (H5L_get_info(&loc, loc_params.loc_data.loc_by_name.name, linfo) < 0)
        ADIOS_VOL_LOG_ERR("unable to find link info by name");
      break;
    } else if (loc_params.type == H5VL_OBJECT_BY_IDX) { /* H5Lget_info_by_idx */
      if (H5L_get_info_by_idx(&loc, loc_params.loc_data.loc_by_idx.name,
                              loc_params.loc_data.loc_by_idx.idx_type,
                              loc_params.loc_data.loc_by_idx.order,
                              loc_params.loc_data.loc_by_idx.n, linfo) < 0)
        ADIOS_VOL_LOG_ERR("unable to find link info by idx");
      break;
    } else
      ADIOS_VOL_LOG_ERR("unable to get link info");
    break;
#endif
  }

  default:
    ADIOS_VOL_NOT_SUPPORTED_ERR("link_get()");
  }
  return ret_value;
}

//
// expects a file in the first input
//
static herr_t H5VL_adios_link_specific(void *file, H5VL_loc_params_t loc_params,
                                       H5VL_link_specific_t specific_type,
                                       hid_t dxpl_id, void H5_ATTR_UNUSED **req,
                                       va_list arguments) {
  herr_t ret_value = SUCCEED; /* Return value */

  /*
  H5VL_adios_t *f = (H5VL_adios_t *)file;


  H5_bp_file_t *bpFile = f->bpFileReader;

  bpFile... could be a corrupted pointer
  if (bpFile == NULL) {
    H5VL_adios_var_t* var = (H5VL_adios_var_t*)file;
    bpFile = var->fileReader;
  }
  */
  H5_bp_file_t *bpFile = NULL;
  if (H5I_GROUP == loc_params.obj_type) {
    H5VL_adios_var_t *var = (H5VL_adios_var_t *)file;
    bpFile = var->fileReader;
  } else if (H5I_FILE == loc_params.obj_type) {
    bpFile = ((H5VL_adios_t *)file)->bpFileReader;
  }

  REQUIRE_NOT_NULL(bpFile);

  switch (specific_type) {
  case H5VL_LINK_EXISTS: {
    if (call_adios_has_var(bpFile, loc_params.loc_data.loc_by_name.name) == 0) {
      ret_value = FAIL;
    }
    return ret_value;
  }
  case H5VL_LINK_ITER: {
    hbool_t recursive = va_arg(arguments, int);
    H5_index_t idx_type = va_arg(arguments, H5_index_t);
    H5_iter_order_t order = va_arg(arguments, H5_iter_order_t);
    hsize_t *idx_p = va_arg(arguments, hsize_t *);
    H5L_iterate_t callback = va_arg(arguments, H5L_iterate_t);
    void *op_data = va_arg(arguments, void *);

    if (callback) {
      if (H5I_FILE == loc_params.obj_type) {
        struct adios_index_comp_struct_v1 *var_ptr = bpFile->vars_root;
        H5L_info_t info;
        info.type = H5L_TYPE_HARD; // for h5ls
        if (!recursive) {
          while (var_ptr != NULL) {
            callback(GLOBAL_FILE_ID, var_ptr->name, &info, op_data);
            var_ptr = var_ptr->next;
          }
          // callback(0, "gndx", NULL, op_data);
          // callback(0, "T", NULL, op_data);
        } else {
          while (var_ptr != NULL) {
            callback(GLOBAL_FILE_ID, var_ptr->name, &info, op_data);
            int numTimeStep = adios_get_var_nsteps(var_ptr);
            char ts[100];
            int i;
            for (i = 0; i < numTimeStep; i++) {
              // sprintf(ts, "%s/TimeStep%d", var_ptr->name, i);
              sprintf(ts, "%s/%d", var_ptr->name, i);
              callback(GLOBAL_FILE_ID, ts, &info, op_data);
            }

            var_ptr = var_ptr->next;
          }
        }
      } else if (H5I_GROUP == loc_params.obj_type) {
#ifndef NEVER
        // not sure how to get gid from obj!!!! -1 is a temperory way of getting
        // the gid of the obj. does not alway work!! TRUE FALSE neither

        H5VL_adios_var_t *var = (H5VL_adios_var_t *)file;
        struct adios_index_comp_struct_v1 *var_ptr = bpFile->vars_root;

        if (op_data == NULL) {
          ADIOS_VOL_WARN("no gid passed to link iter callback ");
          ret_value = FAIL;
        } else {
          hid_t gid = *((hid_t *)op_data);
          if (NULL != var->name) {
            H5L_info_t info;
            info.type = H5L_TYPE_HARD; // for h5ls

            // callback(gid, "0", &info, op_data);
            char timeStepStr[10];
            int ts;
            for (ts = 0; ts < var->ntimestep; ts++) {
              sprintf(timeStepStr, "%d", ts);
              callback(gid, timeStepStr, &info, op_data);
            }

          } else {
            if (!recursive) {
              H5L_info_t info;
              info.type = H5L_TYPE_HARD; // for h5ls
              while (var_ptr != NULL) {
                callback(gid, var_ptr->name, &info, op_data);
                var_ptr = var_ptr->next;
              }
            } else {
              ADIOS_VOL_NOT_SUPPORTED("recursively list a subgroup");
            }
          }
        }

#else
        ADIOS_VOL_NOT_SUPPORTED("not yet list a subgroup");
#endif
      } else {
        // should not reach here
      }
    }

// printf("will browse link  ... just  manually adding var names here for
// now\n");

#ifdef ITERATE_AND_BROWSE_EACH_VAR
    char **data = (char **)op_data;
    data[0] = SAFE_CALLOC(1, sizeof(char));
    data[1] = SAFE_CALLOC(4, sizeof(char));
    data[2] = SAFE_CALLOC(4, sizeof(char));

    sprintf((char *)(data[0]), "T");
    sprintf((char *)(data[1]), "gndx");
    sprintf((char *)(data[2]), "gndx");
#else
// H5VL_adios_varnames_t* data = (H5VL_adios_varnames_t*)op_data;
#endif
    /*
     */
    break;
  }
  default: // e.g. case H5VL_LINK_DELETE:
  {
    ADIOS_VOL_NOT_SUPPORTED_ERR("link type other than exists and iterate ");
  }
  }
  return 1; // succ
}

static void *H5VL_adios_datatype_commit(void *obj, H5VL_loc_params_t loc_params,
                                        const char *name, hid_t type_id,
                                        hid_t lcpl_id, hid_t tcpl_id,
                                        hid_t tapl_id, hid_t dxpl_id,
                                        void **req) {
  H5VL_adios_t *dt;
  H5VL_adios_t *o = (H5VL_adios_t *)obj;

  dt = (H5VL_adios_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_t));

  printf("------- ADIOS H5Tcommit\n");
  return dt;
}
static void *H5VL_adios_datatype_open(void *obj, H5VL_loc_params_t loc_params,
                                      const char *name, hid_t tapl_id,
                                      hid_t dxpl_id, void **req) {
  H5VL_adios_t *dt;
  H5VL_adios_t *o = (H5VL_adios_t *)obj;

  dt = (H5VL_adios_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_t));

  // dt->under_plugin = o->under_plugin;

  printf("------- ADIOS H5Topen\n");
  return (void *)dt;
}

static herr_t H5VL_adios_datatype_get(void *obj, H5VL_datatype_get_t get_type,
                                      hid_t dxpl_id, void **req,
                                      va_list arguments)
// H5VL_adios_datatype_get(void *obj, unsigned char *buf, size_t size, void**
// req)
{

  H5T_t *dt = (H5T_t *)obj;
  herr_t ret_value = SUCCEED; /* Return value */

  switch (get_type) {
  case H5VL_DATATYPE_GET_BINARY: {
    ssize_t *nalloc = va_arg(arguments, ssize_t *);
    void *buf = va_arg(arguments, void *);
    size_t size = va_arg(arguments, size_t);

    if (H5T_encode(dt, (unsigned char *)buf, &size) < 0)
      printf(" oh no \n");

    *nalloc = (ssize_t)size;
    break;
  case H5VL_DATATYPE_GET_TCPL:
  default:
    ADIOS_VOL_NOT_SUPPORTED("adios vol datatype_get only handles binary");
  }
  }

  return ret_value;
}

static herr_t H5VL_adios_datatype_close(void *dt, hid_t dxpl_id, void **req) {
  H5VL_adios_t *type = (H5VL_adios_t *)dt;

  // assert(type->bpFileReader);
  // assert(type->under_plugin);

  free(type);

  printf("------- ADIOS H5Tclose\n");
  return 1;
}

static herr_t H5VL_adios_object_optional(void *obj, hid_t dxpl_id,
                                         void H5_ATTR_UNUSED **req,
                                         va_list arguments) {
  H5VL_object_optional_t optional_type =
      va_arg(arguments, H5VL_object_optional_t);
  H5VL_loc_params_t loc_params = va_arg(arguments, H5VL_loc_params_t);
  H5G_loc_t loc;              /* Location of group */
  herr_t ret_value = SUCCEED; /* Return value */

  switch (optional_type) {
  case H5VL_OBJECT_GET_INFO: {
    H5O_info_t *obj_info = va_arg(arguments, H5O_info_t *);
    if (loc_params.type == H5VL_OBJECT_BY_NAME) {
      const char *name = loc_params.loc_data.loc_by_name.name;
      obj_info->rc = 1; // needed for h5ls to know it is not a link
      if ((name != NULL) && (1 == strlen(name)) && ('/' == name[0])) {
        obj_info->type = H5O_TYPE_GROUP;
        break;
      }
      int val = atoi(name);

      // if (((val == 0) && ('0' == name[0])) || (val > 0)) {
      if (isDatasetName(name) > 0) {
        obj_info->type = H5O_TYPE_DATASET;
      } else {
        obj_info->type = H5O_TYPE_GROUP;
      }
    } else if (loc_params.type == H5VL_OBJECT_BY_SELF) {
      obj_info->rc = 1;
    } else {
      ADIOS_VOL_NOT_SUPPORTED_ERR("COMING SOON");
    }
    break;
  }
  case H5VL_OBJECT_GET_COMMENT: {
    char *comment = va_arg(arguments, char *);
    size_t bufsize = va_arg(arguments, size_t);
    ssize_t *ret = va_arg(arguments, ssize_t *);

    *ret = 0;
    comment = SAFE_MALLOC(1 * sizeof(char *));
    comment[0] = '\0';
    break;
  }
  default:
    ADIOS_VOL_NOT_SUPPORTED_ERR("H5VL_FILE_GET flag");
    return FAIL;
  }
  return ret_value;
}

static void *H5VL_adios_object_open(void *obj, H5VL_loc_params_t loc_params,
                                    H5I_type_t *opened_type, hid_t dxpl_id,
                                    void **req) {
  *opened_type = H5I_GROUP; // h5ls
  if (H5VL_OBJECT_BY_NAME == loc_params.type) {
    const char *name = loc_params.loc_data.loc_by_name.name;

    if (isDatasetName(name) > 0) {
      *opened_type = H5I_DATASET;
      H5VL_adios_t *f = (H5VL_adios_t *)obj;

      return createVar(f->bpFileReader, name);
    }
  }

  H5VL_adios_var_t *dummy =
      (H5VL_adios_var_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_var_t));
  return (void *)dummy;
}
/*
static herr_t H5VL_adios_object_specific(void *obj,
                                         H5VL_loc_params_t loc_params,
                                         H5VL_object_specific_t specific_type,
                                         hid_t dxpl_id, void **req,
                                         va_list arguments) {
  H5VL_adios_t *o = (H5VL_adios_t *)obj;

  printf("------- ADIOS Object specific\n");
  return 1;
}
*/
/*
static void *H5VL_adios_dataset_create(void *obj, H5VL_loc_params_t loc_params,
                                       const char *name, hid_t dcpl_id,
                                       hid_t dapl_id, hid_t dxpl_id,
                                       void **req) {
  H5VL_adios_t *dset;
  H5VL_adios_t *o = (H5VL_adios_t *)obj;

  dset = (H5VL_adios_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_t));

  printf("------- ADIOS H5Dcreate\n");
  return (void *)dset;
}

*/

static herr_t H5VL_adios_dataset_get(void *obj, H5VL_dataset_get_t get_type,
                                     hid_t dxpl_id, void **req,
                                     va_list arguments) {
  H5VL_adios_var_t *var = (H5VL_adios_var_t *)obj;
  herr_t ret_value = SUCCEED; /* Return value */

  hid_t *ret_id = va_arg(arguments, hid_t *);
  switch (get_type) {
  /* H5Dget_space */
  case H5VL_DATASET_GET_SPACE: {
    if (var->dimInfo->ndim == 0) {
      *ret_id = H5Screate(H5S_SCALAR);
    } else
      // if (var->is_global_dim >= 0)
      // not sure how the local variables will behave here
      // right now retrieve the datablock specified todimInfo in createVar()
      *ret_id = H5Screate_simple(var->dimInfo->ndim, var->dimInfo->dims, NULL);

    break;
  }
  case H5VL_DATASET_GET_TYPE: {
    struct adios_index_comp_struct_v1 *varInAdios =
        adios_get_var_byname(var->fileReader, var->name);

    *ret_id = toHDF5type(varInAdios);
    break;
  }
  case H5VL_DATASET_GET_DCPL: {
    *ret_id = H5Pcreate(H5P_DATASET_CREATE);
    break;
  }
  default:
    ADIOS_VOL_NOT_SUPPORTED_ERR(
        "only allows get type or space from dataset, for now");
    return FAIL;
  }
  return ret_value;
}

static void *H5VL_adios_attr_open(void *varGrp, H5VL_loc_params_t loc_params,
                                  const char *attr_name,
                                  hid_t H5_ATTR_UNUSED aapl_id, hid_t dxpl_id,
                                  void H5_ATTR_UNUSED **req)

{

  H5VL_adios_var_t *var = (H5VL_adios_var_t *)varGrp;

  H5VL_adios_attr_t *attr =
      (H5VL_adios_attr_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_attr_t));

  attr->name = (char *)SAFE_CALLOC(strlen(attr_name) + 1, sizeof(char));
  sprintf(attr->name, "%s", attr_name);
  attr->fileReader = var->fileReader;

  // attr->dimInfo = NULL;
  if (strcmp(attr_name, ADIOS_TIMESTEP) == 0) {
    attr->ts = var->ntimestep;
  }
  return attr;
}

// attr obj is actually an adios, as this is for
// getting var timestep
static herr_t H5VL_adios_attr_read(void *varGrp, hid_t dtype_id, void *buf,
                                   hid_t dxpl_id, void H5_ATTR_UNUSED **req) {
  herr_t ret_value = FAIL;

  H5VL_adios_attr_t *attr = (H5VL_adios_attr_t *)varGrp;

  if ((attr->name != NULL) && (strcmp(attr->name, ADIOS_TIMESTEP) == 0)) {
    unsigned int *result = (unsigned int *)buf;
    *result = attr->ts;
    ret_value = SUCCEED;
  } else {
    struct adios_index_comp_struct_v1 *attrInAdios =
        adios_get_attr_byname(attr->fileReader, attr->name);
    if (NULL != attrInAdios) {
      uint64_t typesize = adios_get_type_size(
          attrInAdios->type, attrInAdios->characteristics[0].value);
      const void *val = attrInAdios->characteristics[0].value;
      if (0 == attrInAdios->characteristics->dims.count) {
        memcpy(buf, val, typesize);
      } else {
        uint64_t attrArraySize = attrInAdios->characteristics->dims.charDims[0];

        if (adios_string_array == attrInAdios->type) {
          int N = HDF5_MAX_STR_SIZE_IN_ARRAY;
          int i;
          for (i = 0; i < attrArraySize; i++) {
            char *curr = ((char **)val)[i];
            memcpy(buf + N * i, curr, N);
          }
          // cann't do memcpy, val is char**, buf is expected to be char[N *
          // attrArraySize]
          // N = common size  of strings in char**
          // memcpy(buf, val, typesize*attrArraySize);
        } else {
          memcpy(buf, val, typesize * attrArraySize);
        }
      }
    } else {
      ADIOS_VOL_WARN("NO such attr:");
    }
  }

  return ret_value;
}

static herr_t H5VL_adios_attr_close(void *obj, hid_t H5_ATTR_UNUSED dxpl_id,
                                    void H5_ATTR_UNUSED **req) {
  herr_t ret_value = SUCCEED;

  H5VL_adios_attr_t *attr = (H5VL_adios_attr_t *)obj;
  SAFE_FREE(attr->name);
  // if (attr->dimInfo)
  // SAFE_FREE(attr->dimInfo->dims);
  // SAFE_FREE(attr->dimInfo);
  free(attr);

  return ret_value;
}

static herr_t H5VL_adios_attr_get(void *obj, H5VL_attr_get_t get_type,
                                  hid_t dxpl_id, void H5_ATTR_UNUSED **req,
                                  va_list arguments) {
  herr_t ret_value = SUCCEED; /* Return value */
  hid_t *ret_id = va_arg(arguments, hid_t *);

  H5VL_adios_attr_t *attr = (H5VL_adios_attr_t *)obj;
  struct adios_index_comp_struct_v1 *attrInAdios =
      adios_get_attr_byname(attr->fileReader, attr->name);

  switch (get_type) {
  case H5VL_ATTR_GET_ACPL: {
    *ret_id = H5Pcreate(H5P_ATTRIBUTE_CREATE);
    break;
  }
  case H5VL_ATTR_GET_TYPE: {
    //*ret_id = H5T_NATIVE_UINT;
    if (attrInAdios != NULL)
      *ret_id = toHDF5type(attrInAdios);
    else // the ADIOS_TIMESTEP
      *ret_id = H5T_NATIVE_UINT;
    break;
  }
  case H5VL_ATTR_GET_SPACE: {
    if (NULL == attrInAdios) // the ADIOS_TIMESTEP
      *ret_id = H5Screate(H5S_SCALAR);
    else {
      int attrDim = adios_get_attr_ndim(attrInAdios);
      if (0 == attrDim)
        *ret_id = H5Screate(H5S_SCALAR);
      else {
        hsize_t dims[1];
        dims[0] = attrDim;
        *ret_id = H5Screate_simple(1, dims, NULL);
      }
    }
    break;
  }
  default:
    ADIOS_VOL_WARN("Not supported flag of vol_attr_get");
    ret_value = FAIL;
  }
  return ret_value;
}
static herr_t H5VL_adios_attr_specific(void *obj, H5VL_loc_params_t loc_params,
                                       H5VL_attr_specific_t specific_type,
                                       hid_t dxpl_id, void H5_ATTR_UNUSED **req,
                                       va_list arguments) {
  herr_t ret_value = SUCCEED; /* Return value */

  /*
  if (H5I_FILE == loc_params.obj_type) {
    ADIOS_VOL_WARN("Not yet iteratiing file attributes ");
    return FAIL;
  }
  */
  if (H5I_DATASET == loc_params.obj_type) {
    return ret_value;
  }

  if (H5I_GROUP != loc_params.obj_type) {
    ADIOS_VOL_WARN(
        "does not supoort attribute iteration unless it is at root level");
    return FAIL;
  }

  H5VL_adios_var_t *var = (H5VL_adios_var_t *)obj;

  if (var == NULL) {
    ADIOS_VOL_WARN("null obj. No attrs");
    return FAIL;
  }
  /*
   */
  switch (specific_type) {
  case H5VL_ATTR_ITER: {
    H5_index_t idx_type = va_arg(arguments, H5_index_t);
    H5_iter_order_t order = va_arg(arguments, H5_iter_order_t);
    hsize_t *idx = va_arg(arguments, hsize_t *);
    H5A_operator2_t callback_op = va_arg(arguments, H5A_operator2_t);
    void *op_data = va_arg(arguments, void *);

    if (callback_op) {
      if (op_data == NULL) {
        ADIOS_VOL_WARN("no gid passed to attr iter callback ");
        ret_value = FAIL;
      } else {
        hid_t owner_hid = *((hid_t *)op_data);
        if (var->name == NULL) {
          struct adios_index_comp_struct_v1 *attrs_ptr =
              var->fileReader->attrs_root;

          while (attrs_ptr != NULL) {
            callback_op(owner_hid, attrs_ptr->name, NULL, NULL);
            attrs_ptr = attrs_ptr->next;
          }

          // return ret_value;
        } else { // only return timestep info for variables
          callback_op(owner_hid, ADIOS_TIMESTEP, NULL, NULL);
        }
      }
    }

    break;
  }
  default:
    ADIOS_VOL_NOT_SUPPORTED_ERR("unsupported H5VL_adios_attr_specific flag");
    ret_value = FAIL;
  }
  return ret_value;
}

static void *H5VL_adios_dataset_open(void *obj, H5VL_loc_params_t loc_params,
                                     const char *name, hid_t dapl_id,
                                     hid_t dxpl_id, void **req) {
  REQUIRE_NOT_NULL_ERR(name, NULL);

  struct adios_index_comp_struct_v1 *varInAdios = NULL;
  H5VL_adios_var_t *var = NULL;
  if (H5I_FILE == loc_params.obj_type) {
    /*
     */
    H5VL_adios_t *f = (H5VL_adios_t *)obj;
    int ts = 0;
    getDatasetTimeStep(name, &ts);

    var = createVar(f->bpFileReader, name);
    return var;
  } else if (H5I_GROUP == loc_params.obj_type) {
    // H5VL_adios_t *f = (H5VL_adios_t *)obj;
    var = (H5VL_adios_var_t *)obj;
    varInAdios = adios_get_var_byname(var->fileReader, var->name);
    var->curr_timestep = atoi(name);
  } else {
    ADIOS_VOL_WARN("owner of dataset needs to be either file/group");
    return NULL;
  }

  /*
   */
  assignDim(varInAdios, var);
  return var;
}

void readAll(H5_adios_bb_sel *bb, size_t dataTypeSize,
             H5_bp_request_t *adiosReq, H5VL_adios_var_t *currDataset,
             H5_posMap *sourceSelOrderInC, H5_posMap *targetSelOrderInC,
             void *buf) {
  hsize_t npoints = 1;
  int i = 0;
  uint64_t bb_start[bb->ndim], bb_count[bb->ndim];
  for (i = 0; i < bb->ndim; i++) {
    bb_count[i] = currDataset->dimInfo->dims[i];
    npoints *= bb_count[i];
    bb_start[i] = 0;
  }
  bb->start = bb_start;
  bb->count = bb_count;
  adiosReq->sel = bb;

  if (targetSelOrderInC != NULL) {
    adiosReq->output_data = SAFE_MALLOC(dataTypeSize * npoints);
    read_var_bb(currDataset->fileReader, adiosReq);
    assignToMemSpace(sourceSelOrderInC, targetSelOrderInC, npoints,
                     dataTypeSize, (char *)(adiosReq->output_data),
                     (char *)buf);
    SAFE_FREE(adiosReq->output_data);
  } else {
    read_var_bb(currDataset->fileReader, adiosReq);
  }
}

static herr_t H5VL_adios_dataset_read(void *dset, hid_t mem_type_id,
                                      hid_t mem_space_id, hid_t file_space_id,
                                      hid_t plist_id, void *buf, void **req) {
  H5VL_adios_var_t *currDataset = (H5VL_adios_var_t *)dset;

  REQUIRE_NOT_NULL(currDataset->fileReader);

  size_t dataTypeSize = H5Tget_size(mem_type_id);

  H5S_sel_type memSelType = H5S_SEL_ALL;
  if (mem_space_id != 0) {
    memSelType = H5Sget_select_type(mem_space_id);
    // REQUIRE_SUCC_MSG( (memSelType == H5S_SEL_ALL), -1, "Only handles default
    // memspace for now. ");
  }

  H5S_sel_type h5selType = H5S_SEL_ALL;
  if (file_space_id != 0)
    h5selType = H5Sget_select_type(file_space_id);

  if ((memSelType == H5S_SEL_NONE) || (h5selType == H5S_SEL_NONE)) {
    // no thing was selected, do nothing
    return 0;
  }

  H5_bp_request_t *adiosReq =
      (H5_bp_request_t *)SAFE_CALLOC(1, sizeof(H5_bp_request_t));

  adiosReq->output_data = buf; // buf size >= size of data read.  "=" if
                               // (memSelcType == H5S_SEL_ALL)
  adiosReq->from_step = currDataset->curr_timestep;
  adiosReq->nsteps = 1;
  adiosReq->var =
      adios_get_var_byname(currDataset->fileReader, currDataset->name);

  if (currDataset->dimInfo->ndim == 0) {
    read_scalar(adiosReq, 0, checkSubfile(currDataset->fileReader),
                currDataset->fileReader);
    SAFE_FREE(adiosReq);
    return 1;
  }

  // now ndim >= 1
  H5_adios_bb_sel *bb =
      (H5_adios_bb_sel *)SAFE_CALLOC(1, sizeof(H5_adios_bb_sel));

  bb->ndim = currDataset->dimInfo->ndim;

  H5_posMap *sourceSelOrderInC = NULL;
  GetSelOrder(file_space_id, h5selType, &sourceSelOrderInC);

  // H5_posMap* sourceSelCopy = sourceSelOrderInC;
  H5_posMap *targetSelOrderInC = NULL;
  GetSelOrder(mem_space_id, memSelType, &targetSelOrderInC);
  // H5_posMap* targetSelCopy = targetSelOrderInC;

  if (h5selType == H5S_SEL_ALL) {
    readAll(bb, dataTypeSize, adiosReq, currDataset, sourceSelOrderInC,
            targetSelOrderInC, buf);
  } else if (h5selType == H5S_SEL_HYPERSLABS) {
    hsize_t npoints = H5Sget_select_npoints(file_space_id);
    int i = 0, n = 0;

    hsize_t varSize = 1;
    for (i = 0; i < bb->ndim; i++)
      varSize *= currDataset->dimInfo->dims[i];

    if (varSize == npoints) {
      readAll(bb, dataTypeSize, adiosReq, currDataset, sourceSelOrderInC,
              targetSelOrderInC, buf);
    } else {
      // create adios block selection:
      hsize_t nblocks = H5Sget_select_hyper_nblocks(file_space_id);

      H5S_t *space = (H5S_t *)H5I_object_verify(file_space_id, H5I_DATASPACE);
      int ndims = bb->ndim;

      hsize_t *blockinfo =
          (hsize_t *)SAFE_MALLOC(sizeof(hsize_t) * 2 * ndims * nblocks);
      herr_t status = H5Sget_select_hyper_blocklist(file_space_id, (hsize_t)0,
                                                    nblocks, blockinfo);

      adiosReq->output_data = SAFE_MALLOC(dataTypeSize * npoints);
      void *dup = adiosReq->output_data;

      for (n = 0; n < nblocks; n++) {
        uint64_t blockSize = 1;
        uint64_t h_start[ndims], h_count[ndims];
        for (i = 0; i < ndims; i++) {
          int pos = 2 * ndims * n;
          h_start[i] = blockinfo[pos + i];
          h_count[i] = blockinfo[pos + ndims + i] - h_start[i] + 1;
          blockSize *= h_count[i];
        }

        bb->start = h_start;
        bb->count = h_count;
        adiosReq->sel = bb;

        if (n > 0) {
          adiosReq->output_data += blockSize * dataTypeSize; // move up
        }

        read_var_bb(currDataset->fileReader, adiosReq);
      }

      // assignToMemSpace(sourceSelOrderInC, targetSelOrderInC, npoints,
      // dataTypeSize, (char*)(adiosReq->output_data), (char*)buf);

      assignToMemSpace(sourceSelOrderInC, targetSelOrderInC, npoints,
                       dataTypeSize, (char *)(dup), (char *)buf);
      SAFE_FREE(dup);
      // qsort();
      SAFE_FREE(blockinfo);
    }
  } else if (h5selType == H5S_SEL_POINTS) {
    REQUIRE_SUCC_MSG((memSelType != H5S_SEL_POINTS), -1,
                     "Only handles bulk memspace for now. ");
  }

  SAFE_FREE(sourceSelOrderInC);
  SAFE_FREE(targetSelOrderInC);

  free(bb);
  free(adiosReq);
  return 1;
}

/*
static herr_t H5VL_adios_dataset_write(void *dset, hid_t mem_type_id,
                                       hid_t mem_space_id, hid_t file_space_id,
                                       hid_t plist_id, const void *buf,
                                       void **req) {
  H5VL_adios_t *d = (H5VL_adios_t *)dset;

  printf("------- ADIOS H5Dwrite\n");
  return 1;
}
*/

static herr_t
// H5VL_adios_dataset_close(void *dset, void** req)
H5VL_adios_dataset_close(void *dset, hid_t dxpl_id, void **req) {
  H5VL_adios_var_t *d = (H5VL_adios_var_t *)dset;

  // free(d);

  // no need to free d, it is freed in H5Gclose()
  return 1;
}

//
//
//

static herr_t H5VL_native_datatype_get(void *obj, H5VL_datatype_get_t get_type,
                                       hid_t dxpl_id, void H5_ATTR_UNUSED **req,
                                       va_list arguments) {
  // needed by h5ls
  // calling h5ls => H5Oopen => ... H5VL_register_id( only accepts H5I_DATATYPE)
  // and which in turn asked to construct data type
  /* frame #0: 0x00000001004954d5
    h5ls_bp`H5T_construct_datatype(dt_obj=0x0000000103e00000) + 421 at
    H5Tcommit.c:927 frame #1: 0x00000001005e16c5
    h5ls_bp`H5VL_register_id(type=H5I_DATATYPE, object=0x0000000100f02920,
    vol_plugin=0x0000000101948760, app_ref=true) + 789 at H5VLint.c:169 frame
    #2: 0x00000001002ecf4d h5ls_bp`H5Oopen(loc_id=72057594037927936, name="/",
    lapl_id=792633534417207296) + 1565 at H5O.c:266
   */

  return 1;
}

//
//
//

static const H5VL_class_t H5VL_adios_def = {
    0,
    ADIOS,
    VOLNAME,         /* name */
    H5VL_adios_init, /* initialize */
    H5VL_adios_term, /* terminate */
    sizeof(hid_t),   /* size of vol info */
    NULL,            /* copy vol info */
    NULL,            /* free vol info copy*/
    {
        /* attribute_cls */
        NULL, // H5VL_adios_attr_create,                /* create */
        H5VL_adios_attr_open, /* open */
        H5VL_adios_attr_read, /* read */
        NULL, // H5VL_adios_attr_write,                 /* write */
        H5VL_adios_attr_get,      /* get */
        H5VL_adios_attr_specific, /* specific */
        NULL,                     // H5VL_adios_attr_optional,  /*optional*/
        H5VL_adios_attr_close     /* close */
    },
    {
        /* dataset_cls */
        NULL,                    // H5VL_adios_dataset_create, /* create */
        H5VL_adios_dataset_open, /* open */
        H5VL_adios_dataset_read, /* read */
        NULL,                    // H5VL_adios_dataset_write,  /* write */
        H5VL_adios_dataset_get,  /* get properties*/
        NULL,                    // H5VL_adios_dataset_specific
        NULL,                    // optional
        H5VL_adios_dataset_close /* close */
    },
    {
        /* datatype_cls */
        NULL, // H5VL_adios_datatype_commit,           /* commit */
        NULL, // H5VL_adios_datatype_open,             /* open */
        H5VL_adios_datatype_get, /* get_size */
        NULL,                    // H5VL_adios_datatype_close /* close */
    },
    {
        /* file_cls */
        NULL,                 // no need H5VL_adios_file_create, /* create */
                              // NULL,
        H5VL_adios_file_open, /* open */
        H5VL_adios_file_get,  /* get */
        NULL, // H5VL_adios_file_misc,                  /* misc */
        H5VL_adios_file_optional, /* optional */
        H5VL_adios_file_close     /* close */
    },
    {
        /* group_cls */
        NULL,                  // H5VL_adios_group_create, /* create */
        H5VL_adios_group_open, /* open */
        H5VL_adios_group_get,  /* get  */
        NULL,                  /* specific */
        NULL,                  /* optional */
        H5VL_adios_group_close /* close */
    },
    {
        /* link_cls */
        NULL, // H5VL_adios_link_create,                /* create */
        NULL, // H5VL_adios_link_copy,                  /* copy   */
        NULL, // H5VL_adios_link_move,                  /* move   */
        H5VL_adios_link_get,      /* get    */
        H5VL_adios_link_specific, /* iterate etc*/
        NULL, // H5VL_adios_link_remove                 /* remove */
    },
    {
        /* object_cls */
        H5VL_adios_object_open, /* open */
        NULL, // H5VL_adios_object_copy,                /* copy */
        NULL, // H5VL_adios_object_get,                 /* get */
        NULL, // H5VL_adios_object_specific, /* specific */
        H5VL_adios_object_optional, /* optional */
    }};

herr_t H5Pset_fapl_adios(hid_t *acc_tpl, hid_t *fapl_id) {
  H5VL_ADIOS_g = H5VLregister_driver(&H5VL_adios_def);
  hid_t plugin_id = H5VLget_driver_id(VOLNAME);
  H5VLinitialize(plugin_id, H5P_DEFAULT);
  H5VLclose(plugin_id);

  hid_t ret_value = H5Pset_vol(*acc_tpl, H5VL_ADIOS_g, fapl_id);

#ifdef NEVER
  if (H5Pget_driver(*acc_tpl) == H5FD_MPIO) {
    MPI_Comm comm;
    MPI_Info info;
    H5Pget_fapl_mpio(*acc_tpl, &comm, &info);
    int rank;
    MPI_Comm_rank(comm, &rank);
    printf(" ... parallel ..rank=%d\n", rank);
  } else {
    printf(" ... serial ..\n");
  }
#else
  if (H5Pget_driver(*acc_tpl) == H5FD_MPIO) {
    H5Pget_fapl_mpio(*acc_tpl, &H5VL_MPI_COMM, &H5VL_MPI_INFO);
    int rank, size;
    MPI_Comm_rank(H5VL_MPI_COMM, &rank);
    MPI_Comm_size(H5VL_MPI_COMM, &size);
    if (rank == 0)
      printf(" ... parallel ..size=%d\n", size);
  } else {
    // printf(" ... serial ..\n");
  }
#endif

  return ret_value;
}

herr_t H5P_unset_adios() { return H5VLunregister_driver(H5VL_ADIOS_g); }
