#include "reader.h"
#include <errno.h> // errno
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> // read()

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#if defined(__APPLE__) || defined(__CYGWIN__)
#define O_LARGEFILE 0
#endif

void not_supported(const char *msg) {
  printf("%s is not yet supported in ADIOS VOL..\n", msg);
}

void swap_endian(void *data, int numBytes, enum ADIOS_YES_NOFLAG changeMe) {
  if (changeMe != ADIOS_YES)
    return;

  if (numBytes == 4) {
    uint32_t d = *(uint32_t *)data;
    *(uint32_t *)data = ((d & 0x000000FF) << 24) + ((d & 0x0000FF00) << 8) +
                        ((d & 0x00FF0000) >> 8) + ((d & 0xFF000000) >> 24);
  } else if (numBytes == 8) {
    uint64_t d = *(uint64_t *)data;
    *(uint64_t *)data =
        ((d & 0x00000000000000FF) << 56) + ((d & 0x000000000000FF00) << 40) +
        ((d & 0x0000000000FF0000) << 24) + ((d & 0x00000000FF000000) << 8) +
        ((d & 0x000000FF00000000LL) >> 8) + ((d & 0x0000FF0000000000LL) >> 24) +
        ((d & 0x00FF000000000000LL) >> 40) + ((d & 0xFF00000000000000LL) >> 56);
  } else if (numBytes == 2) {
    uint16_t d = *(uint16_t *)data;
    *(uint16_t *)data = d >> 8 | d << 8;
  } else if (numBytes == 16) // adios_long_double
  {
    uint64_t d = *(uint64_t *)data;
    *(uint64_t *)data =
        ((d & 0x00000000000000FF) << 56) + ((d & 0x000000000000FF00) << 40) +
        ((d & 0x0000000000FF0000) << 24) + ((d & 0x00000000FF000000) << 8) +
        ((d & 0x000000FF00000000LL) >> 8) + ((d & 0x0000FF0000000000LL) >> 24) +
        ((d & 0x00FF000000000000LL) >> 40) + ((d & 0xFF00000000000000LL) >> 56);

    d = *((uint64_t *)data + 1);
    d = ((d & 0x00000000000000FF) << 56) + ((d & 0x000000000000FF00) << 40) +
        ((d & 0x0000000000FF0000) << 24) + ((d & 0x00000000FF000000) << 8) +
        ((d & 0x000000FF00000000LL) >> 8) + ((d & 0x0000FF0000000000LL) >> 24) +
        ((d & 0x00FF000000000000LL) >> 40) + ((d & 0xFF00000000000000LL) >> 56);

    *((uint64_t *)data + 1) = *(uint64_t *)data;
    *(uint64_t *)data = d;
  }
}

uint64_t adios_get_type_size(enum ADIOS_DATATYPES type, const void *var) {
  switch (type) {
  case adios_byte:
  case adios_unsigned_byte:
    return 1;

  case adios_string:
    if (!var)
      return 0;
    else
      return strlen((const char *)var);

  case adios_string_array:
    return sizeof(char *);

  case adios_short:
  case adios_unsigned_short:
    return 2;

  case adios_integer:
  case adios_unsigned_integer:
    return 4;

  case adios_long:
  case adios_unsigned_long:
    return 8;

  case adios_real:
    return 4;

  case adios_double:
    return 8;

  case adios_long_double:
    return 16;

  case adios_complex:
    return 2 * 4;

  case adios_double_complex:
    return 2 * 8;

  default:
    SHOW_ERROR_MSG("unknown adios type %d found.", type);
    return 0;
  }
}

const char *bp_value_to_string(enum ADIOS_DATATYPES type, void *data) {
  static char s[100];
  s[0] = 0;

  switch (type) {
  case adios_unsigned_byte:
    sprintf(s, "%u", *(((uint8_t *)data)));
    break;

  case adios_byte:
    sprintf(s, "%d", *(((int8_t *)data)));
    break;

  case adios_short:
    sprintf(s, "%hd", *(((int16_t *)data)));
    break;

  case adios_unsigned_short:
    sprintf(s, "%uh", *(((uint16_t *)data)));
    break;
  case adios_integer:
    sprintf(s, "%d", *(((int32_t *)data)));
    break;

  case adios_unsigned_integer:
    sprintf(s, "%u", *(((uint32_t *)data)));
    break;

  case adios_long:
    sprintf(s, "%" PRId64 "", *(((int64_t *)data)));
    break;

  case adios_unsigned_long:
    sprintf(s, "%" PRIu64 "", *(((uint64_t *)data)));
    break;

  case adios_real:
    sprintf(s, "%f", *(((float *)data)));
    break;

  case adios_double:
    sprintf(s, "%le", *(((double *)data)));
    break;

  case adios_long_double:
    sprintf(s, "%Le", *(((long double *)data)));
    break;
  case adios_string:
    sprintf(s, "\"%s\"", ((char *)data));
    break;

  case adios_string_array:
    // we expect here one of the array elements, which has char * type
    sprintf(s, "\"%s\"", (*(char **)data));
    break;

  case adios_complex:
    /*
    sprintf (s, "(%f %f)", *(((float *) data) + 0)
             , *(((float *) data) + 1)
             );
    */
    break;

  case adios_double_complex:
    /*
    sprintf (s, "(%lf %lf)", *(((double *) data) + 0)
             , *(((double *) data) + 1)
             );
    */
    break;

  default:
    break;
  }

  return s;
}

// printing out

const char *print_adios_type(int type) {
  switch (type) {
  case adios_unsigned_byte:
    return "unsigned byte";
  case adios_unsigned_short:
    return "unsigned short";
  case adios_unsigned_integer:
    return "unsigned integer";
  case adios_unsigned_long:
    return "unsigned long long";

  case adios_byte:
    return "byte";
  case adios_short:
    return "short";
  case adios_integer:
    return "integer";
  case adios_long:
    return "long long";

  case adios_real:
    return "real";
  case adios_double:
    return "double";
  case adios_long_double:
    return "long double";

  case adios_string:
    return "string";
  case adios_string_array:
    return "string array";
  case adios_complex:
    return "complex";
  case adios_double_complex:
    return "double complex";

  default: {
    SHOW_ERROR_MSG("(unknown type: %d)", type);

    static char buf[50];
    sprintf(buf, "(unknown type: %d)", type);
    return buf;
  }
  }
}

//
// returns bytes read
// var hold the value that was read
//
int readChar(H5_bp_file_t *bpFile, char *var) {
  //*var =  (uint8_t) *(bpFile->scratch + bpFile->offset);
  *var = *(bpFile->scratch + bpFile->offset);
  bpFile->offset += 1;
  return 1;
}

// int readInt(int bytes, H5_bp_file_t* bpFile, uint64_t* pos, uint64_t* var)
void readInt8(H5_bp_file_t *bpFile, uint64_t *var) {
  *var = *(uint64_t *)(bpFile->scratch + bpFile->offset);
  bpFile->offset += 8;

  swap_endian(var, 8, bpFile->change_endianness);
}

void readInt4(H5_bp_file_t *bpFile, uint32_t *var) {
  *var = *(uint32_t *)(bpFile->scratch + bpFile->offset);
  bpFile->offset += 4;

  swap_endian(var, 4, bpFile->change_endianness);
}

void readInt2(H5_bp_file_t *bpFile, uint16_t *var) {
  *var = *(uint16_t *)(bpFile->scratch + bpFile->offset);
  bpFile->offset += 2;

  swap_endian(var, 2, bpFile->change_endianness);
}

void readInt1(H5_bp_file_t *bpFile, uint8_t *var) {
  *var = *(uint8_t *)(bpFile->scratch + bpFile->offset);
  bpFile->offset += 1;
}

void readString(char **strVal, int strLen, H5_bp_file_t *bpFile) {
  // read group_name;
  *strVal = (char *)SAFE_MALLOC(strLen + 1);
  (*strVal)[strLen] = '\0';
  memcpy(*strVal, bpFile->scratch + bpFile->offset, strLen);

  bpFile->offset += strLen;
}

// simplied bp_read_data_from_buffer
double read_double(H5_bp_file_t *bpFile) {
  int doubleSize = 8;
  double result = 0;
  memcpy(&result, (bpFile->scratch + bpFile->offset), doubleSize);
  bpFile->offset += doubleSize;
  return result;
}

uint64_t adios_get_stat_size(void *data, enum ADIOS_DATATYPES type,
                             enum ADIOS_STAT stat_id) {
  uint64_t size = 0;

  switch (type) {
  case adios_complex:
    switch (stat_id) {
    case ADIOS_STATISTIC_min:
    case ADIOS_STATISTIC_max:
    case ADIOS_STATISTIC_sum:
    case ADIOS_STATISTIC_sum_square:
      return adios_get_type_size(adios_double, "");

    case ADIOS_STATISTIC_finite:
      return adios_get_type_size(adios_byte, "");

    case ADIOS_STATISTIC_cnt:
      return adios_get_type_size(adios_unsigned_integer, "");

    case ADIOS_STATISTIC_hist:
      // this is not supported
      return 0;
    }
  case adios_double_complex:
    switch (stat_id) {
    case ADIOS_STATISTIC_min:
    case ADIOS_STATISTIC_max:
    case ADIOS_STATISTIC_sum:
    case ADIOS_STATISTIC_sum_square:
      return adios_get_type_size(adios_long_double, "");

    case ADIOS_STATISTIC_finite:
      return adios_get_type_size(adios_byte, "");

    case ADIOS_STATISTIC_cnt:
      return adios_get_type_size(adios_unsigned_integer, "");

    case ADIOS_STATISTIC_hist:
      // this is not supported
      return 0;
    }
  default: {
    switch (stat_id) {
    case ADIOS_STATISTIC_min:
    case ADIOS_STATISTIC_max:
      return adios_get_type_size(type, "");
    case ADIOS_STATISTIC_sum:
    case ADIOS_STATISTIC_sum_square:
      return adios_get_type_size(adios_double, "");

    case ADIOS_STATISTIC_finite:
      return adios_get_type_size(adios_byte, "");

    case ADIOS_STATISTIC_cnt:
      return adios_get_type_size(adios_unsigned_integer, "");

    case ADIOS_STATISTIC_hist: {
      struct adios_hist_struct *hist = (struct adios_hist_struct *)data;
      size +=
          adios_get_type_size(adios_unsigned_integer, ""); // Number of breaks
      size += adios_get_type_size(adios_double, "");       // Min
      size += adios_get_type_size(adios_double, "");       // Max
      size += ((hist->num_breaks + 1) *
               adios_get_type_size(adios_unsigned_integer, ""));
      size += (hist->num_breaks * adios_get_type_size(adios_double, ""));
      return size;
    }
    default:
      return 0;
    }
  }
  }
}

int adios_get_var_nsteps(struct adios_index_comp_struct_v1 *var_root) {
  uint64_t i;
  int nsteps = 0;
  int prev_step = -1;

  for (i = 0; i < var_root->characteristics_count; i++) {
    if (var_root->characteristics[i].time_index != prev_step) {
      prev_step = var_root->characteristics[i].time_index;
      nsteps++;
    }
  }

  return nsteps;
}

int adios_get_ch_idx(struct adios_index_comp_struct_v1 *var,
                     int curr_timestep) {
  uint64_t i;
  for (i = 0; i < var->characteristics_count; i++) {
    // time_index is always starting at 1, so need to add 1 before comparing.
    if (var->characteristics[i].time_index == curr_timestep + 1)
      return i;
  }
  return -1;
}

int adios_get_var_ndim(H5_bp_file_t *fp, int idx,
                       struct adios_index_comp_struct_v1 *varInAdios,
                       unsigned int *inqIsGlobal) {
  struct adios_index_characteristic_dims_struct_v1 *dims =
      &(varInAdios->characteristics[idx].dims);
  int ndim = dims->count;

  if (ndim == 0) {
    return 0;
  }

  int has_time_index_characteristic =
      fp->bp_version & ADIOS_VERSION_HAVE_TIME_INDEX_CHARACTERISTIC;
  uint64_t ldims[ndim], gdims[ndim], offsets[ndim];
  int is_global = adios_get_var_dimensions(dims, -1, ldims, gdims, offsets,
                                           fp->adios_host_language_fortran);

  if (!is_global) { // local array
    int j = 0, i = 0;

    int n = ndim;
    for (i = 0; i < n; i++) {
      // size of time dimension is always registered as 1 for an array
      if (ldims[i] == 1 && varInAdios->characteristics_count > 1) {
        ndim = ndim - 1;
      }
    }
  } else {
    // global array:
    //     time dimension: ldims=1, gdims=0
    //     in C array, it can only be the first dim
    //     in Fortran array, it can only be the last dim
    //     (always the slowest changing dim)
    //
    if (gdims[ndim - 1] == 0) { // with time
      enum ADIOS_YES_NOFLAG file_is_fortran = fp->adios_host_language_fortran;
      if (file_is_fortran == ADIOS_YES) {
        not_supported("Fortran conversion");
      } else {
        // first dimension is the time (C array)
        //   * ldims[0] = 1 but gdims does not contain time info and
        //   * gdims[0] is 1st data dimension and
        //   * gdims is shorter by one value than ldims in case of C.
        //   * Therefore, gdims[*ndim-1] = 0 if there is a time dimension.
        //
        // error check
        if (ndim > 1 && ldims[0] != 1) {
          SHOW_ERROR_MSG("Inconsistent indication on local dimension");
        }
      }

      ndim = ndim - 1;
    }
  }

  if (inqIsGlobal)
    *inqIsGlobal = is_global;

  return ndim;
}

int adios_get_var_dimensions(
    const struct adios_index_characteristic_dims_struct_v1 *dims, int ndim,
    uint64_t *ldims, uint64_t *gdims, uint64_t *offsets,
    enum ADIOS_YES_NOFLAG is_fortran_order) {
  int is_global = 0; // global array or just an array written by one process?
  int k;

  // int ndim = dims->count; // see adios_get_var_ndim
  if (ndim == -1) {
    ndim = dims->count; // default
  }

  for (k = 0; k < ndim; k++) {
    ldims[k] = dims->charDims[k * 3];
    gdims[k] = dims->charDims[k * 3 + 1];
    offsets[k] = dims->charDims[k * 3 + 2];
    is_global = is_global || gdims[k];
  }

  if ((ndim > 0) && (ndim < dims->count)) {
    // adjust the timed characteristic storage ..
    for (k = 0; k < ndim; k++) {
      ldims[k] = dims->charDims[k * 3 + 3];
    }
  }

  //
  // swap
  //
  if (is_fortran_order == ADIOS_YES) {
    uint64_t tmp;
    for (k = 0; k < ndim / 2; k++) {
      tmp = ldims[k];
      ldims[k] = ldims[ndim - 1 - k];
      ldims[ndim - 1 - k] = tmp;

      tmp = gdims[k];
      gdims[k] = gdims[ndim - 1 - k];
      gdims[ndim - 1 - k] = tmp;

      tmp = offsets[k];
      offsets[k] = offsets[ndim - 1 - k];
      offsets[ndim - 1 - k] = tmp;
    }
  }

  return is_global;
}

/* Convert 'step' to time, which is used in ADIOS internals.
 * 'step' should start from 0.
 */
int adios_get_time(struct adios_index_comp_struct_v1 *v, int step) {
  int i = 0;
  int prev_ti = 0, counter = 0;

  while (i < v->characteristics_count) {
    if (v->characteristics[i].time_index != prev_ti) {
      counter++;
      if (counter == (step + 1)) {
        return v->characteristics[i].time_index;
      }
      prev_ti = v->characteristics[i].time_index;
    }

    i++;
  }

  return -1;
}

// Search for the start var index.
int64_t adios_get_var_start_index(struct adios_index_comp_struct_v1 *v, int t) {
  int64_t i = 0;

  while (i < v->characteristics_count) {
    if (v->characteristics[i].time_index == t) {
      return i;
    }

    i++;
  }
  return -1;
}

// Search for the stop var index
int64_t adios_get_var_stop_index(struct adios_index_comp_struct_v1 *v, int t) {
  int64_t i = v->characteristics_count - 1;

  while (i > -1) {
    if (v->characteristics[i].time_index == t) {
      return i;
    }

    i--;
  }

  return -1;
}

int adios_util_check_selection(int ndim, H5_bp_request_t *req, uint64_t *gdims,
                               uint64_t *ldims, uint64_t *offsets,
                               uint64_t *payload_size) {
  int idx_useful = 1;
  uint64_t *start = req->sel->start;
  uint64_t *count = req->sel->count;

  int j = 0;
  int flag;
  for (j = 0; j < ndim; j++) {
    *payload_size *= ldims[j];

    if ((count[j] > gdims[j]) || (start[j] > gdims[j]) ||
        (start[j] + count[j] > gdims[j])) {
      SHOW_ERROR_MSG("Error: Variable (id=%d) out of bound 1("
                     "the data in dimension %d to read is %" PRIu64
                     " elements from index %" PRIu64
                     " but the actual data is [0,%" PRId64 "])\n",
                     req->var->id, j + 1, count[j], start[j], gdims[j] - 1);
      return 0;
    }

    /* check if there is any data in this pg and this dimension to read in */
    flag = (offsets[j] >= start[j] && offsets[j] < (start[j] + count[j])) ||
           (offsets[j] < start[j] &&
            (offsets[j] + ldims[j]) > (start[j] + count[j])) ||
           (((offsets[j] + ldims[j]) > start[j]) &&
            (offsets[j] + ldims[j]) <= (start[j] + count[j]));

    idx_useful = idx_useful && flag;
  }

  return idx_useful;
}

void adios_util_count_partial_pg(const int ndim, uint64_t *size_in_dset,
                                 uint64_t *offset_in_dset,
                                 uint64_t *offset_in_var, const uint64_t *start,
                                 const uint64_t *count, const uint64_t *ldims,
                                 const uint64_t *offsets) {
  memset(size_in_dset, 0, 10 * 8);
  memset(offset_in_dset, 0, 10 * 8);
  memset(offset_in_var, 0, 10 * 8);

  uint64_t isize;
  int i;
  for (i = 0; i < ndim; i++) {
    isize = offsets[i] + ldims[i];
    if (start[i] >= offsets[i]) {
      // head is in
      if (start[i] < isize) {
        if (start[i] + count[i] > isize)
          size_in_dset[i] = isize - start[i];
        else
          size_in_dset[i] = count[i];

        offset_in_dset[i] = start[i] - offsets[i];
        offset_in_var[i] = 0;
      }
    } else {
      // middle is in
      if (isize < start[i] + count[i]) {
        size_in_dset[i] = ldims[i];
      } else {
        // tail is in
        size_in_dset[i] = count[i] + start[i] - offsets[i];
      }

      offset_in_dset[i] = 0;
      offset_in_var[i] = offsets[i] - start[i];
    }
  }
}

//////////
//
// adios supports up to  1 D attributes only
// so: return 0 if scalar
//     otherwise N for the size of 1D dim
int adios_get_attr_ndim(struct adios_index_comp_struct_v1 *attrInAdios) {
  struct adios_index_characteristic_dims_struct_v1 *dims =
      &(attrInAdios->characteristics[0].dims);
  int ndim = dims->count;

  // scalar
  if (ndim == 0) {
    return ndim;
  }

  return dims->charDims[0];
}
