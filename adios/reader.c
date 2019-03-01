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

// this is added for cori
#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif

#define BYTE_ALIGN 8

void change_endianness(void *data, uint64_t slice_size,
                       enum ADIOS_DATATYPES type) {
  printf(" -- TODO: change_endianness() -- \n");
}

int handle_char_value(H5_bp_file_t *b, struct adios_index_comp_struct_v1 **root,
                      enum ADIOS_CHARACTERISTICS c, uint64_t j);
/*
 */

uint8_t adios_get_stat_set_count(enum ADIOS_DATATYPES type) {
  // NCSU - Complex numbers have the statistic related values repeated.
  // NCSU - Holds values for real, imaginary, complex parts
  if (type == adios_complex || type == adios_double_complex)
    return 3;
  return 1;
}

struct adios_index_comp_struct_v1 *adios_get_var_byname(H5_bp_file_t *file,
                                                        const char *varname) {
  struct adios_index_comp_struct_v1 *curr = file->vars_root;
  while (curr) {
    if (strcmp(curr->name, varname) == 0) {
      return curr;
    }
    curr = curr->next;
  }
  return NULL;
}

struct adios_index_comp_struct_v1 *adios_get_attr_byname(H5_bp_file_t *file,
                                                         const char *attrname) {
  struct adios_index_comp_struct_v1 *curr = file->attrs_root;
  while (curr) {
    if (strcmp(curr->name, attrname) == 0) {
      return curr;
    }
    curr = curr->next;
  }
  return NULL;
}

enum ADIOS_DATATYPES adios_transform_get_var_original_type_var_header(
    H5_bp_var_data_header_t *var_header) {
  if (var_header->characteristics.transform.transform_type !=
      /*adios_transform_none*/ 0)
    return var_header->characteristics.transform.pre_transform_type;
  else
    return var_header->type;
}

enum ADIOS_DATATYPES adios_transform_get_var_original_type_index(
    struct adios_index_comp_struct_v1 *var) {
  if (var->characteristics[0].transform.transform_type !=
      /*adios_transform_none*/ 0)
    return var->characteristics[0].transform.pre_transform_type;
  else
    return var->type;
}

// void readHist(H5_bp_file_t* b, struct adios_index_var_struct_v1 ** root,
// uint64_t j)
void readHist(H5_bp_file_t *b,
              struct adios_index_characteristics_hist_struct *hist) {
  uint32_t bi, num_breaks;

  // Getting the number of breaks of histogram
  readInt4(b, &(hist->num_breaks));
  num_breaks = hist->num_breaks;

  hist->min = read_double(b);
  hist->max = read_double(b);

  // Getting the frequencies of the histogram
  hist->frequencies = SAFE_MALLOC(
      (num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
  memcpy(hist->frequencies, (b->scratch + b->offset),
         (num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));

  if (b->change_endianness == ADIOS_YES) {
    for (bi = 0; bi <= num_breaks; bi++) {
      // swap_32(hist->frequencies[bi]);
      not_supported("swap_32\n");
    }
  }
  b->offset += 4 * (num_breaks + 1);

  // Getting the breaks of the histogram
  hist->breaks =
      SAFE_MALLOC(num_breaks * adios_get_type_size(adios_double, ""));
  memcpy(hist->breaks, (b->scratch + b->offset),
         num_breaks * adios_get_type_size(adios_double, ""));
  if (b->change_endianness == ADIOS_YES) {
    for (bi = 0; bi < num_breaks; bi++)
      // swap_64(hist->breaks[bi]);
      not_supported("swap_64\n");
  }
  b->offset += 8 * num_breaks;
}

static uint64_t read64(int f, char *buff, uint64_t readsize) {
  // if we need to read > 2 GB, need to do it in parts
  // since count is limited to MAX_MPIWRITE_SIZE (signed 32-bit max).
  uint64_t bytes_read = 0;
  int32_t to_read = 0;
  int err = 0;
  const size_t MAX_READ_SIZE = 0x7ffff000; // = 2,147,479,552

  while (bytes_read < readsize && !err) {
    if (readsize - bytes_read > MAX_READ_SIZE)
      to_read = MAX_READ_SIZE;
    else
      to_read = readsize - bytes_read;

    ssize_t actual_read_bytes = read(f, buff + bytes_read, to_read);
    REQUIRE_SUCC_MSG((actual_read_bytes >= 0), 0,
                     "Error while reading from file. %d bytes:  '%s'\n",
                     to_read, strerror(errno));
    /*
     */
    REQUIRE_SUCC_MSG((actual_read_bytes == to_read), 0,
                     "Error while reading from file tried to read %d bytes but "
                     "only got %ld bytes\n",
                     to_read, actual_read_bytes);
    /*
     */
    bytes_read += actual_read_bytes;
  }
  return bytes_read;
}

static void realloc_aligned(H5_bp_file_t *b, uint64_t size) {
  // b->allocated_buff_ptr = realloc (b->allocated_buff_ptr, size + BYTE_ALIGN -
  // 1);
  b->allocated_buff_ptr =
      (char *)(SAFE_REALLOC(b->allocated_buff_ptr, size + BYTE_ALIGN - 1));

  if (!b->allocated_buff_ptr) {
    printf("BP_V1: Cannot allocate %" PRIu64 "\n", size);

    b->scratch = 0;
    b->length = 0;

    return;
  }
  uint64_t p = (uint64_t)b->allocated_buff_ptr;

  b->scratch = (char *)((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
  b->length = size;
}

int adios_parse_var_data_header_v1(H5_bp_file_t *b,
                                   H5_bp_var_data_header_t *var_header) {
  uint64_t diff = b->length - b->offset;

  REQUIRE_SUCC_MSG((diff < 23), 1,
                   "adios_parse_vars_header_v1() requires a buffer of at least "
                   "23 bytes Only %llu  were provided\n",
                   diff);

  uint64_t initial_offset = b->offset; // save to calc payload size
  uint64_t length_of_var;
  uint16_t len;
  uint8_t flag;

  readInt8(b, &length_of_var);

  // validate remaining length
  readInt4(b, &(var_header->id));

  readInt2(b, &len);
  readString(&(var_header->name), len, b);

  readInt2(b, &len);
  readString(&(var_header->path), len, b);

  readInt1(b, &flag);
  var_header->type = (enum ADIOS_DATATYPES)flag;

  char yn_flag;
  readChar(b, &yn_flag);
  var_header->is_dim = (yn_flag == 'y' ? ADIOS_YES : ADIOS_NO);

  int i;
  uint8_t dims_count;
  uint16_t dims_length;
  uint8_t characteristics_count;
  uint32_t characteristics_length;

  // validate remaining length

  readInt1(b, &dims_count);

  readInt2(b, &dims_length);

  var_header->dims = 0;
  struct adios_dimension_struct_v1 **root = &var_header->dims;

  for (i = 0; i < dims_count; i++) {
    if (!*root) {
      *root = (struct adios_dimension_struct_v1 *)malloc(
          sizeof(struct adios_dimension_struct_v1));
      (*root)->next = 0;
    }

    (*root)->dimension.rank = 0;
    (*root)->dimension.var_id = 0;
    (*root)->dimension.is_time_index = ADIOS_NO;

    (*root)->global_dimension.rank = 0;
    (*root)->global_dimension.var_id = 0;
    (*root)->global_dimension.is_time_index = ADIOS_NO;

    (*root)->local_offset.rank = 0;
    (*root)->local_offset.var_id = 0;
    (*root)->local_offset.is_time_index = ADIOS_NO;

    (*root)->local_offset.rank = 0;
    (*root)->local_offset.var_id = 0;
    (*root)->local_offset.is_time_index = ADIOS_NO;

    readInt1(b, &flag);
    if (flag == 'y') {
      (*root)->dimension.rank = 0;
      readInt4(b, &((*root)->dimension.var_id));
      if ((*root)->dimension.var_id == 0)
        (*root)->dimension.is_time_index = ADIOS_YES;
    } else {
      readInt8(b, &((*root)->dimension.rank));
      (*root)->dimension.var_id = 0;
    }

    readInt1(b, &flag);
    if (flag == 'y') {
      (*root)->global_dimension.rank = 0;
      readInt4(b, &((*root)->global_dimension.var_id));
      if ((*root)->global_dimension.var_id == 0)
        (*root)->global_dimension.is_time_index = ADIOS_YES;
    } else {
      readInt8(b, &((*root)->global_dimension.rank));
      (*root)->global_dimension.var_id = 0;
    }

    readInt1(b, &flag);
    if (flag == 'y') {
      (*root)->local_offset.rank = 0;
      readInt4(b, &((*root)->local_offset.var_id));
      if ((*root)->local_offset.var_id == 0)
        (*root)->local_offset.is_time_index = ADIOS_YES;
    } else {
      readInt8(b, &((*root)->local_offset.rank));
      (*root)->local_offset.var_id = 0;
    }

    root = &(*root)->next;
  }

  readInt1(b, &characteristics_count);
  readInt4(b, &characteristics_length);

  uint64_t size = adios_get_type_size(var_header->type, "");

  var_header->characteristics.offset = 0;
  var_header->characteristics.payload_offset = 0;
  var_header->characteristics.stats = 0;
  var_header->characteristics.value = 0;
  var_header->characteristics.dims.count = 0;
  var_header->characteristics.dims.charDims = 0;
  // NCSU - Initialize statistics fields
  var_header->characteristics.bitmap = 0;
  var_header->characteristics.stats = 0;
  // NCSU ALACRITY-ADIOS - Initialize transform field
  // not supporting transform - Junmin
  // adios_transform_init_transform_characteristic(&var_header->characteristics.transform);

  for (i = 0; i < characteristics_count; i++) {
    uint8_t flag;
    readInt1(b, &flag);
    enum ADIOS_CHARACTERISTICS c;

    c = (enum ADIOS_CHARACTERISTICS)flag;

    switch (c) {
    case ADIOS_CHARACTERISTIC_offset:
      readInt8(b, &var_header->characteristics.offset);
      break;

    case ADIOS_CHARACTERISTIC_payload_offset:
      readInt8(b, &var_header->characteristics.payload_offset);
      break;

    // NCSU ALACRITY-ADIOS - Read in transform type field
    case ADIOS_CHARACTERISTIC_transform_type:
      // adios_transform_deserialize_transform_characteristic(&var_header->characteristics.transform,
      // b);
      not_supported("adios_transform_deserialize_transform_characteristic");
      break;

    // NCSU - Read in bitmap
    case ADIOS_CHARACTERISTIC_bitmap:
      readInt4(b, &var_header->characteristics.bitmap);
      break;

    // NCSU - Read in the statistics
    case ADIOS_CHARACTERISTIC_stat: {
      enum ADIOS_DATATYPES original_var_type =
          adios_transform_get_var_original_type_var_header(var_header);
      handle_char_stat(b, original_var_type, &var_header->characteristics);

      break;
    }

    case ADIOS_CHARACTERISTIC_dimensions: {
      // uint8_t dim_count;
      b->offset += 1; // why?
      uint16_t dim_length;
      readInt2(b, &dim_length);

      var_header->characteristics.dims.charDims = malloc(dim_length);
      memcpy(var_header->characteristics.dims.charDims,
             (b->scratch + b->offset), dim_length);
      if (b->change_endianness == ADIOS_YES) {
        uint16_t di = 0;
        uint16_t dim_num = dim_length / 8;
        for (di = 0; di < dim_num; di++) {
          // swap_64((var_header->characteristics.dims.dims)[di]);
          not_supported("swap_64()");
        }
      }
      b->offset += dim_length;
      break;
    }

    // NCSU - Adding backward compatibility
    // Reading min and max here implies older bp file
    case ADIOS_CHARACTERISTIC_min:
      readStat(&var_header->characteristics, ADIOS_STATISTIC_min, size, b);
      break;
    case ADIOS_CHARACTERISTIC_max:
      readStat(&var_header->characteristics, ADIOS_STATISTIC_max, size, b);
      break;

    case ADIOS_CHARACTERISTIC_value: {
      uint16_t val_size;
      readInt2(b, &val_size);
      char *val = (char *)(var_header->characteristics.value);
      readString(&val, val_size, b);
      break;
    }

    case ADIOS_CHARACTERISTIC_var_id:
    case ADIOS_CHARACTERISTIC_file_index:
    case ADIOS_CHARACTERISTIC_time_index:
      // these are not in the var header in the PG
      break;
    }
  }

  var_header->payload_size = length_of_var - (b->offset - initial_offset);

  /* At this point we are at the payload, processed all metadata.
 If offset and payload_offset was not read in, let's calculate them here
       b->offset was 0 when reading in the whole PG and start processing with
       parsing the process group.
        initial_offset points to the beginning of this variable's metadata.
        b->offset points to the beginning of this variable's payload.
  */
  if (var_header->characteristics.offset == 0)
    var_header->characteristics.offset = b->read_pg_offset + initial_offset;
  if (var_header->characteristics.payload_offset == 0)
    var_header->characteristics.payload_offset = b->read_pg_offset + b->offset;
  return 0;
}

/*
int adios_parse_vars_header_v1 (H5_bp_file_t* b,
                                H5_bp_vars_header_t* vardata_header)
{
  vars_header->count=0; vars_header->length = 0;

  uint64_t diff = b->length - b->offset;

  REQUIRE_SUCC_MSG((diff< 12), 1, "adios_parse_vars_header_v1() requires a
buffer of at least 12 bytes", "Only %" PRId64, " were provided\n", diff);

  readInt4(b, &vars_header->count);
  readInt8(b, &vars_header->length);

  vars_header->count = *(uint32_t *) (b->buff + b->offset);
  if(b->change_endianness == adios_flag_yes) {
    swap_32(vars_header->count);
  }
  b->offset += 4;
  vars_header->length = *(uint64_t *) (b->buff + b->offset);
  if(b->change_endianness == adios_flag_yes) {
    swap_64(vars_header->length);
  }
  b->offset += 8;

  return 0;
}
*/

/* from adios: emulate MPI_File_read with support for data >2GiB */
void mpi_FILE_READ64(H5_bp_file_t *fh, uint64_t size, MPI_Datatype type) {
  int read_len = 0;
  uint64_t total_read = 0;
  uint64_t to_read = size;
  int count;

  MPI_Status status;
  while (to_read > 0) {
    read_len = (to_read > MAX_MPIWRITE_SIZE) ? MAX_MPIWRITE_SIZE : to_read;
    MPI_File_read(fh->mpi_handle, (char *)fh->scratch + total_read, read_len,
                  type, &status);
    MPI_Get_count(&status, type, &count);
    int err;
    if (count != read_len) {
      SHOW_ERROR_MSG("Need to do multi-read (tried: "
                     "%d read: %d) errno %d\n",
                     read_len, count, errno);
      err = count;
      break;
    }
    total_read += count;
    to_read -= count;
  }
}

MPI_File *adios_get_BP_subfile_handle(H5_bp_file_t *fh, uint32_t file_index) {
  BP_file_handle_list *lst = &fh->subfile_handles; // just for simplifying
                                                   // typing
  // printf ("%s # of handles=%d, search for file_index=%d, fh=%p\n", __func__,
  // lst->n_handles, file_index, fh);
  if (!lst->head)
    return 0;

  struct BP_file_handle *curr = lst->head;

  while (curr) {
    if (curr->file_index == file_index)
      return &curr->fh;

    curr = curr->next;
  }

  return 0;
}

static void close_BP_subfile(struct BP_file_handle *sfh) {
  if (sfh)
    MPI_File_close(&sfh->fh);
}

void add_BP_subfile_handle(H5_bp_file_t *fh, struct BP_file_handle *n) {
  if (!n)
    return;

  BP_file_handle_list *lst = &fh->subfile_handles; // just for simplifying
                                                   // typing

  // link the newcomer as first (head)
  n->next = lst->head;
  if (lst->head) {
    lst->head->prev = n;
  }
  lst->head = n;

  // set tail if this is first element
  if (!lst->tail) {
    lst->tail = n;
  }
  lst->n_handles++;

  // printf ("%s # of handles=%d, now added file_index=%d, fh=%p\n", __func__,
  // lst->n_handles, n->file_index, fh);

  // Don't run out of file descriptiors by keeping the number in check
  if (lst->n_handles > ADIOS_MAX_HANDLES) {
    if (!lst->warning_printed) {
      SHOW_ERROR_MSG("Number of subfiles of file %s opened in a single process "
                     "reached %d which indicates an inefficient reading "
                     "pattern.\n",
                     fh->file_name, lst->n_handles);
      lst->warning_printed = 1;
    }
    // remove the tail
    lst->tail->prev->next = NULL;
    struct BP_file_handle *oldest = lst->tail;
    lst->tail = lst->tail->prev;
    close_BP_subfile(oldest);
    free(oldest);
    lst->n_handles--;
  }
}

int mpi_OPS2_readvalue(H5_bp_file_t *fh, struct adios_index_comp_struct_v1 *v,
                       uint64_t slice_size, uint64_t slice_offset,
                       int64_t currIdx) {
  int succ = 1;

  realloc_aligned(fh, slice_size);
  fh->offset = 0;

  MPI_File *sfh =
      adios_get_BP_subfile_handle(fh, v->characteristics[currIdx].file_index);

  if (!sfh) {
    char *ch, *name_no_path, *name;
    struct BP_file_handle *new_h =
        (struct BP_file_handle *)malloc(sizeof(struct BP_file_handle));

    new_h->file_index = v->characteristics[currIdx].file_index;
    new_h->next = 0;
    if ((ch = strrchr(fh->file_name, '/'))) {
      name_no_path = (char *)malloc(strlen(ch + 1) + 1);
      strcpy(name_no_path, ch + 1);
    } else {
      name_no_path = (char *)malloc(strlen(fh->file_name) + 1);
      strcpy(name_no_path, fh->file_name);
    }

    name = (char *)malloc(strlen(fh->file_name) + 5 + strlen(name_no_path) + 1 +
                          10 + 1);
    sprintf(name, "%s.dir/%s.%d", fh->file_name, name_no_path,
            new_h->file_index);
    int err = MPI_File_open(MPI_COMM_SELF, name, MPI_MODE_RDONLY, MPI_INFO_NULL,
                            &new_h->fh);

    REQUIRE_SUCC_MSG((err == 0), !succ, "can not open file %s\n", name);

    add_BP_subfile_handle(fh, new_h);
    sfh = &new_h->fh;

    free(name_no_path);
    free(name);
  }

  MPI_File_seek(*sfh, (MPI_Offset)slice_offset, MPI_SEEK_SET);

  H5_bp_file_t temp;
  temp.mpi_handle = *sfh;
  temp.scratch = fh->scratch;
  mpi_FILE_READ64(&temp, slice_size, MPI_BYTE);
  // mpi_FILE_READ64 (*sfh, fh->scratch, slice_size, MPI_BYTE);

  fh->offset = 0;

  return succ;
}

void mpi_OPS3_readvalue(int64_t currIdx, H5_bp_file_t *fh,
                        struct adios_index_comp_struct_v1 *v) {
  MPI_Status status;
  MPI_File_seek(fh->mpi_handle, (MPI_Offset)v->characteristics[currIdx].offset,
                MPI_SEEK_SET);
  MPI_File_read(fh->mpi_handle, fh->scratch, 8, MPI_BYTE, &status);
  uint64_t tmpcount = *((uint64_t *)fh->scratch);

  realloc_aligned(fh, tmpcount + 8);
  fh->offset = 0;

  MPI_File_seek(fh->mpi_handle,
                (MPI_Offset)(v->characteristics[currIdx].offset), MPI_SEEK_SET);
  mpi_FILE_READ64(fh, tmpcount + 8, MPI_BYTE);
  fh->offset = 0;
  H5_bp_var_data_header_t var_header;
  adios_parse_var_data_header_v1(fh, &var_header);
}

void mpi_OPS1_readvalue(uint64_t slice_size, uint64_t slice_offset,
                        H5_bp_file_t *fh) {
  // MPI_Status status;
  realloc_aligned(fh, slice_size);
  fh->offset = 0;

  MPI_File_seek(fh->mpi_handle, (MPI_Offset)slice_offset, MPI_SEEK_SET);
  mpi_FILE_READ64(fh, slice_size, MPI_BYTE);

  fh->offset = 0;
}

int adios_posix_read_vars_index(H5_bp_file_t *b) {
  uint64_t vars_size = b->attr_index - b->var_index;

  realloc_aligned(b, vars_size);
  b->offset = 0;

  uint64_t r;

  lseek(b->serial_handle, b->var_index, SEEK_SET);
  r = read64(b->serial_handle, b->scratch, vars_size);

  REQUIRE_SUCC_MSG((r == vars_size), -1,
                   "reading vars_index: wanted %" PRIu64 ", read: %" PRIu64
                   "\n",
                   vars_size, r);
  return 0;
}

int adios_posix_read_attributes_index(H5_bp_file_t *b) {
  uint64_t attr_size = b->attrs_end - b->attr_index;

  realloc_aligned(b, attr_size);
  b->offset = 0;
  uint64_t r;

  lseek(b->serial_handle, b->attr_index, SEEK_SET);
  r = read64(b->serial_handle, b->scratch, attr_size);

  REQUIRE_SUCC_MSG((r == attr_size), -1,
                   "reading attributess_index: wanted %" PRIu64
                   ", read: %" PRIu64 "\n",
                   attr_size, r);

  return 0;
}

int adios_parse_attributes_index_v1(
    H5_bp_file_t *b, struct adios_index_comp_struct_v1 **attrs_root) {
  struct adios_index_comp_struct_v1 **root;

  REQUIRE_SUCC_MSG((b->length - b->offset >= 10), 1,
                   "adios_parse_attr_index_v1 requires a buffer of at least 10 "
                   "bytes. Only %" PRId64 " were provided\n",
                   b->length - b->offset);

  root = attrs_root;

  /* BP Format v1: vars_count and attrs_count was 16bit, now it's 32 bit */
  uint32_t attrs_count;
  readInt4(b, &attrs_count);
  uint64_t attrs_length;
  readInt8(b, &attrs_length);

  // validate remaining length
  uint32_t i;
  for (i = 0; i < attrs_count; i++) {
    if (!*root) {
      *root = (struct adios_index_comp_struct_v1 *)SAFE_MALLOC(
          sizeof(struct adios_index_comp_struct_v1));
      (*root)->next = 0;
    }
    (*root)->nelems =
        1; // initialize to 1 in case there will be no dimension characteristic
    uint16_t len;
    uint32_t attr_entry_length;

    readInt4(b, &attr_entry_length);
    readInt4(b, &((*root)->id));

    readInt2(b, &len);
    readString(&((*root)->group_name), len, b);

    readInt2(b, &len);
    readString(&((*root)->name), len, b);

    readInt2(b, &len);
    readString(&((*root)->path), len, b);

    uint8_t flag;
    readInt1(b, &flag);
    (*root)->type = (enum ADIOS_DATATYPES)flag;

    uint64_t characteristics_sets_count;
    readInt8(b, &characteristics_sets_count);
    (*root)->characteristics_count = characteristics_sets_count;
    (*root)->characteristics_allocated = characteristics_sets_count;

    uint64_t j;
    (*root)->characteristics =
        SAFE_MALLOC(characteristics_sets_count *
                    sizeof(struct adios_index_characteristic_struct_v1));

    memset((*root)->characteristics, 0,
           characteristics_sets_count *
               sizeof(struct adios_index_characteristic_struct_v1));

    for (j = 0; j < characteristics_sets_count; j++) {
      uint8_t characteristic_set_count;
      readInt1(b, &characteristic_set_count);
      uint32_t characteristic_set_length;
      readInt4(b, &characteristic_set_length);

      uint8_t item = 0;
      while (item < characteristic_set_count) {
        uint8_t flag;
        readInt1(b, &flag);
        enum ADIOS_CHARACTERISTICS c = (enum ADIOS_CHARACTERISTICS)flag;
        switch (c) {
        case ADIOS_CHARACTERISTIC_value: {
          handle_char_value(b, root, c, j);
          break;
        }
        case ADIOS_CHARACTERISTIC_offset: {
          readInt8(b, &((*root)->characteristics[j].offset));
          break;
        }
        case ADIOS_CHARACTERISTIC_bitmap: {
          readInt4(b, &((*root)->characteristics[j].bitmap));
          break;
        }
        case ADIOS_CHARACTERISTIC_payload_offset: {
          readInt8(b, &((*root)->characteristics[j].payload_offset));
          break;
        }
        case ADIOS_CHARACTERISTIC_file_index: {
          readInt4(b, &((*root)->characteristics[j].file_index));
          break;
        }
        case ADIOS_CHARACTERISTIC_time_index: {
          readInt4(b, &((*root)->characteristics[j].time_index));
          break;
        }
        case ADIOS_CHARACTERISTIC_var_id: {
          /* BP Format v1: varid/attrid was 16bit, now it's 32 bit */
          readInt4(b, &((*root)->characteristics[j].var_id));
          break;
        }
        case ADIOS_CHARACTERISTIC_transform_type: // NCSU ALACRITY-ADIOS -
                                                  // Deserialize transform
                                                  // characteristic
        {
          // adios_transform_deserialize_transform_characteristic(&(*root)->characteristics[j].transform,
          // b);
          not_supported("NCSU related");
          break;
        }
        case ADIOS_CHARACTERISTIC_dimensions: {
          readInt1(b, &((*root)->characteristics[j].dims.count));
          uint16_t dims_length;
          readInt2(b, &dims_length);

          (*root)->characteristics[j].dims.charDims =
              (uint64_t *)malloc(dims_length);
          memcpy((*root)->characteristics[j].dims.charDims,
                 (b->scratch + b->offset), dims_length);
          if (b->change_endianness == ADIOS_YES) {
            int di = 0;
            int dims_num = dims_length / 8;
            for (di = 0; di < dims_num; di++) {
              // swap_64(((*root)->characteristics [j].dims.dims)[di]);
              not_supported("swap_64");
            }
          }
          b->offset += dims_length;
          (*root)->nelems = (*root)->characteristics[j].dims.charDims[0];
          break;
        }

        default:
          break;
        }
        item++;
      }
    }
    root = &(*root)->next;
  }
  return 0;
}

int handle_char_stat(H5_bp_file_t *b, enum ADIOS_DATATYPES original_var_type,
                     struct adios_index_characteristic_struct_v1 *chaStruc) {
  uint8_t k, c, idx;
  uint64_t count = adios_get_stat_set_count(original_var_type);

  chaStruc->stats = SAFE_MALLOC(
      count * sizeof(struct adios_index_characteristics_stat_struct *));

  for (c = 0; c < count; c++) {
    chaStruc->stats[c] =
        SAFE_CALLOC(ADIOS_STAT_LENGTH,
                    sizeof(struct adios_index_characteristics_stat_struct));
    k = idx = 0;
    while (chaStruc->bitmap >> k) {
      chaStruc->stats[c][k].data = 0;

      if ((chaStruc->bitmap >> k) & 1) {
        if (k == ADIOS_STATISTIC_hist) {
          struct adios_index_characteristics_hist_struct *hist = SAFE_MALLOC(
              sizeof(struct adios_index_characteristics_hist_struct));
          chaStruc->stats[c][idx].data = hist;

          readHist(b, hist);
        } else {
          // NCSU - Generic for non-histogram data
          uint16_t characteristic_size;
          characteristic_size = adios_get_stat_size(
              chaStruc->stats[c][idx].data, original_var_type, k);
          chaStruc->stats[c][idx].data = SAFE_MALLOC(characteristic_size);

          void *data = chaStruc->stats[c][idx].data;
          memcpy(data, (b->scratch + b->offset), characteristic_size);
          b->offset += characteristic_size;

          if (b->change_endianness == ADIOS_YES)
            // swap_ptr(data, characteristic_size * 8);
            not_supported("swap_ptr\n");
        }
        idx++;
      }
      k++;
    }
  }

  return 0;
}

void assignStat(struct adios_index_characteristic_struct_v1 *cc, void *data,
                enum ADIOS_STAT statType) {
  if (!(cc->stats)) {
    cc->stats =
        SAFE_MALLOC(sizeof(struct adios_index_characteristics_stat_struct *));
    cc->stats[0] =
        SAFE_MALLOC(2 * sizeof(struct adios_index_characteristics_stat_struct));
    cc->bitmap = 0;
  }
  cc->stats[0][statType].data = data;
  cc->bitmap |= (1 << statType);
}

void readStat(struct adios_index_characteristic_struct_v1 *cc,
              enum ADIOS_STAT statType, uint32_t size, H5_bp_file_t *b) {
  assignStat(cc, NULL, statType);

  cc->stats[0][statType].data = malloc(size);

  memcpy(cc->stats[0][statType].data, (b->scratch + b->offset), size);

  if (b->change_endianness == ADIOS_YES) {
    // swap_adios_type(var_header->characteristics.stats[0][statType].data,
    // var_header->type);
    not_supported("swap stat");
  }
  b->offset += size;
}

void handleNumbers(H5_bp_file_t *b, enum ADIOS_DATATYPES type) {
  uint16_t data_size;
  void *data = SAFE_MALLOC(data_size);

  memcpy(data, (b->scratch + b->offset), data_size);
  if (b->change_endianness == ADIOS_YES) {
    if (type == adios_complex) {
      // TODO
    } else if (type == adios_double_complex) {
      // TODO
    } else {
      switch (data_size) {
      case 2:
        // swap_16_ptr(data);
        not_supported("swap_16\n");
        break;
      case 4:
        // swap_32_ptr(data);
        not_supported("swap_32\n");
        break;
      case 8:
        // swap_64_ptr(data);
        not_supported("swap_64\n");
        break;
      case 16:
        // swap_128_ptr(data);
        not_supported("swap_128\n");
        break;
      }
    }
  }
  b->offset += data_size;
}

// this is just like handle_char_min_max
// except:
// (1) does not bother to backward compatible with min/max
// (2) handles nelemes > 1 for string arrays
//
int handle_char_value(H5_bp_file_t *b, struct adios_index_comp_struct_v1 **root,
                      enum ADIOS_CHARACTERISTICS c, uint64_t j) {
  uint16_t data_size;
  void *data = 0;

  if ((*root)->type == adios_string_array) {
    // there is no length info at this point
    data_size = 0;
    data = SAFE_MALLOC((*root)->nelems * sizeof(char *));
  } else if ((*root)->type == adios_string) {
    readInt2(b, &data_size);
    data = SAFE_MALLOC(data_size + 1);
    ((char *)data)[data_size] = '\0';
  } else {
    data_size = adios_get_type_size((*root)->type, "");
    data = SAFE_MALLOC((*root)->nelems * data_size);
  }

  switch ((*root)->type) {
  case adios_byte:
  case adios_short:
  case adios_integer:
  case adios_long:
  case adios_unsigned_byte:
  case adios_unsigned_short:
  case adios_unsigned_integer:
  case adios_unsigned_long:
  case adios_real:
  case adios_double:
  case adios_long_double:
  case adios_complex:
  case adios_double_complex:
    memcpy(data, (b->scratch + b->offset), (*root)->nelems * data_size);

    if (b->change_endianness == ADIOS_YES && (data_size > 1)) {
      if ((*root)->type == adios_complex) {
        // TODO
      } else if ((*root)->type == adios_double_complex) {
        // TODO
      } else {
        int k;
        char *p = (char *)data;
        for (k = 0; k < (*root)->nelems; k++) {
          switch (data_size) {
          case 2:
            // swap_16_ptr(p);
            not_supported("swap_16_ptr");
            break;
          case 4:
            // swap_32_ptr(p);
            not_supported("swap_32_ptr");
            break;
          case 8:
            // swap_64_ptr(p);
            not_supported("swap_64_ptr");
            break;
          case 16:
            // swap_128_ptr(p);
            not_supported("swap_128_ptr");
            break;
          }
          p += data_size;
        }
      }
    }

    b->offset += (*root)->nelems * data_size;
    break;

  case adios_string:
    memcpy(data, (b->scratch + b->offset), data_size);
    b->offset += data_size;
    break;

  case adios_string_array: {
    char **p = (char **)data;
    int k;
    for (k = 0; k < (*root)->nelems; k++) {
      data_size = *(uint16_t *)(b->scratch + b->offset);
      if (b->change_endianness == ADIOS_YES) {
        // swap_16(data_size);
        not_supported("swap_16");
      }
      b->offset += 2;
      p[k] = malloc(data_size + 1);
      p[k][data_size] = '\0';
      memcpy(p[k], (b->scratch + b->offset), data_size);
      b->offset += data_size;
    }
  } break;

  default:
    free(data);
    data = 0;
    break;
  }

  SAFE_FREE((*root)->characteristics[j].value);
  (*root)->characteristics[j].value = data;
  return 0;
}

int handle_char_min_max_vars(H5_bp_file_t *b,
                             struct adios_index_comp_struct_v1 **root,
                             enum ADIOS_CHARACTERISTICS c, uint64_t j) {
  uint16_t data_size;
  void *data = 0;

  if ((*root)->type == adios_string) {
    readInt2(b, &data_size);
  } else {
    data_size = adios_get_type_size((*root)->type, "");
  }

  switch ((*root)->type) {
  case adios_byte:
  case adios_short:
  case adios_integer:
  case adios_long:
  case adios_unsigned_byte:
  case adios_unsigned_short:
  case adios_unsigned_integer:
  case adios_unsigned_long:
  case adios_real:
  case adios_double:
  case adios_long_double:
  case adios_complex:
  case adios_double_complex:
    handleNumbers(b, (*root)->type);
    break;
  case adios_string:
    data = SAFE_MALLOC(data_size + 1);
    ((char *)data)[data_size] = '\0';
    memcpy(data, (b->scratch + b->offset), data_size);
    b->offset += data_size;
    break;

  default:
    data = 0;
    break;
  } // root type

  switch (c) {
  case ADIOS_CHARACTERISTIC_value:
    (*root)->characteristics[j].value = data;
    break;

  // NCSU - reading older bp files
  // ADIOS_CHARACTERISTIC_min, max are not used anymore. If this is encountered
  // it is an older bp file format Code below reads min and min, and sets the
  // bitmap for those 2 alone
  case ADIOS_CHARACTERISTIC_min:
    assignStat(&((*root)->characteristics[j]), data, ADIOS_STATISTIC_min);
    break;
  case ADIOS_CHARACTERISTIC_max:
    assignStat(&((*root)->characteristics[j]), data, ADIOS_STATISTIC_max);
    break;
  default:
    break;
  }
  // break;
  return 0;
}

int adios_parse_vars_index_v1(H5_bp_file_t *b,
                              struct adios_index_comp_struct_v1 **vars_root,
                              /*qhashtbl_t*/ void *hashtbl_vars,
                              struct adios_index_comp_struct_v1 **vars_tail) {
  struct adios_index_comp_struct_v1 **root;

  REQUIRE_SUCC_MSG((b->length - b->offset >= 10), 1,
                   "adios_parse_vars_index_v1 requires a buffer of at least 10 "
                   "bytes. Only %" PRId64 " were provided\n",
                   b->length - b->offset);
  root = vars_root;

  /* BP Format v1: vars_count and attrs_count was 16bit, now it's 32 bit */
  uint32_t vars_count;
  readInt4(b, &vars_count);
  uint64_t vars_length;
  readInt8(b, &vars_length);

  // validate remaining length
  uint32_t i;
  for (i = 0; i < vars_count; i++) {
    if (!*root) {
      *root = (struct adios_index_comp_struct_v1 *)SAFE_MALLOC(
          sizeof(struct adios_index_comp_struct_v1));
      (*root)->next = 0;
      (*root)->nelems = 1; // always one for nelem;
    }

    uint32_t var_entry_length;
    readInt4(b, &var_entry_length);
    /* BP Format v1: varid/attrid was 16bit, now it's 32 bit */
    readInt4(b, &((*root)->id));
    uint16_t len;
    readInt2(b, &len);

    readString(&((*root)->group_name), len, b);

    readInt2(b, &len);
    readString(&((*root)->name), len, b);

    readInt2(b, &len);
    readString(&((*root)->path), len, b);

    uint8_t flag;
    readInt1(b, &flag);
    (*root)->type = (enum ADIOS_DATATYPES)flag;

    uint64_t csets_count;
    readInt8(b, &csets_count);
    (*root)->characteristics_count = csets_count;
    (*root)->characteristics_allocated = csets_count;

    // validate remaining length: offsets_count * (8 + 2 * (size of type))
    (*root)->characteristics = SAFE_MALLOC(
        csets_count * sizeof(struct adios_index_characteristic_struct_v1));
    memset((*root)->characteristics, 0,
           csets_count * sizeof(struct adios_index_characteristic_struct_v1));

    uint64_t j;
    for (j = 0; j < csets_count; j++) {
      // NCSU - Clear stats structure (Drew: probably redundant with memset
      // above, but leave it to be safe)
      (*root)->characteristics[j].stats = 0;

      uint8_t characteristic_set_count;
      readInt1(b, &characteristic_set_count);
      uint32_t characteristic_set_length;
      readInt4(b, &characteristic_set_length);

      uint8_t item = 0;
      while (item < characteristic_set_count) {
        uint8_t flag;
        readInt1(b, &flag);
        enum ADIOS_CHARACTERISTICS c = (enum ADIOS_CHARACTERISTICS)flag;

        switch (c) {
        case ADIOS_CHARACTERISTIC_min:
        case ADIOS_CHARACTERISTIC_max:
        case ADIOS_CHARACTERISTIC_value: {
          // handle_char_min_max_vars(b, root, c,j);
          // NOTE: min max is not handled in handle_char_value() for simplicity,
          handle_char_value(b, root, c, j);
          break;
        }
        // NCSU - Statistics - Parsing stat related info from bp file based on
        // the bitmap
        case ADIOS_CHARACTERISTIC_stat: {
          // handle_char_stat(b,root,j);
          enum ADIOS_DATATYPES original_var_type =
              adios_transform_get_var_original_type_index(*root);
          handle_char_stat(b, original_var_type,
                           &((*root)->characteristics[j]));
          break;
        }
        // NCSU - Reading bitmap value
        case ADIOS_CHARACTERISTIC_bitmap: {
          readInt4(b, &((*root)->characteristics[j].bitmap));
          break;
        }
        case ADIOS_CHARACTERISTIC_offset: {
          readInt8(b, &((*root)->characteristics[j].offset));
          break;
        }
        case ADIOS_CHARACTERISTIC_payload_offset: {
          readInt8(b, &((*root)->characteristics[j].payload_offset));
          break;
        }
        case ADIOS_CHARACTERISTIC_file_index: {
          readInt4(b, &((*root)->characteristics[j].file_index));
          break;
        }
        case ADIOS_CHARACTERISTIC_time_index: {
          readInt4(b, &((*root)->characteristics[j].time_index));
          break;
        }
        case ADIOS_CHARACTERISTIC_dimensions: {
          readInt1(b, &((*root)->characteristics[j].dims.count));
          uint16_t dims_length;
          readInt2(b, &dims_length);

          (*root)->characteristics[j].dims.charDims =
              (uint64_t *)SAFE_MALLOC(dims_length);
          memcpy((*root)->characteristics[j].dims.charDims,
                 (b->scratch + b->offset), dims_length);
          if (b->change_endianness == ADIOS_YES) {
            uint16_t di = 0;
            uint16_t dims_num = dims_length / 8;
            for (di = 0; di < dims_num; di++) {
              // swap_64(((*root)->characteristics [j].dims.dims)[di]);
              not_supported("swap_64\n");
            }
          }
          b->offset += dims_length;
        } break;

        // NCSU ALACRITY-ADIOS - Reading variable transformation type
        case ADIOS_CHARACTERISTIC_transform_type: {
          // adios_transform_deserialize_transform_characteristic(&(*root)->characteristics[j].transform,
          // b);
          not_supported("adios_transform_deserialize_transform_characteristic");
          break;
        }
        case ADIOS_CHARACTERISTIC_var_id: {
          // this cannot happen, only attributes have variable references
          break;
        }
        }
        item++;
      }
    }

    // Add variable to the hash table too

    if (hashtbl_vars) {
      not_supported("hash table vars\n");
    }

    // record this as the tail
    if (vars_tail)
      *vars_tail = (*root);

    root = &(*root)->next;
  }
  return 0;
}

int adios_posix_read_process_group_index(H5_bp_file_t *b) {
  uint64_t pg_size = b->var_index - b->pg_index;
  realloc_aligned(b, pg_size);
  b->offset = 0;

  if (-1 == b->serial_handle) {
    b->serial_handle = open(b->file_name, O_RDONLY | O_LARGEFILE);
  }

  REQUIRE_SUCC_MSG((-1 != b->serial_handle), -1, "file cannt be open\n");

  lseek(b->serial_handle, b->pg_index, SEEK_SET);

  uint64_t r = read64(b->serial_handle, b->scratch, pg_size);
  REQUIRE_SUCC_MSG((r == pg_size), -1,
                   "reading group_index: wanted %" PRIu64 ", read: %" PRIu64
                   "\n",
                   pg_size, r);

  return 0;
}

int adios_parse_process_group_index_v1(
    H5_bp_file_t *b, struct adios_index_process_group_struct_v1 **pg_root,
    struct adios_index_process_group_struct_v1 **pg_tail) {
  struct adios_index_process_group_struct_v1 **root;

  REQUIRE_SUCC_MSG((b->length - b->offset >= 16), 1,
                   "adios_parse_process_group_index_v1 requires a buffer of at "
                   "least 16 bytes."
                   "Only %" PRId64 " were provided\n",
                   b->length - b->offset);

  root = pg_root;

  uint64_t process_groups_count;
  readInt8(b, &process_groups_count);
  uint64_t process_groups_length;
  readInt8(b, &process_groups_length);

  // validate remaining length
  uint64_t i;
  for (i = 0; i < process_groups_count; i++) {
    uint16_t length_of_group;
    readInt2(b, &length_of_group);

    if (!*root) {
      *root = (struct adios_index_process_group_struct_v1 *)SAFE_MALLOC(
          sizeof(struct adios_index_process_group_struct_v1));
      (*root)->is_time_aggregated = 0;
      (*root)->next = 0;
    }

    uint16_t length_of_name;
    readInt2(b, &length_of_name);
    //
    readString(&((*root)->group_name), length_of_name, b);

    char fortran_flag;
    readChar(b, &fortran_flag);
    (*root)->adios_host_language_fortran =
        (fortran_flag == 'y' ? ADIOS_YES : ADIOS_NO);

    readInt4(b, &((*root)->process_id));

    readInt2(b, &length_of_name);
    readString(&((*root)->time_index_name), length_of_name, b);

    readInt4(b, &((*root)->pg_time_index));
    readInt8(b, &((*root)->offset_in_file));

    // record this as the tail
    if (pg_tail)
      *pg_tail = (*root);

    root = &(*root)->next;
  }

  b->adios_host_language_fortran = (*pg_root)->adios_host_language_fortran;
  return 0;
}

void print_process_group_index(
    struct adios_index_process_group_struct_v1 *pg_root) {
  while (pg_root) {
    printf("Group: %s\n", pg_root->group_name);
    printf("\tProcess ID: %d\n", pg_root->process_id);
    printf("\tTime Name: %s\n", pg_root->time_index_name);
    printf("\tTime: %d\n", pg_root->pg_time_index);
    printf("\tOffset in File: %" PRIu64 "\n", pg_root->offset_in_file);

    pg_root = pg_root->next;
  }
}

void print_vars_index(struct adios_index_comp_struct_v1 *vars_root) {
  while (vars_root) {
    if (!strcmp(vars_root->path, "/")) {
      printf("Var (Group) [ID]: /%s (%s) [%d]\n", vars_root->name,
             vars_root->group_name, vars_root->id);
    } else {
      printf("Var (Group) [ID]:: %s/%s (%s) [%d]\n", vars_root->path,
             vars_root->name, vars_root->group_name, vars_root->id);
    }

    const char *typestr = print_adios_type(vars_root->type);
    printf("\tDatatype: %s\n", typestr);
    printf("\tVars Characteristics: %" PRIu64 "\n",
           vars_root->characteristics_count);

    int i;
    for (i = 0; i < vars_root->characteristics_count; i++) {
      struct adios_index_characteristic_struct_v1 *curr =
          &(vars_root->characteristics[i]);
      printf("\tOffset(%" PRIu64 ")", curr->offset);
      printf("\tPayload Offset(%" PRIu64 ")", curr->payload_offset);
      printf("\tFile Index(%d)", curr->file_index);
      printf("\tTime Index(%d)", curr->time_index);

      if (NULL != curr->value) {
        printf("\tValue(%s)", bp_value_to_string(vars_root->type, curr->value));
      }

      if (curr->stats && curr->stats[0][ADIOS_STATISTIC_min].data) {
        printf("\tMin(%s)",
               bp_value_to_string(vars_root->type,
                                  curr->stats[0][ADIOS_STATISTIC_min].data));
      }
      if (curr->stats && curr->stats[0][ADIOS_STATISTIC_max].data) {
        printf("\tMax(%s)",
               bp_value_to_string(vars_root->type,
                                  curr->stats[0][ADIOS_STATISTIC_max].data));
      }

      if (curr->dims.count != 0) {
        int j;
        printf("\tDims (l:g:o): (");

        for (j = 0; j < curr->dims.count; j++) {
          if (j != 0)
            printf(",");

          if (curr->dims.charDims[j * 3 + 1] != 0) {
            printf("%" PRIu64 ":%" PRIu64 ":%" PRIu64,
                   curr->dims.charDims[j * 3 + 0],
                   curr->dims.charDims[j * 3 + 1],
                   curr->dims.charDims[j * 3 + 2]);
          } else {
            printf("%" PRIu64, curr->dims.charDims[j * 3 + 0]);
          }
        }
        printf(")");
      }
      printf("\n");
    }

    vars_root = vars_root->next;
  }
}

/*
 */

int call_init(H5_bp_file_t *bpFile, MPI_File fh) {
  bpFile->scratch = NULL;
  bpFile->mpi_handle = fh;
  bpFile->serial_handle = -1;

  bpFile->pg_root = NULL;
  bpFile->vars_root = NULL;
  MPI_Offset file_size;
  int err = MPI_File_get_size(fh, &file_size);
  // H5Adios_list_MPI_Error(err);
  REQUIRE_MPI_SUCC(err);

  bpFile->file_size = file_size;
  bpFile->change_endianness = ADIOS_YES; // change it
  bpFile->offset = 0;
  return 0;
}

static void adios_clear_comp_index_v1(struct adios_index_comp_struct_v1 *root) {
  while (root) {
    int i;
    struct adios_index_comp_struct_v1 *temp = root->next;
    enum ADIOS_DATATYPES original_var_type =
        adios_transform_get_var_original_type_index(root);

    if (root->group_name)
      free(root->group_name);
    if (root->name)
      free(root->name);
    if (root->path)
      free(root->path);

    for (i = 0; i < root->characteristics_count; i++) {
      if (root->characteristics[i].dims.count != 0)
        free(root->characteristics[i].dims.charDims);
      if (root->characteristics[i].value)
        free(root->characteristics[i].value);

      // NCSU - Clears up the statistical data, based on bitmap
      if (root->characteristics[i].stats != 0) {
        uint8_t j = 0, idx = 0;
        uint8_t c = 0, count = adios_get_stat_set_count(original_var_type);

        for (c = 0; c < count; c++) {
          while (root->characteristics[i].bitmap >> j) {
            if ((root->characteristics[i].bitmap >> j) & 1) {
              if (j == ADIOS_STATISTIC_hist) {
                struct adios_index_characteristics_hist_struct *hist =
                    (struct adios_index_characteristics_hist_struct *)root
                        ->characteristics[i]
                        .stats[c][idx]
                        .data;
                free(hist->breaks);
                free(hist->frequencies);
              } else
                free(root->characteristics[i].stats[c][idx].data);
              idx++;
            }
            j++;
          }
          free(root->characteristics[i].stats[c]);
        }

        free(root->characteristics[i].stats);
      }
    }
    if (root->characteristics)
      free(root->characteristics);

    free(root);
    root = temp;
  }
}

void clean_scratch(H5_bp_file_t *bpFile) {
  if (bpFile->scratch != NULL) {
    free(bpFile->scratch);
  }
  bpFile->scratch = NULL;

  struct adios_index_process_group_struct_v1 *root = bpFile->pg_root;
  if (root != NULL) {
    while (root) {
      struct adios_index_process_group_struct_v1 *temp = root->next;
      if (root->group_name)
        free(root->group_name);
      if (root->time_index_name)
        free(root->time_index_name);
      free(root);
      root = temp;
    }
  }

  adios_clear_comp_index_v1(bpFile->vars_root);
  adios_clear_comp_index_v1(bpFile->attrs_root);
}

void close_BP_subfiles(H5_bp_file_t *bpFile) {
  BP_file_handle_list *lst = &bpFile->subfile_handles;

  struct BP_file_handle *l = lst->head;
  struct BP_file_handle *n;

  while (l) {
    n = l->next;

    MPI_File_close(&l->fh);
    free(l);

    l = n;
  }

  lst->n_handles = 0;
  lst->head = NULL;
  lst->tail = NULL;
}

void call_close(H5_bp_file_t *bpFile) {
  clean_scratch(bpFile);
  if (bpFile->serial_handle > -1) {
    close(bpFile->serial_handle);
  }
  close_BP_subfiles(bpFile);
}

int call_seek_back(H5_bp_file_t *bpFile, hsize_t start) {
  // printf("====todo: seek to pg_index and read the rest (iterate over
  // MAX_MPIWRITE if necessary) === \n");

  MPI_Status status;
  // bpFile->scratch = (char*) malloc(start+1);
  if (bpFile->scratch != NULL)
    free(bpFile->scratch);

  bpFile->scratch = (char *)SAFE_CALLOC(start + 1, sizeof(char));
  int err =
      MPI_File_seek(bpFile->mpi_handle, (MPI_Offset)(bpFile->file_size) - start,
                    MPI_SEEK_SET);

  err = MPI_File_read(bpFile->mpi_handle, bpFile->scratch, start, MPI_BYTE,
                      &status);

  // char buff[start+1];
  // MPI_File_read (bpFile->mpi_handle, buff, start, MPI_BYTE, &status);

  if (err)
    // H5Adios_list_MPI_Error(status.MPI_ERROR);
    REQUIRE_MPI_SUCC(status.MPI_ERROR);

  return 0;
}

void parse_characteristics(H5_bp_file_t *bpFile) {
  uint8_t flag;
  readInt1(bpFile, &flag);

  switch (flag) {
  case 1: // ADIOS_CHARACTERISTIC_min
  {
    double min = read_double(bpFile);
    break;
  }
  case 2: // ADIOS_CHARACTERISTIC_max
  {
    double max = read_double(bpFile);
    break;
  }
  case 3: // ADIOS_CHARACTERISTIC_offset
  {
    uint64_t val;
    readInt8(bpFile, &val);
    break;
  }
  case 4: // ADIOS_CHARACTERISTIC_dimensions     = 4
  {
    uint8_t dimCount;
    readInt1(bpFile, &dimCount);
    uint16_t dimLength;
    readInt2(bpFile, &dimLength);

    uint64_t currDimSize;
    uint16_t dims_num = dimLength / 8;
    uint16_t di;
    for (di = 0; di < dims_num; di++) {
      readInt8(bpFile, &currDimSize);
    }
    break;
  }
  case 6: // ADIOS_CHARACTERISTIC_payload_offset = 6
  {
    uint64_t val;
    readInt8(bpFile, &val);
    break;
  }
  case 7: // ADIOS_CHARACTERISTIC_file_index     = 7
  {
    uint32_t val;
    readInt4(bpFile, &val);
    break;
  }
  case 8: // ADIOS_CHARACTERISTIC_time_index     = 8
  {
    uint32_t val;
    readInt4(bpFile, &val);
    break;
  }
  case 0:  // ADIOS_CHARACTERISTIC_value
  case 5:  // ADIOS_CHARACTERISTIC_var_id        = 5
  case 9:  // ADIOS_CHARACTERISTIC_bitmap         = 9
  case 10: // ADIOS_CHARACTERISTIC_stat           = 10
  case 11: // ADIOS_CHARACTERISTIC_transform_type = 11
  default:
    printf(" Error: No SUCH support for charactereristic=%d\n", flag);
    break;
  }
}
/*
int call_adios_get_vars(H5_bp_file_t* bpFile)
{
  int bpversion = bpFile->bp_version & ADIOS_VERSION_NUM_MASK;

  REQUIRE_SUCC_MSG((bpFile->file_size - bpFile->offset >=
ADIOS_VARS_MINIHEADER_SIZE), 1, "insufficient bytes for reading var.");

  uint32_t var_count;
  if (bpversion > 1) {
    readInt4(bpFile, &var_count);
  } else {
    uint16_t bp1vc; readInt2(bpFile, &bp1vc);
    var_count = bp1vc;
  }

  uint64_t var_len; readInt8(bpFile, &var_len);
  // To speed find_var_byid(). Q. Liu, 11-2013.
  //fh->vars_table = (struct adios_index_var_struct_v1 **) malloc
(8*(size_t)mh->vars_count);
  // validate remaining length
  int i;
  uint64_t var_entry_length, vid, len, flag;
  uint64_t characteristics_sets_count;
  for (i = 0; i < var_count; i++) {
    // read one var:
    char* group_name=NULL; char* var_name=NULL; char* var_path=NULL;
    readInt(4, bpFile, &(bpFile->offset), &var_entry_length);
    if (bpversion > 1) {
      readInt(4, bpFile, &(bpFile->offset), &vid);
    } else {
      readInt(2, bpFile, &(bpFile->offset), &vid);
    }

    readInt(2, bpFile, &(bpFile->offset), &len);
    readString(&group_name, len, bpFile);

    readInt(2, bpFile, &(bpFile->offset), &len);
    readString(&var_name, len, bpFile);

    readInt(2, bpFile, &(bpFile->offset), &len);
    readString(&var_path, len, bpFile);

    readInt(1, bpFile, &(bpFile->offset), &flag);
    //(*root)->type = (enum ADIOS_DATATYPES) flag;

    readInt8(bpFile, &characteristics_sets_count);
    //readObj(characteristics, characteristics_sets_count, bpFile, offset);
    // validate remaining length: offsets_count *
    // (8 + 2 * (size of type))

    // NOTE: Above memset assumes that all 0's is a valid initialization.
    //       This is true, currently, but be careful in the future.

    uint64_t j;
    for (j = 0; j < characteristics_sets_count; j++)
    {
      uint8_t item = 0;
      uint8_t  cc; readInt1(bpFile, &cc);
      uint32_t ll; readInt4(bpFile, &ll);

      while (item < cc) {
        parse_characteristics (bpFile);
        item++;
      }
      //Old BP files do not have time_index characteristics, so we
        //       set it here automatically: j div # of pgs per timestep
          //     Assumed that in old BP files, all pgs write each variable in
each timestep.
        //
        //if ((*root)->characteristics [j].time_index == 0) {
        //  (*root)->characteristics [j].time_index =
        //    j / (mh->pgs_count / (fh->tidx_stop - fh->tidx_start + 1)) + 1;
    }
  }
  //process_joined_array((*root));
  //root = &(*root)->next;
  return 0;
}


void call_adios_get_pg(H5_bp_file_t* bpFile)
{
  //int offset=0;
  bpFile->offset = 0; // read from start of bpFile->scratch

  uint64_t pg_count; readInt8(bpFile, &pg_count);
  uint64_t skip; readInt8(bpFile, &skip);

  char ** namelist = (char **) SAFE_MALLOC(sizeof(char *) * pg_count);
  uint64_t * group_id_list = (uint64_t *) SAFE_MALLOC(sizeof(uint64_t) *
pg_count);

  int i,j;
  uint64_t group_count = 0;
  char fortran_flag;

  uint64_t tidx_start, tidx_stop; // Determine first and last timestep in file

  uint16_t strLen, skipme;
  for (i = 0; i < pg_count; i++) {
    namelist[i] = 0;

    readInt2(bpFile, &skipme); // length of group not used
    readInt2(bpFile, &strLen);

    char* group_name = NULL; readString(&group_name, strLen, bpFile);

    if ( group_count == 0 ) {
      namelist[group_count] = (char *) SAFE_MALLOC (strLen + 1);
      strcpy (namelist[group_count], group_name);
      ++group_count;
      group_id_list[i] = group_count-1;
    } else {
      for (j=0; j<group_count; j++) {
        if (!strcmp(namelist[j], group_name)) {
          break;
        }
      }
      if (j==group_count) {
        namelist[group_count] = (char *) SAFE_MALLOC (strLen + 1);
        strcpy (namelist[group_count], group_name);
        ++group_count;
        group_id_list[i] = group_count - 1;
      } else {
        group_id_list[i] = j;
      }
    }

    //uint64_t fortran_flag;
    readChar(bpFile, &fortran_flag);
    bpFile->adios_host_language_fortran = (fortran_flag == 'y' ? 1 : 0);

    uint32_t process_id; readInt4(bpFile, &process_id);

    readInt2(bpFile, &strLen);

    char* time_index_name = NULL; readString(&time_index_name, strLen, bpFile);

    uint32_t time_index; readInt4(bpFile, &time_index);
    uint64_t offset_in_file; readInt8(bpFile, &offset_in_file);

    if (i == 0)
      tidx_start = time_index;
    if (i == pg_count-1) {
      tidx_stop =  time_index;
    }

    // now read the next group
    printf("Now ready to read next group\n");
  }
  bpFile->time_steps = tidx_stop - tidx_start + 1;
}
*/

int call_seek(H5_bp_file_t *bpFile, hsize_t start) {
  MPI_Status status;
  hsize_t toRead = (bpFile->file_size) - start;
  if (bpFile->scratch != NULL) {
    free(bpFile->scratch);
    bpFile->scratch = 0;
  }

  bpFile->scratch = (char *)SAFE_CALLOC(toRead + 1, sizeof(char));
  int err = MPI_File_seek(bpFile->mpi_handle, (MPI_Offset)start, MPI_SEEK_SET);
  err = MPI_File_read(bpFile->mpi_handle, bpFile->scratch, toRead, MPI_BYTE,
                      &status);

  if (err)
    // H5Adios_list_MPI_Error(status.MPI_ERROR);
    REQUIRE_MPI_SUCC(status.MPI_ERROR);

  //
  // scratch contents will be read later from realloc_aligned()
  // and no malloc is needed. so free here.
  free(bpFile->scratch);

  REQUIRE_SUCC((0 == adios_posix_read_process_group_index(bpFile)), -1);
  REQUIRE_SUCC((0 == adios_parse_process_group_index_v1(
                         bpFile, &(bpFile->pg_root), NULL)),
               -1);

  REQUIRE_SUCC((0 == adios_posix_read_vars_index(bpFile)), -1);
  REQUIRE_SUCC((0 == adios_parse_vars_index_v1(bpFile, &(bpFile->vars_root),
                                               NULL, NULL)),
               -1);

  // print_process_group_index(bpFile->pg_root);
  // print_vars_index(bpFile->vars_root);

  // call_adios_get_pg(bpFile);
  // call_adios_get_vars(bpFile);

  REQUIRE_SUCC((0 == adios_posix_read_attributes_index(bpFile)), -1);

  REQUIRE_SUCC(
      (0 == adios_parse_attributes_index_v1(bpFile, &(bpFile->attrs_root))),
      -1);

  // print_attributes_index (bpFile->attrs_root);

  return 0;
}

int call_adios_has_var(H5_bp_file_t *bpFile, const char *name) {
  REQUIRE_NOT_NULL(bpFile);
  REQUIRE_NOT_NULL(name);

  struct adios_index_comp_struct_v1 *vars_root = bpFile->vars_root;

  while (vars_root) {
    if (!strcmp(vars_root->name, name)) {
      return 1;
    }
    vars_root = vars_root->next;
  }

  return 0;
}

int call_adios_get_version(/*char* buff,*/ H5_bp_file_t *bpFile, int offset) {
  // char buff[ADIOS_BP_MINIFOOTER_SIZE+1];
  char *buff = bpFile->scratch;
  uint32_t test = 1; // high bit big indian
  bpFile->bp_version = ntohl(*(uint32_t *)(buff + offset));
  char *v = (char *)&(bpFile->bp_version);
  if ((*v && !*(char *)&test) // both writer and this machine are big endian
      || (!*(v + 3) && *(char *)&test) // both are little endian
  ) {
    bpFile->change_endianness = ADIOS_NO;
  }

  bpFile->bp_version = bpFile->bp_version & 0x7fffffff;

  bpFile->pg_index = *(uint64_t *)(buff);
  // printf("pg_index: %lld\n", bpFile->pg_index);

  // assert(pg_index+28 < file_size)
  bpFile->var_index = *(uint64_t *)(buff + 8);
  // printf("var_index: %lld\n", bpFile->var_index);

  // assert(var_index > pg_dex)
  // assert(var_index+28 < file_size)

  bpFile->attr_index = *(uint64_t *)(buff + 16);
  // printf("attr_index: %lld\n", bpFile->attr_index);
  // assert(attr_index+28 < file_size)
  // assert(attr_index > var_index)

  // uint64_t attrs_end = b->file_size - MINIFOOTER_SIZE;
  bpFile->attrs_end = bpFile->file_size - ADIOS_BP_MINIFOOTER_SIZE;

  uint64_t pg_size = bpFile->var_index - bpFile->pg_index;
  uint64_t var_size = bpFile->attr_index - bpFile->var_index;
  uint64_t attr_size = bpFile->attrs_end - bpFile->attr_index;

  // printf("Checking: pg_size=%lld, var_size=%lld, attr_size=%lld\n", pg_size,
  // var_size, attr_size);

  return 0;
}

void read_scalar(H5_bp_request_t *req, int64_t start_idx,
                 enum ADIOS_YES_NOFLAG has_subfile, H5_bp_file_t *fh) {
  struct adios_index_comp_struct_v1 *v = req->var;
  int size_of_type =
      adios_get_type_size(v->type, v->characteristics[start_idx].value);
  void *data = req->output_data;
  /* READ A SCALAR VARIABLE */
  /* Prepare slice_offset, slice_size and idx for the later macro:
     MPI_FILE_READ_OPS1 and MPI_FILE_READ_OPS2
  */
  int idx = 0;
  uint64_t slice_offset = v->characteristics[start_idx + idx].payload_offset;
  uint64_t slice_size = size_of_type;

  if (v->type == adios_string) {
    // Note: strings are stored without \0 in file
    // size_of_type above includes \0 so decrease by one
    // size_of_type--;
    // JG: I comment this out b/c it does not seem \0 is at the end when I
    // tested scalar.bp generated from adios1
  }

  if (v->characteristics[start_idx + idx].payload_offset > 0) {
    if (has_subfile == ADIOS_NO) {
      // MPI_FILE_READ_OPS1;
      mpi_OPS1_readvalue(slice_size, slice_offset, fh);
    } else {
      if (mpi_OPS2_readvalue(fh, v, slice_size, slice_offset,
                             start_idx + idx) == 0)
        return;
    }
  } else {
    slice_offset = 0;
    mpi_OPS3_readvalue(start_idx + idx, fh, v);
  }

  memcpy((char *)data, fh->scratch + fh->offset, size_of_type);

  if (fh->change_endianness == ADIOS_YES) {
    change_endianness((char *)data, size_of_type, v->type);
  }

  if (v->type == adios_string) {
    // add \0 to the end of string
    // size_of_type here is the length of string
    // FIXME: how would this work for strings written over time?
    ((char *)data)[size_of_type] = '\0';
  }

  data = (char *)data +
         (v->type == adios_string ? size_of_type + 1 : size_of_type);
}

void read_pg(uint64_t payload_size, int64_t currIdx,
             enum ADIOS_YES_NOFLAG has_subfile, H5_bp_file_t *fp,
             H5_bp_request_t *req) {
  uint64_t slice_offset = 0, slice_size = 0;
  void *data = req->output_data;

  /* The complete read happens to be exactly one pg, and the entire pg */
  /* This means we enter this only once, and npg=1 at the end */
  /* This is a rare case. FIXME: cannot eliminate this? */

  slice_size = payload_size;
  slice_offset = req->var->characteristics[currIdx].payload_offset;
  if (req->var->characteristics[currIdx].payload_offset > 0) {
    if (has_subfile == ADIOS_NO) {
      mpi_OPS1_readvalue(slice_size, slice_offset, fp);
    } else {
      if (mpi_OPS2_readvalue(fp, req->var, slice_size, slice_offset, currIdx) ==
          0)
        return;
    }
  } else {
    slice_offset = 0;
    mpi_OPS3_readvalue(currIdx, fp, req->var);
  }

  memcpy((char *)data, fp->scratch + fp->offset, slice_size);

  if (fp->change_endianness == ADIOS_YES) {
    change_endianness(data, slice_size, req->var->type);
  }
}

void read_block(uint64_t datasize, const uint64_t *ldims,
                const uint64_t *offsets, int size_of_type, int currIdx,
                enum ADIOS_YES_NOFLAG has_subfile, H5_bp_file_t *fp,
                H5_bp_request_t *req) {
  // JG: so, only the slowest dimension needs to be truncated
  // the datasize is the size of dims1->ndim.
  // so just read from offset[0]*datasize*typesize
  //
  /* The slowest changing dimensions should not be read completely but
     we still need to read only one block */

  int ndim = req->sel->ndim;
  uint64_t *count = req->sel->count;
  uint64_t *start = req->sel->start;

  uint64_t isize;
  uint64_t size_in_dset = 0, offset_in_dset = 0, offset_in_var = 0;

  isize = offsets[0] + ldims[0];
  if (start[0] >= offsets[0]) {
    // head is in
    if (start[0] < isize) {
      if (start[0] + count[0] > isize)
        size_in_dset = isize - start[0];
      else
        size_in_dset = count[0];
      offset_in_dset = start[0] - offsets[0];
      offset_in_var = 0;
    }
  } else {
    // middle is in
    if (isize < start[0] + count[0])
      size_in_dset = ldims[0];
    else
      // tail is in
      size_in_dset = count[0] + start[0] - offsets[0];

    offset_in_dset = 0;
    offset_in_var = offsets[0] - start[0];
  }

  uint64_t slice_size = size_in_dset * datasize * size_of_type;
  uint64_t write_offset = offset_in_var * datasize * size_of_type;
  uint64_t slice_offset = 0;

  if (req->var->characteristics[currIdx].payload_offset > 0) {
    slice_offset = req->var->characteristics[currIdx].payload_offset +
                   offset_in_dset * datasize * size_of_type;
    if (has_subfile == ADIOS_NO) {
      mpi_OPS1_readvalue(slice_size, slice_offset, fp);
    } else {
      if (mpi_OPS2_readvalue(fp, req->var, slice_size, slice_offset, currIdx) ==
          0)
        return;
    }
  } else {
    slice_offset = 0;
    mpi_OPS3_readvalue(currIdx, fp, req->var);
  }

  void *data = req->output_data;
  memcpy((char *)data + write_offset, fp->scratch + fp->offset, slice_size);
  if (fp->change_endianness == ADIOS_YES) {
    not_supported("change endianness");
    // change_endianness((char *)data + write_offset, slice_size, v->type);
  }

  // write_offset +=  slice_size;
}

void adios_util_copy_data(void *dst, void *src, int idim, int ndim,
                          uint64_t *size_in_dset, uint64_t *ldims,
                          const uint64_t *readsize, uint64_t dst_stride,
                          uint64_t src_stride, uint64_t dst_offset,
                          uint64_t src_offset, uint64_t ele_num,
                          int size_of_type,
                          enum ADIOS_YES_NOFLAG change_endianess,
                          enum ADIOS_DATATYPES type)

{
  unsigned int i, j;

  if (ndim - 1 == idim) {
    for (i = 0; i < size_in_dset[idim]; i++) {
      memcpy((char *)dst + (i * dst_stride + dst_offset) * size_of_type,
             (char *)src + (i * src_stride + src_offset) * size_of_type,
             ele_num * size_of_type);
      if (change_endianess == ADIOS_YES) {
        not_supported("change_endianess\n");
        // change_endianness ((char *)dst +
        // (i*dst_stride+dst_offset)*size_of_type, ele_num*size_of_type, type);
      }
    }
    return;
  }

  uint64_t dst_offset_new = 0;
  uint64_t src_offset_new = 0;
  uint64_t src_step, dst_step;

  for (i = 0; i < size_in_dset[idim]; i++) {
    // get the different step granularity
    // for each different reading pattern broke
    src_step = 1;
    dst_step = 1;
    for (j = idim + 1; j <= ndim - 1; j++) {
      src_step *= ldims[j];
      dst_step *= readsize[j];
    }
    src_offset_new = src_offset + i * src_stride * src_step;
    dst_offset_new = dst_offset + i * dst_stride * dst_step;
    adios_util_copy_data(dst, src, idim + 1, ndim, size_in_dset, ldims,
                         readsize, dst_stride, src_stride, dst_offset_new,
                         src_offset_new, ele_num, size_of_type,
                         change_endianess, type);
  }
}

/*
void adios_util_copy_data (void *dst, void *src,
                           int idim,
                           int ndim,
                           uint64_t* size_in_dset,
                           uint64_t* ldims,
                           const uint64_t * readsize,
                           uint64_t dst_stride,
                           uint64_t src_stride,
                           uint64_t dst_offset,
                           uint64_t src_offset,
                           uint64_t ele_num,
                           int      size_of_type,
                           enum ADIOS_YES_NOFLAG change_endianess,
                           enum ADIOS_DATATYPES type)

{
  unsigned int i, j;

  if (ndim-1==idim) {
      for (i=0;i<size_in_dset[idim];i++) {
           memcpy ((char *)dst + (i*dst_stride+dst_offset)*size_of_type,
                   (char *)src + (i*src_stride+src_offset)*size_of_type,
                   ele_num*size_of_type);
           if (change_endianess == ADIOS_YES) {
               not_supported("change_endianess\n");
               //change_endianness ((char *)dst +
(i*dst_stride+dst_offset)*size_of_type, ele_num*size_of_type, type);
           }
      }
      return;
  }

  uint64_t dst_offset_new=0;
  uint64_t src_offset_new=0;
  uint64_t src_step, dst_step;

  for (i = 0; i<size_in_dset[idim];i++) {
    // get the different step granularity
    // for each different reading pattern broke
    src_step = 1;
    dst_step = 1;
    for (j = idim+1; j <= ndim-1;j++) {
      src_step *= ldims[j];
      dst_step *= readsize[j];
    }
    src_offset_new =src_offset + i * src_stride * src_step;
    dst_offset_new = dst_offset + i * dst_stride * dst_step;
    adios_util_copy_data ( dst, src, idim+1, ndim, size_in_dset,
                           ldims,readsize,
                           dst_stride, src_stride,
                           dst_offset_new, src_offset_new,
                           ele_num, size_of_type, change_endianess, type);
  }
}
*/
void read_partial(int hole_break, uint64_t *ldims, uint64_t *offsets,
                  int size_of_type, int64_t currIdx,
                  enum ADIOS_YES_NOFLAG has_subfile, H5_bp_file_t *fp,
                  H5_bp_request_t *req)

{
  int ndim = req->sel->ndim;
  uint64_t *count = req->sel->count;
  uint64_t *start = req->sel->start;

  uint64_t isize, datasize, var_stride, dset_stride, slice_size, slice_offset;
  uint64_t size_in_dset[10], offset_in_dset[10], offset_in_var[10];

  adios_util_count_partial_pg(ndim, size_in_dset, offset_in_dset, offset_in_var,
                              start, count, ldims, offsets);

  datasize = 1;
  var_stride = 1;
  dset_stride = 1;

  int i;
  for (i = ndim - 1; i >= hole_break; i--) {
    datasize *= size_in_dset[i];
    dset_stride *= ldims[i];
    var_stride *= count[i];
  }

  uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
  for (i = ndim - 1; i > -1; i--) {
    start_in_payload += s * offset_in_dset[i] * size_of_type;
    end_in_payload +=
        s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
    s *= ldims[i];
  }

  slice_size = end_in_payload - start_in_payload + 1 * size_of_type;
  if (req->var->characteristics[currIdx].payload_offset > 0) {
    slice_offset =
        req->var->characteristics[currIdx].payload_offset + start_in_payload;
    if (has_subfile == ADIOS_NO) {
      mpi_OPS1_readvalue(slice_size, slice_offset, fp);
    } else {
      if (mpi_OPS2_readvalue(fp, req->var, slice_size, slice_offset, currIdx) ==
          0)
        return;
    }
  } else {
    slice_offset = start_in_payload;
    mpi_OPS3_readvalue(currIdx, fp, req->var);
  }

  for (i = 0; i < ndim; i++) {
    offset_in_dset[i] = 0;
  }

  uint64_t var_offset = 0;
  uint64_t dset_offset = 0;

  for (i = 0; i < ndim; i++) {
    var_offset = offset_in_var[i] + var_offset * count[i];
    dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
  }

  void *data = req->output_data;

  adios_util_copy_data(data, fp->scratch + fp->offset, 0, hole_break,
                       size_in_dset, ldims, count, var_stride, dset_stride,
                       var_offset, dset_offset, datasize, size_of_type,
                       fp->change_endianness, req->var->type);
  return;
}

void read_array(int64_t start_idx, int64_t stop_idx,
                enum ADIOS_YES_NOFLAG has_subfile, H5_bp_file_t *fp,
                H5_bp_request_t *req, int ndim) {
  uint64_t *start = req->sel->start;
  uint64_t *count = req->sel->count;

  int64_t idx = 0;
  // int ndim = req->var->characteristics[0].dims.count; // this may variable
  // according to versions, see adios_get_var_ndim()
  int i, j;
  uint64_t ldims[ndim], gdims[ndim], offsets[ndim];
  int size_of_type =
      adios_get_type_size(req->var->type, req->var->characteristics[0].value);
  /* READ AN ARRAY VARIABLE */
  uint64_t write_offset = 0;

  // loop over the list of pgs to read from one-by-one
  for (idx = 0; idx < stop_idx - start_idx + 1; idx++) {
    // int flag;
    uint64_t datasize = 1;
    // var_stride = 1; dset_stride = 1;
    // int idx_useful=1;
    int currIdx = start_idx + idx;
    int is_global = adios_get_var_dimensions(
        &(req->var->characteristics[currIdx].dims), ndim, ldims, gdims, offsets,
        fp->adios_host_language_fortran);

    if (!is_global) {
      // we use gdims below, which is 0 for a local array; set to ldims here
      for (j = 0; j < ndim; j++) {
        gdims[j] = ldims[j];
      }
      // we need to read only the first PG, not all, so let's prevent a second
      // loop
      stop_idx = start_idx;
    }

    uint64_t payload_size = size_of_type;
    int idx_useful = adios_util_check_selection(ndim, req, gdims, ldims,
                                                offsets, &payload_size);

    if (!idx_useful)
      continue;

    /* determined how many (fastest changing) dimensions can we read in in one
     * read */
    for (i = ndim - 1; i > -1; i--) {
      if (offsets[i] == start[i] && ldims[i] == count[i]) {
        datasize *= ldims[i];
      } else
        break;
    }

    int hole_break = i;

    if (hole_break == -1) {
      read_pg(payload_size, start_idx + idx, has_subfile, fp, req);
    } else if (hole_break == 0) {
      read_block(datasize, ldims, offsets, size_of_type, start_idx + idx,
                 has_subfile, fp, req);
    } else {
      read_partial(hole_break, ldims, offsets, size_of_type, start_idx + idx,
                   has_subfile, fp, req);
    }
  } // end for (idx ... loop over pgs
}

// adios: bp_util.c  has_subfiles
enum ADIOS_YES_NOFLAG checkSubfile(H5_bp_file_t *fp) {
  uint32_t result = (fp->bp_version & ADIOS_VERSION_HAVE_SUBFILE);
  if (result == 0) {
    return ADIOS_NO;
  } else {
    return ADIOS_YES;
  }
}

//
// read_request: provides:
// bbsel;
//
int read_var_bb(H5_bp_file_t *fp, H5_bp_request_t *req)
// const ADIOS_FILE *fp, read_request * r)
{
  // BP_PROC * p = GET_BP_PROC (fp);
  // BP_FILE * fh = GET_BP_FILE (fp);

  void *data = req->output_data;
  // uint64_t* start = req->sel->start;
  uint64_t *count = req->sel->count;

  struct adios_index_comp_struct_v1 *v = req->var;

  enum ADIOS_YES_NOFLAG file_is_fortran = fp->adios_host_language_fortran;
  enum ADIOS_YES_NOFLAG has_subfile = checkSubfile(fp);

  int i, j, t, time;
// uint64_t* tmpcount;
/* Get dimensions and flip if caller != writer language */
/* Note: ndim below doesn't include time if there is any */
// NCSU ALACRITY-ADIOS - Note: this function has been modified to return
//   the "raw" dimensions (i.e., 1D byte array)
#ifdef NEVER
  int ndim, nsteps;
  bp_get_and_swap_dimensions(fp, v, file_is_fortran, &ndim, &dims, &nsteps,
                             file_is_fortran);
#else // get ndim and nsteps directly
  // int ndim = req->var->characteristics[0].dims.count;
  int idx = adios_get_ch_idx(req->var, req->from_step);
  int ndim = adios_get_var_ndim(fp, idx, req->var, NULL);
  int nsteps = adios_get_var_nsteps(req->var);
#endif

  REQUIRE_SUCC_MSG((ndim == req->sel->ndim), 0,
                   "Error: variable ndim mismatch");

  // JG:no fortran call here so skip fortran order swap

  uint64_t datasize, dset_stride, var_stride, total_size = 0, items_read = 1;

  int dummy = -1, is_global = 0;

  MPI_Status status;
  H5_adios_varchunk_t *chunk;
  // not used struct adios_var_header_struct_v1 var_header;
  // struct adios_var_payload_struct_v1 var_payload;
  int size_of_type = adios_get_type_size(v->type, v->characteristics[0].value);

  for (i = 0; i < ndim; i++) {
    items_read *= count[i];
  }

  /* Note fp->current_step is always 0 for file mode. */
  int startTimeStep = fp->current_step + req->from_step;
  for (t = startTimeStep; t < startTimeStep + req->nsteps; t++) {
    if (!fp->streaming) {
      time = adios_get_time(v, t);
    } else {
#ifdef NEVER
      // Fix: the assumption that for streaming mode, the time in file
      // always starts from 1 is not correct. So here we add fh->tidx_start to
      // adjust Q. Liu, 06/2013
      time = fp->tidx_start + t;
#else
      SHOW_ERROR_MSG(
          "Don't think ADIOS streaming mode will be specified from HDF5.");
#endif
    }

    int64_t start_idx = adios_get_var_start_index(v, time);
    int64_t stop_idx = adios_get_var_stop_index(v, time);

    if (start_idx < 0 || stop_idx < 0) {
      SHOW_ERROR_MSG("Variable %s has no data at %d time step\n", v->name, t);
      continue;
    }

    if (ndim == 0) {
      read_scalar(req, start_idx, has_subfile, fp);
    } else {
      read_array(start_idx, stop_idx, has_subfile, fp, req, ndim);
      total_size += items_read * size_of_type;
      // shift target pointer for next read in
      data = (char *)data + (items_read * size_of_type);
    }

  } // end timesteps "t" loop

  // free (dims);

  return 1;
}
