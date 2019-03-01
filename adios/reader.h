#ifndef ADIOS_READER_H
#define ADIOS_READER_H

#include "error.h"
#include "hdf5.h"
#include "mpi.h"

#define ADIOS_BP_MINIFOOTER_SIZE 28
#define ADIOS_VERSION_NUM_MASK 0x000000FF
#define ADIOS_VERSION_HAVE_SUBFILE 0x00000100
#define ADIOS_VARS_MINIHEADER_SIZE 10
#define ADIOS_VERSION_HAVE_TIME_INDEX_CHARACTERISTIC 0x00000200

#define ADIOS_MAX_HANDLES 512

#ifndef ADIOS_STAT_LENGTH
#define ADIOS_STAT_LENGTH 7
#endif

#define MAX_MPIWRITE_SIZE 2130706432 // from adios

// NCSU - Adding statistics
enum ADIOS_STAT {
  ADIOS_STATISTIC_min = 0,
  ADIOS_STATISTIC_max = 1,
  ADIOS_STATISTIC_cnt = 2,
  ADIOS_STATISTIC_sum = 3,
  ADIOS_STATISTIC_sum_square = 4,
  ADIOS_STATISTIC_hist = 5,
  ADIOS_STATISTIC_finite = 6
};

typedef struct {
  int ndim;
  uint64_t *start;
  uint64_t *count;

} H5_adios_bb_sel;

struct adios_hist_struct {
  double min; // minimum value of histogram ** for when we use complex variables
  double max; // maximum value of histogram
  uint32_t num_breaks;   // number of break points for the histogram
  uint32_t *frequencies; // array of frequencies for the histogram
  double *breaks; // breaks array for the histogram, output this to gnuplot
};

enum ADIOS_YES_NOFLAG { ADIOS_UNKNOWN = 0, ADIOS_YES = 1, ADIOS_NO = 2 };

enum ADIOS_DATATYPES {
  adios_unknown = -1 /* (size) */
  ,
  adios_byte = 0 /* (1) */
  ,
  adios_short = 1 /* (2) */
  ,
  adios_integer = 2 /* (4) */
  ,
  adios_long = 4 /* (8) */

  ,
  adios_unsigned_byte = 50 /* (1) */
  ,
  adios_unsigned_short = 51 /* (2) */
  ,
  adios_unsigned_integer = 52 /* (4) */
  ,
  adios_unsigned_long = 54 /* (8) */

  ,
  adios_real = 5 /* (4) */
  ,
  adios_double = 6 /* (8) */
  ,
  adios_long_double = 7 /* (16) */

  ,
  adios_string = 9 /* (?) */
  ,
  adios_complex = 10 /* (8) */
  ,
  adios_double_complex = 11 /* (16) */

  /* Only for attributes: char** array of strings.
     Number of elements must be known externally */
  ,
  adios_string_array = 12 /* (sizeof(char*)) usually 4 */
};

enum ADIOS_CHARACTERISTICS {
  ADIOS_CHARACTERISTIC_value = 0,
  ADIOS_CHARACTERISTIC_min =
      1 // This is no longer used. Used to read in older bp file format
  ,
  ADIOS_CHARACTERISTIC_max =
      2 // This is no longer used. Used to read in older bp file format
  ,
  ADIOS_CHARACTERISTIC_offset = 3,
  ADIOS_CHARACTERISTIC_dimensions = 4,
  ADIOS_CHARACTERISTIC_var_id = 5,
  ADIOS_CHARACTERISTIC_payload_offset = 6,
  ADIOS_CHARACTERISTIC_file_index = 7,
  ADIOS_CHARACTERISTIC_time_index = 8,
  ADIOS_CHARACTERISTIC_bitmap = 9,
  ADIOS_CHARACTERISTIC_stat = 10,
  ADIOS_CHARACTERISTIC_transform_type = 11
};

/*
typedef struct
{
  uint32_t count;
  uint64_t length;
} H5_bp_vars_header_t; // i.e. struct adios_vars_header_struct_v1; and it seems
not actively used
*/

struct adios_dimension_item_struct_v1 {
  uint64_t rank;
  uint32_t var_id;
  enum ADIOS_YES_NOFLAG is_time_index;
};

struct adios_dimension_struct_v1 {
  struct adios_dimension_item_struct_v1 dimension;
  struct adios_dimension_item_struct_v1 global_dimension;
  struct adios_dimension_item_struct_v1 local_offset;
  struct adios_dimension_struct_v1 *next;
};

struct BP_file_handle {
  uint32_t file_index;
  MPI_File fh;
  struct BP_file_handle *next;
  struct BP_file_handle *prev;
};

struct BP_file_handle_header {
  uint32_t n_handles; // number of handles open
  struct BP_file_handle *head;
  struct BP_file_handle *tail;
  int warning_printed; // 1 if already printed the warning
};

typedef struct BP_file_handle_header BP_file_handle_list;

typedef struct {
  MPI_File mpi_handle;
  int serial_handle;
  hsize_t file_size;
  struct BP_file_handle_header subfile_handles;

  //
  // index indicates the location
  // written in the bpfile
  //
  hsize_t pg_index;
  hsize_t var_index;
  hsize_t attr_index;

  hsize_t attrs_end;

  uint32_t bp_version;

  char *scratch;
  uint64_t length;
  uint64_t offset;          // for scratch
  char *allocated_buff_ptr; // initial alloc for aligning on 8-byte boundary

  int time_steps;

  enum ADIOS_YES_NOFLAG change_endianness;
  enum ADIOS_YES_NOFLAG adios_host_language_fortran;

  char *file_name;

  struct adios_index_process_group_struct_v1 *pg_root;
  struct adios_index_comp_struct_v1 *vars_root;
  struct adios_index_comp_struct_v1 *attrs_root;

  uint64_t read_pg_offset;
  uint64_t read_pg_size;

  int current_step;
  // from BP_PROC
  ////BP_FILE * fh;
  int streaming;
  int *varid_mapping;
  // read_request * local_read_request_list;
  // void * b; //internal buffer for chunk reading
  // void * priv;

} H5_bp_file_t;
// H5_bp_file_t is a mix of original
// adios_bp_buffer_struct_v1, and others, ..
//

struct adios_index_comp_struct_v1 *adios_get_var_byname(H5_bp_file_t *file,
                                                        const char *varname);
struct adios_index_comp_struct_v1 *adios_get_attr_byname(H5_bp_file_t *file,
                                                         const char *attrname);

int adios_get_var_nsteps(struct adios_index_comp_struct_v1 *var_root);
int adios_get_ch_idx(struct adios_index_comp_struct_v1 *var, int curr_timestep);
int adios_get_time(struct adios_index_comp_struct_v1 *v, int step);
int64_t adios_get_var_start_index(struct adios_index_comp_struct_v1 *v, int t);
int64_t adios_get_var_stop_index(struct adios_index_comp_struct_v1 *v, int t);

typedef struct {
  // int ndim;
  // uint64_t* start;
  // uint64_t* count;
  H5_adios_bb_sel *sel;
  int from_step;
  int nsteps;
  void *output_data;

  struct adios_index_comp_struct_v1 *var;
} H5_bp_request_t;

struct adios_index_process_group_struct_v1 {
  char *group_name;
  enum ADIOS_YES_NOFLAG adios_host_language_fortran;
  uint32_t process_id;
  char *time_index_name;
  uint32_t pg_time_index;
  uint64_t offset_in_file;
  int is_time_aggregated;

  struct adios_index_process_group_struct_v1 *next;
};

//
// comp = var or attr
// only difference is attr can take string array
//
struct adios_index_comp_struct_v1 {
  uint32_t id;
  char *group_name;
  char *name;
  char *path;
  enum ADIOS_DATATYPES type;

  uint64_t characteristics_count;
  uint64_t characteristics_allocated;

  struct adios_index_characteristic_struct_v1 *characteristics;

  struct adios_index_comp_struct_v1 *next;

  int nelems; // for array of values for attributes
};

struct adios_index_characteristic_dims_struct_v1 {
  uint8_t count;
  uint64_t *charDims; // each 3 uint64_t represents one dimension (l, g, o)
};

struct adios_index_characteristic_transform_struct {
  uint8_t transform_type;

  enum ADIOS_DATATYPES pre_transform_type;
  struct adios_index_characteristic_dims_struct_v1 pre_transform_dimensions;

  uint16_t transform_metadata_len;
  void *transform_metadata;
};

struct adios_index_characteristics_stat_struct {
  void *data;
};

// NCSU - Structure for histogram
struct adios_index_characteristics_hist_struct {
  double min; // minimum value of histogram ** for when we use complex variables
  double max; // maximum value of histogram
  uint32_t num_breaks;   // number of break points for the histogram
  uint32_t *frequencies; // array of frequencies for the histogram
  double *breaks; // breaks array for the histogram, output this to gnuplot
};

struct adios_index_characteristic_struct_v1 {
  uint64_t offset; // beginning of the var or attr entry
  struct adios_index_characteristic_dims_struct_v1 dims;
  uint32_t var_id;
  void *value;
  uint64_t payload_offset; // beginning of the var or attr payload
  uint32_t file_index;     // subfile index
  uint32_t time_index;

  uint32_t bitmap;

  struct adios_index_characteristics_stat_struct **stats;

  // NCSU ALACRITY-ADIOS - Adding transform-related fields
  /*
    uint8_t transform_type;
    enum ADIOS_DATATYPES pre_transform_type;
    struct adios_index_characteristic_dims_struct_v1 pre_transform_dimensions;
    uint16_t transform_metadata_len;
    void *transform_metadata;
  */
  struct adios_index_characteristic_transform_struct transform;
};

typedef struct {
  uint32_t id;
  char *name;
  char *path;
  enum ADIOS_DATATYPES type;
  enum ADIOS_YES_NOFLAG is_dim;

  struct adios_dimension_struct_v1 *dims;
  struct adios_index_characteristic_struct_v1 characteristics;
  uint64_t payload_size;
} H5_bp_var_data_header_t; // i.e. struct adios_var(no s)_eader_struct_v1;;

typedef struct {
  int varid; /* variable index (0..ADIOS_FILE.nvars-1)              */
  enum ADIOS_DATATYPES type; /* type of variable */

  int from_steps; /* the first timestep in the returned data             */
  int nsteps;     /* the number of timesteps in the returned data        */
  H5_adios_bb_sel *sel; /* sub-selection of requested selection */
  void *data;           /* pointer to data, at next adios_read_check() memory
                           will likely be overwritten                          */
} H5_adios_varchunk_t;

enum ADIOS_DATATYPES adios_transform_get_var_original_type_index(
    struct adios_index_comp_struct_v1 *var);

enum ADIOS_DATATYPES adios_transform_get_var_original_type_var_header(
    H5_bp_var_data_header_t *var_header);

enum ADIOS_YES_NOFLAG checkSubfile(H5_bp_file_t *fp);

//
// util functions
//
const char *print_adios_type(int type);
void not_supported(const char *msg);
uint64_t adios_get_type_size(enum ADIOS_DATATYPES type, const void *var);

uint64_t adios_get_stat_size(void *data, enum ADIOS_DATATYPES type,
                             enum ADIOS_STAT stat_id);

const char *bp_value_to_string(enum ADIOS_DATATYPES type, void *data);

int adios_get_var_ndim(H5_bp_file_t *fp, int idx,
                       struct adios_index_comp_struct_v1 *varInAdios,
                       unsigned int *inqIsGlobal);

int adios_get_attr_ndim(struct adios_index_comp_struct_v1 *attrInAdios);

int adios_get_var_dimensions(
    const struct adios_index_characteristic_dims_struct_v1 *dims, int ndim,
    uint64_t *ldims, uint64_t *gdims, uint64_t *offsets,
    enum ADIOS_YES_NOFLAG is_fortran_order);

int adios_util_check_selection(int ndim, H5_bp_request_t *req, uint64_t *gdims,
                               uint64_t *ldims, uint64_t *offsets,
                               uint64_t *payload_size);
void adios_util_count_partial_pg(const int ndim, uint64_t *size_in_dset,
                                 uint64_t *offset_in_dset,
                                 uint64_t *offset_in_var, const uint64_t *start,
                                 const uint64_t *count, const uint64_t *ldims,
                                 const uint64_t *offset);

//
// adios functions
//
void mpi_OPS1_readvalue(uint64_t slice_size, uint64_t slice_offset,
                        H5_bp_file_t *fh);

// const char * adios_type_to_string_int (int type);

// NCSU - Returns count of the set of characteristcs for a variable

uint8_t adios_get_stat_set_count(enum ADIOS_DATATYPES type);

// void H5Adios_list_MPI_Error(int err);
int handle_char_stat(H5_bp_file_t *b, enum ADIOS_DATATYPES original_var_type,
                     struct adios_index_characteristic_struct_v1 *chaStruc);

// void call_adios_get_version(char* buff, H5_bp_file_t* bpFile, int offset);
int call_adios_get_version(H5_bp_file_t *bpFile, int offset);

int call_init(H5_bp_file_t *bpFile, MPI_File fh);

void call_close(H5_bp_file_t *bpFile);

int call_seek(H5_bp_file_t *bpFile, hsize_t start);

int call_adios_has_var(H5_bp_file_t *bpFile, const char *name);

int call_seek_back(H5_bp_file_t *bpFile, hsize_t start);

// helpers:
void readStat(struct adios_index_characteristic_struct_v1 *cc,
              enum ADIOS_STAT statType, uint32_t size, H5_bp_file_t *b);
int readChar(H5_bp_file_t *bpFile, char *var);
// int readInt(int bytes, H5_bp_file_t* bpFile, uint64_t* pos, uint64_t* var);
void readInt8(H5_bp_file_t *bpFile, uint64_t *var);
void readInt4(H5_bp_file_t *bpFile, uint32_t *var);
void readInt2(H5_bp_file_t *bpFile, uint16_t *var);
void readInt1(H5_bp_file_t *bpFile, uint8_t *var);

void readString(char **strVal, int strLen, H5_bp_file_t *bpFile);
double read_double(H5_bp_file_t *bpFile);

int read_var_bb(H5_bp_file_t *fp, H5_bp_request_t *req);
void read_scalar(H5_bp_request_t *req,
                 int64_t start_idx, // int size_of_type,
                 enum ADIOS_YES_NOFLAG has_subfile, H5_bp_file_t *fh);

#endif // ADIOS_READER_H
