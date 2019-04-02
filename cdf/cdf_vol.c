/*
 * Purpose:     This is a CDF file-format VOL connector.
 */

/* Header files needed */
/* (Public HDF5 and standard C / POSIX only) */
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include "hdf5.h"
#include "cdf_vol.h"

/**********/
/* Macros */
/**********/

/* Whether to display log messages when callback is invoked */
/* (Uncomment to enable) */
//#define ENABLE_CDF_LOGGING
//#define ENABLE_CDF_VERBOSE
#ifdef H5_HAVE_PARALLEL
#define H5_VOL_HAVE_PARALLEL  /* Comment if you want to turn off MPI-IO */
#endif

/************************* */
/* Borrowed from ADIOS VOL */
/************************* */

/* Structure for data flattening */
typedef struct {
	hsize_t posCompact;  // Element index in intermediate compact flattened buffer
	hsize_t posInSource; // Element index in HDF5 file/memory
} H5_posMap;

/* Helper comparison function for sorting */
int cmpfuncWobble(const void *a, const void *b) {
	H5_posMap *wa = (H5_posMap *)a;
	H5_posMap *wb = (H5_posMap *)b;
	return (wa->posInSource - wb->posInSource);
}
void traverseBlock(int ndim, int currDim, hsize_t *start, hsize_t *count, hsize_t *m, hsize_t previous, H5_posMap **result);
void getMultiplier(int ndim, hsize_t *dims, hsize_t *result);
void assignToMemSpace(H5_posMap *sourceSelOrderInC, H5_posMap *targetSelOrderInC, hsize_t npoints, size_t dataTypeSize, char *compactData, char *buf);
static void GetSelOrder(hid_t space_id, H5S_sel_type space_type, H5_posMap **result);

/************/
/* Typedefs */
/************/

typedef enum {
	NC_BYTE=1,	//=	\x00 \x00 \x00 \x01	//8-bit signed integers
	NC_CHAR,		//=	\x00 \x00 \x00 \x02	//text characters
	NC_SHORT,		//=	\x00 \x00 \x00 \x03	//16-bit signed integers
	NC_INT,			//=	\x00 \x00 \x00 \x04	//32-bit signed integers
	NC_FLOAT,		//=	\x00 \x00 \x00 \x05	//IEEE single precision floats
	NC_DOUBLE,	//=	\x00 \x00 \x00 \x06	//IEEE double precision floats
	NC_UBYTE,		//=	\x00 \x00 \x00 \x07	//unsigned 1 byte integer
	NC_USHORT,	//=	\x00 \x00 \x00 \x08	//unsigned 2-byte integer
	NC_UINT,		//=	\x00 \x00 \x00 \x09	//unsigned 4-byte integer
	NC_INT64,		//=	\x00 \x00 \x00 \x0A	//signed 8-byte integer
	NC_UINT64,	//=	\x00 \x00 \x00 \x0B	//unsigned 8-byte integer
} cdf_nc_type_t;

/* Assume non_neg and offsets are unsigned 64-bit int */
typedef uint64_t cdf_non_neg_t;
typedef uint64_t cdf_offset_t;

/* CDF name object (structure) */
typedef struct cdf_name_t {
	cdf_non_neg_t npadded; /* Length of name string aligned to 4-byte boundary */
	cdf_non_neg_t nelems; /* Length (number of bytes) of name string */
	char* string;
} cdf_name_t;

/* CDF dimension object (structure) */
typedef struct cdf_dim_t {
	cdf_name_t *name;
	cdf_non_neg_t length;
} cdf_dim_t;

/* CDF attribute object (structure) */
typedef struct cdf_att_t {
	cdf_name_t *name;
	void *file; /* "void" Pointer back to parent file object */
	cdf_nc_type_t nc_type;
	cdf_non_neg_t nvals;
	void *values;
	int nc_type_size;
} cdf_att_t;

/* CDF variable object (structure) */
typedef struct cdf_var_t {
	cdf_name_t *name;
	void *file; /* "void" Pointer back to parent file object */
	cdf_non_neg_t offset; /* Offset byte in file for this variable (begin) */
	cdf_non_neg_t vsize; /* vsize */
	cdf_nc_type_t nc_type; /* var data type */
	int nc_type_size;
	cdf_non_neg_t dimrank; /* Dimensionality (rank) of this variable */
	cdf_non_neg_t *dimids; /* Dimension ID (index into dim_list) for variable
                          * shape. We say this is a "record variable" if and only
                          * if the first dimension is the record dimension.  */
	hsize_t *dimlens; /* Size of each dimension */
	cdf_non_neg_t natts; /* number of attributes for this var */
	cdf_att_t *atts; /* list of cdf_att_t objects for this var */
} cdf_var_t;

/* The CDF VOL info object */
typedef struct H5VL_cdf_t {
	char *fname; /* Store the file name */
	uint8_t fmt; /* CDF File Spec (1,2,5) */
	cdf_non_neg_t numrecs; /* length of record dimension */
	cdf_non_neg_t ndims; /* number of dimensions */
	cdf_non_neg_t natts; /* number of attributes */
	cdf_non_neg_t nvars; /* number of variables */
	cdf_dim_t *dims; /* list of cdf_dim_t objects */
	cdf_att_t *atts; /* list of cdf_att_t objects */
	cdf_var_t *vars; /* list of cdf_var_t objects */
	int rank; /* MPI Rank */
#ifdef H5_VOL_HAVE_PARALLEL
	MPI_File fh; /* MPI File Handle */
	MPI_Comm comm; /* MPI Communicator */
	MPI_Info info; /* MPI Info */
	MPI_Offset size; /* Total Size of read-only file */
	char *cache; /* In-memory cache for parsing header (size = H5VL_CDF_CACHE_SIZE) */
	MPI_Offset cache_offset; /* What is the file offset for first byte of cache */
#else
	int fd; /* posix file object */
#endif
} H5VL_cdf_t;

/* VOL Info Object */
typedef struct H5VL_cdf_info_t {
	hid_t vol_id; /* VOL ID */
} H5VL_cdf_info_t;

/* Keep track of open files */
int n_file_objects_open = 0;
H5VL_cdf_t *file_objects_open[MAX_FILES_OPEN];

/********************* */
/* Function prototypes */
/********************* */

/* Helper routines */
static void cdf_handle_error(char* errmsg, int lineno);
#ifdef H5_VOL_HAVE_PARALLEL
static H5VL_cdf_t *H5VL_cdf_new_obj_mpio(MPI_File fh, hid_t acc_tpl);
#else
static H5VL_cdf_t *H5VL_cdf_new_obj(int fd);
#endif
static herr_t H5VL_cdf_free_obj(H5VL_cdf_t *obj);
static void bytestr_rev(char *p, int size);
static uint32_t bytestr_to_uint32_t(uint8_t *bytestr);
static uint64_t bytestr_to_uint64_t(uint8_t *bytestr);
static uint32_t get_uint32_t(H5VL_cdf_t* file, off_t *rdoff);
static uint64_t get_uint64_t(H5VL_cdf_t* file, off_t *rdoff);
static cdf_non_neg_t get_non_neg(H5VL_cdf_t* file, off_t *rdoff);
static cdf_offset_t get_offset(H5VL_cdf_t* file, off_t *rdoff);
static uint8_t get_cdf_spec(H5VL_cdf_t* file, off_t *rdoff);
static cdf_name_t *get_cdf_name_t(H5VL_cdf_t* file, off_t *rdoff);
static int get_nc_type_size(cdf_nc_type_t nc_type);
static void *get_att_values(H5VL_cdf_t* file, off_t *rdoff, cdf_att_t *atts, int iatt);
static cdf_dim_t *get_cdf_dims_t(H5VL_cdf_t* file, off_t *rdoff);
static cdf_att_t *get_cdf_atts_t(H5VL_cdf_t* file, off_t *rdoff, cdf_non_neg_t natts);
static void get_header_list_item(H5VL_cdf_t* file, off_t *rdoff, cdf_var_t* var);
static herr_t cdf_read_header(H5VL_cdf_t* file);
static herr_t cdf_cache_read(H5VL_cdf_t* file, off_t *rdoff, char *buf, MPI_Offset count);

/* "Management" callbacks */
static herr_t H5VL_cdf_init(hid_t vipl_id);
static herr_t H5VL_cdf_term(void);

/* Attribute callbacks */
void *H5VL_cdf_attr_open(void *obj, const H5VL_loc_params_t *loc_params, const char *attr_name, hid_t aapl_id, hid_t dxpl_id, void **req);
herr_t H5VL_cdf_attr_read(void *obj, hid_t dtype_id, void *buf, hid_t dxpl_id, void **req);
herr_t H5VL_cdf_attr_close(void *obj, hid_t dxpl_id, void **req);

/* Dataset callbacks */
void *H5VL_cdf_dataset_open(void *obj, const H5VL_loc_params_t *loc_params, const char *name, hid_t dapl_id, hid_t dxpl_id, void **req);
static herr_t H5VL_cdf_dataset_read(void *dset, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t plist_id, void *buf, void **req);
static herr_t H5VL_cdf_dataset_get(void *obj, H5VL_dataset_get_t get_type, hid_t dxpl_id, void **req, va_list arguments);
static herr_t H5VL_cdf_dataset_close(void *dset, hid_t dxpl_id, void **req);

/* Datatype callbacks */

/* File callbacks */
static void *H5VL_cdf_file_open(const char *name, unsigned flags, hid_t fapl_id, hid_t dxpl_id, void **req);
static herr_t H5VL_cdf_file_close(void *file, hid_t dxpl_id, void **req);

/* Group callbacks */
static void *H5VL_cdf_group_open(void *obj, const H5VL_loc_params_t *loc_params, const char *varName, hid_t gapl_id, hid_t dxpl_id, void **req);
static herr_t H5VL_cdf_group_get(void *obj, H5VL_group_get_t get_type, hid_t dxpl_id, void **req, va_list arguments);

/* Link callbacks */

/* Object callbacks */

/* Async request callbacks */

/*******************/
/* Local variables */
/*******************/

/* Pass through VOL connector class struct */
const H5VL_class_t H5VL_cdf_g = {
    H5VL_CDF_VERSION,                      /* version      */
    (H5VL_class_value_t)H5VL_CDF_VALUE,    /* value        */
    H5VL_CDF_NAME,                         /* name         */
    0,                                      /* capability flags */
    H5VL_cdf_init,                         /* initialize   */
    H5VL_cdf_term,                         /* terminate    */
    sizeof(H5VL_cdf_info_t),               /* info size    */
    NULL, //H5VL_cdf_info_copy,                    /* info copy    */
    NULL, //H5VL_cdf_info_cmp,                     /* info compare */
    NULL, //H5VL_cdf_info_free,                    /* info free    */
    NULL, //H5VL_cdf_info_to_str,                  /* info to str  */
    NULL, //H5VL_cdf_str_to_info,                  /* str to info  */
    NULL, //H5VL_cdf_get_object,                   /* get_object   */
    NULL, //H5VL_cdf_get_wrap_ctx,                 /* get_wrap_ctx */
    NULL, //H5VL_cdf_wrap_object,                  /* wrap_object  */
    NULL, //H5VL_cdf_free_wrap_ctx,                /* free_wrap_ctx */
    {                                           /* attribute_cls */
        NULL, //H5VL_cdf_attr_create,                       /* create */
        H5VL_cdf_attr_open,                         /* open */
        H5VL_cdf_attr_read,                         /* read */
        NULL, //H5VL_cdf_attr_write,                        /* write */
        NULL, //H5VL_cdf_attr_get,                          /* get */
        NULL, //H5VL_cdf_attr_specific,                     /* specific */
        NULL, //H5VL_cdf_attr_optional,                     /* optional */
        H5VL_cdf_attr_close                         /* close */
    },
    {                                           /* dataset_cls */
        NULL, //H5VL_cdf_dataset_create,                    /* create */
        H5VL_cdf_dataset_open,                      /* open */
        H5VL_cdf_dataset_read,                      /* read */
        NULL, //H5VL_cdf_dataset_write,                     /* write */
        H5VL_cdf_dataset_get,                       /* get */
        NULL, //H5VL_cdf_dataset_specific,                  /* specific */
        NULL, //H5VL_cdf_dataset_optional,                  /* optional */
        H5VL_cdf_dataset_close                      /* close */
    },
    {                                               /* datatype_cls */
        NULL, //H5VL_cdf_datatype_commit,                   /* commit */
        NULL, //H5VL_cdf_datatype_open,                     /* open */
        NULL, //H5VL_cdf_datatype_get,                      /* get_size */
        NULL, //H5VL_cdf_datatype_specific,                 /* specific */
        NULL, //H5VL_cdf_datatype_optional,                 /* optional */
        NULL, //H5VL_cdf_datatype_close                     /* close */
    },
    {                                           /* file_cls */
        NULL, //H5VL_cdf_file_create,                       /* create */
        H5VL_cdf_file_open,                         /* open */
        NULL, //H5VL_cdf_file_get,                          /* get */
        NULL, //H5VL_cdf_file_specific,                     /* specific */
        NULL, //H5VL_cdf_file_optional,                     /* optional */
        H5VL_cdf_file_close                         /* close */
    },
    {                                           /* group_cls */
        NULL, //H5VL_cdf_group_create,                      /* create */
        H5VL_cdf_group_open,                        /* open */
        H5VL_cdf_group_get,                         /* get */
        NULL, //H5VL_cdf_group_specific,                    /* specific */
        NULL, //H5VL_cdf_group_optional,                    /* optional */
        NULL, //H5VL_cdf_group_close                        /* close */
    },
    {                                           /* link_cls */
        NULL, //H5VL_cdf_link_create,                       /* create */
        NULL, //H5VL_cdf_link_copy,                         /* copy */
        NULL, //H5VL_cdf_link_move,                         /* move */
        NULL, //H5VL_cdf_link_get,                          /* get */
        NULL, //H5VL_cdf_link_specific,                     /* specific */
        NULL, //H5VL_cdf_link_optional,                     /* optional */
    },
    {                                           /* object_cls */
        NULL, //H5VL_cdf_object_open,                       /* open */
        NULL, //H5VL_cdf_object_copy,                       /* copy */
        NULL, //H5VL_cdf_object_get,                        /* get */
        NULL, //H5VL_cdf_object_specific,                   /* specific */
        NULL, //H5VL_cdf_object_optional,                   /* optional */
    },
    {                                           /* request_cls */
        NULL, //H5VL_cdf_request_wait,                      /* wait */
        NULL, //H5VL_cdf_request_notify,                    /* notify */
        NULL, //H5VL_cdf_request_cancel,                    /* cancel */
        NULL, //H5VL_cdf_request_specific,                  /* specific */
        NULL, //H5VL_cdf_request_optional,                  /* optional */
        NULL, //H5VL_cdf_request_free                       /* free */
    },
    NULL                                        /* optional */
};

/* The connector identification number, initialized at runtime */
static hid_t H5VL_CDF_g = H5I_INVALID_HID;

/*-------------------------------------------------------------------------
 * Function:    traverseBlock
 *
 * Purpose:     Recursive function used to generate a flattened representation
 *              of the hyperslab selection. (Borrowed from ADIOS VOL)
 *
 * Return:      Void
 *-------------------------------------------------------------------------
 */
void
traverseBlock(int ndim, int currDim, hsize_t *start, hsize_t *count,
                   hsize_t *m, hsize_t previous, H5_posMap **result) {
	hsize_t k, n, pos, nextCompactPos, base;
	if (ndim == currDim + 1) {
		for (k = 0; k < count[currDim]; k++) {
			pos = start[currDim] + k;
			(**result).posInSource = previous + pos * m[currDim];
			nextCompactPos = (**result).posCompact + 1;
			(*result)++;
			(*result)->posCompact = nextCompactPos;
		}
		return;
	}

	for (k = start[currDim]; k < start[currDim] + count[currDim]; k++) {
		base = k * m[currDim] + previous;
		traverseBlock(ndim, currDim + 1, start, count, m, base, result);
	}
}

/*-------------------------------------------------------------------------
 * Function:    getMultiplier
 *
 * Purpose:     Calculate the simple "multiplier" array for generating a
 *              flat offset from array coordinates. (Borrowed from ADIOS VOL)
 *
 * Return:      Void
 *-------------------------------------------------------------------------
 */
void
getMultiplier(int ndim, hsize_t *dims, hsize_t *result) {
	int n, k;
	for (n = 0; n < ndim; n++) {
		result[n] = 1;
		for (k = n + 1; k < ndim; k++) {
			result[n] *= dims[k];
		}
	}
}

/*-------------------------------------------------------------------------
 * Function:    assignToMemSpace
 *
 * Purpose:     Copy data from the contiguous/compact array (compactData)
 *              to the target selection (buf). (Borrowed from ADIOS VOL)
 *
 * Return:      Void
 *-------------------------------------------------------------------------
 */
void assignToMemSpace(H5_posMap *sourceSelOrderInC,
                      H5_posMap *targetSelOrderInC, hsize_t npoints,
                      size_t dataTypeSize, char *compactData, char *buf) {
  size_t k;
  hsize_t n;

  if (sourceSelOrderInC == NULL) { // fileSpace=ALL memSpace=partial
    for (n = 0; n < npoints; n++) {
      hsize_t dd = n;
      hsize_t ss = targetSelOrderInC[n].posInSource;
      for (k = 0; k < dataTypeSize; k++) {
        buf[ss * dataTypeSize + k] = compactData[dd * dataTypeSize + k];
      }
    }
    return;
  }
  if (targetSelOrderInC == NULL) {
    for (n = 0; n < npoints; n++) {
      hsize_t dd = sourceSelOrderInC[n].posCompact;
      hsize_t ss = n;
      for (k = 0; k < dataTypeSize; k++) {
        buf[ss * dataTypeSize + k] = compactData[dd * dataTypeSize + k];
      }
    }
    return;
  }

  for (n = 0; n < npoints; n++) {
    hsize_t ss = targetSelOrderInC[n].posInSource;
    hsize_t dd = sourceSelOrderInC[n].posCompact;
    //printf("Copying index %llu in Compact to %llu in buf.\n",dd,ss);
    for (k = 0; k < dataTypeSize; k++) {
      buf[ss * dataTypeSize + k] = compactData[dd * dataTypeSize + k];
    }
  }
}

/*-------------------------------------------------------------------------
 * Function:    GetSelOrder
 *
 * Purpose:     Populates the H5_posMap position-map array
 *              (Function borrowed from ADIOS VOL)
 *
 * Return:      Populates pointer to *result
 *-------------------------------------------------------------------------
 */
static void
GetSelOrder(hid_t space_id, H5S_sel_type space_type, H5_posMap **result) {

	int i, n, ndim;
	hsize_t nblocks, total;
	H5_posMap *startPointer;
	hsize_t *blockinfo;
	herr_t status;

	if (space_type != H5S_SEL_HYPERSLABS) {
		/* Not set up for point selection yet.. */
		return;
	}

	ndim = H5Sget_simple_extent_ndims(space_id);
	hsize_t dims[ndim];
	H5Sget_simple_extent_dims(space_id, dims, NULL);

	hsize_t mm[ndim];
	getMultiplier(ndim, dims, mm);

	nblocks = H5Sget_select_hyper_nblocks(space_id);
	total = H5Sget_select_npoints(space_id);
	*result = (H5_posMap *)malloc(sizeof(H5_posMap) * (total + 1));
	(*result)->posCompact = 0;

	blockinfo = (hsize_t *) malloc(sizeof(hsize_t) * 2 * ndim * nblocks);
	status = H5Sget_select_hyper_blocklist(space_id, (hsize_t)0, nblocks, blockinfo);
	startPointer = *result;

	i = 0; n = 0;
	for (n = 0; n < nblocks; n++) {
		uint64_t blockSize = 1;
		uint64_t h_start[ndim], h_count[ndim];
		for (i = 0; i < ndim; i++) {
			int pos = 2 * ndim * n;
			h_start[i] = blockinfo[pos + i];
			h_count[i] = blockinfo[pos + ndim + i] - h_start[i] + 1;
			blockSize *= h_count[i];
		}
		//printf(" got block: start[%llu, %llu, %llu], count[%llu, %llu, %llu], mm[%llu, %llu, %llu]\n", h_start[0], h_start[1], h_start[2], h_count[0], h_count[1], h_count[2], mm[0], mm[1], mm[2]);
		traverseBlock(ndim, 0, h_start, h_count, mm, 0, result);
	}

	/* Need to reorder the positions */
	qsort(startPointer, total, sizeof(H5_posMap), cmpfuncWobble);

	*result = startPointer;
	free(blockinfo);

}


#ifdef H5_VOL_HAVE_PARALLEL
/*-------------------------------------------------------------------------
 * Function:    cdf_cache_read
 *
 * Purpose:     Get required data from cache - If needed, advance the cache by
 *              reading to the cache on rank 0 and broadcasting the data to
 *              the other ranks.
 *
 * Return:      Success:    0
 *              Failure:    -1
 *-------------------------------------------------------------------------
 */
static herr_t
cdf_cache_read(H5VL_cdf_t* file, off_t *rdoff, char *buf, MPI_Offset count)
{
	MPI_Status status;
	MPI_Offset extent_want, extent_have;
	MPI_Offset offset;
	MPI_Offset rdsize;
	herr_t retval=0;

	if (count > H5VL_CDF_CACHE_SIZE) {
		printf("ERROR!! Cannot handle count = %llu.\n",count);
	}

	if (file->cache_offset < 0) {

			rdsize = H5VL_CFD_MIN( H5VL_CDF_CACHE_SIZE, (file->size)-(file->cache_offset) );
			//printf("Initializing cache to size = %llu\n", rdsize);
			file->cache_offset = 0;
			if (file->rank==0)
				MPI_File_read_at( file->fh, file->cache_offset, file->cache, rdsize, MPI_BYTE, &status );
			MPI_Bcast( file->cache, rdsize, MPI_BYTE, 0, file->comm );

	} else {

		extent_want = ((MPI_Offset)(*rdoff)) + count;
		extent_have = file->cache_offset + H5VL_CDF_CACHE_SIZE;
		if( extent_want > extent_have ) {

			rdsize = H5VL_CFD_MIN( H5VL_CDF_CACHE_SIZE, (file->size)-(file->cache_offset) );
			printf("Advancing cache by %llu\n", rdsize);
			file->cache_offset = file->cache_offset + H5VL_CDF_CACHE_SIZE;

			/* Advance cache by reading first H5VL_CDF_CACHE_SIZE from file,
			 * and bcasting the data to the other ranks.
			 */
			if (file->rank==0)
				MPI_File_read_at( file->fh, file->cache_offset, file->cache, rdsize, MPI_BYTE, &status );
			MPI_Bcast( file->cache, rdsize, MPI_BYTE, 0, file->comm );

		}

	}

	/* "Read" (memcopy) necessary data */
	offset = (MPI_Offset)(*rdoff) - file->cache_offset;
	memcpy(buf, &file->cache[offset], count);

	return retval;
}
#endif

/*-------------------------------------------------------------------------
 * Function:    cdf_handle_error
 *
 * Purpose:     Error-handler helper routine.
 *
 * Return:      Void.
 *-------------------------------------------------------------------------
 */
static void
cdf_handle_error(char* errmsg, int lineno)
{
	fprintf(stderr, "CDF VOL Error at line %d of %s: %s\n", lineno, __FILE__, errmsg);
	//MPI_Abort(MPI_COMM_WORLD, 1);
}

/*-------------------------------------------------------------------------
 * Function:    bytestr_rev
 *
 * Purpose:     Reverse the order of a byte string in memory
 *
 * Return:      void
 *-------------------------------------------------------------------------
 */
static void
bytestr_rev(char *p, int size)
{
  char *q = p;
  //while(q && *q) ++q;
	q+=size;
  for(--q; p < q; ++p, --q)
    *p = *p ^ *q,
    *q = *p ^ *q,
    *p = *p ^ *q;
}

/*-------------------------------------------------------------------------
 * Function:    bytestr_to_uint32_t
 *
 * Purpose:     Convert 4-character byte string to uint32_t (Big Endian)
 *
 * Return:      uint32_t value.
 *-------------------------------------------------------------------------
 */
static uint32_t
bytestr_to_uint32_t(uint8_t *bytestr)
{
	/* Little Endian Case */
	//uint32_t myInt1 = (uint32_t)bytestr[0] + ((uint32_t)bytestr[1] << 8) + ((uint32_t)bytestr[2] << 16) + ((uint32_t)bytestr[3] << 24);
	/* Big Endian Case */
	uint32_t myInt1 = ((uint32_t)bytestr[0] << 24) +
	                  ((uint32_t)bytestr[1] << 16) +
	                  ((uint32_t)bytestr[2] << 8) +
	                  (uint32_t)bytestr[3];
	return myInt1;
}

/*-------------------------------------------------------------------------
 * Function:    bytestr_to_uint64_t
 *
 * Purpose:     Convert 8-character byte string to uint64_t (Big Endian)
 *
 * Return:      uint64_t value.
 *-------------------------------------------------------------------------
 */
static uint64_t
bytestr_to_uint64_t(uint8_t *bytestr)
{
	/* Big Endian Case */
	uint64_t myInt1 = ((uint64_t)(bytestr[0]) << 56) +
	                  ((uint64_t)(bytestr[1]) << 48) +
	                  ((uint64_t)(bytestr[2]) << 40) +
	                  ((uint64_t)(bytestr[3]) << 32) +
	                  ((uint64_t)(bytestr[4]) << 24) +
	                  ((uint64_t)(bytestr[5]) << 16) +
	                  ((uint64_t)(bytestr[6]) << 8) +
	                  (uint64_t)(bytestr[7]);
	return myInt1;
}

/*-------------------------------------------------------------------------
 * Function:    get_uint32_t
 *
 * Purpose:     Parse file header to get an uint32_t value from 4 consecutive
 *              bytes. File storage is "Big Endian"
 *
 * Return:      uint32_t value.
 *-------------------------------------------------------------------------
 */
static uint32_t
get_uint32_t(H5VL_cdf_t* file, off_t *rdoff)
{
	char bytestr[4];
#ifdef H5_VOL_HAVE_PARALLEL
		cdf_cache_read(file, rdoff, &bytestr[0], 4);
#else
		pread(file->fd, &bytestr[0], 4, *rdoff);
#endif
	*rdoff = *rdoff + 4;
#ifdef ENABLE_CDF_VERBOSE
	printf(" bytes: %d %d %d %d\n", (uint8_t)bytestr[0], (uint8_t)bytestr[1], (uint8_t)bytestr[2], (uint8_t)bytestr[3]);
#endif
	return (uint64_t)bytestr_to_uint32_t( (uint8_t *)bytestr );
}

/*-------------------------------------------------------------------------
 * Function:    get_uint64_t
 *
 * Purpose:     Parse file header to get an uint64_t value from 8 consecutive
 *              bytes. File storage is "Big Endian"
 *
 * Return:      uint64_t value.
 *-------------------------------------------------------------------------
 */
static uint64_t
get_uint64_t(H5VL_cdf_t* file, off_t *rdoff)
{
	char bytestr[8];
#ifdef H5_VOL_HAVE_PARALLEL
		cdf_cache_read(file, rdoff, &bytestr[0], 8);
#else
		pread(file->fd, &bytestr[0], 8, *rdoff);
#endif
	*rdoff = *rdoff + 8;
	return bytestr_to_uint64_t( (uint8_t *)bytestr );
}

/*-------------------------------------------------------------------------
 * Function:    get_cdf_spec
 *
 * Purpose:     Parse file header to get the CDF-? specification version.
 *
 * Return:      uint8_t Specification version.
 *-------------------------------------------------------------------------
 */
static uint8_t
get_cdf_spec(H5VL_cdf_t* file, off_t *rdoff)
{
	char bytestr[4];
#ifdef H5_VOL_HAVE_PARALLEL
		cdf_cache_read(file, rdoff, &bytestr[0], 4);
#else
		pread(file->fd, &bytestr[0], 4, *rdoff);
#endif
	*rdoff = *rdoff + 4;
	file->fmt = (uint8_t) bytestr[3];
#ifdef ENABLE_CDF_LOGGING
	printf("------- (FILE FORMAT = %c%c%c-%d) \n", bytestr[0], bytestr[1], bytestr[2], (uint8_t) bytestr[3]);
#endif
	return (uint8_t) bytestr[3];
}

/*-------------------------------------------------------------------------
 * Function:    get_non_neg
 *
 * Purpose:     Parse file header to get a NON_NEG value.
 *
 * Return:      cdf_non_neg_t value.
 *-------------------------------------------------------------------------
 */
static cdf_non_neg_t
get_non_neg(H5VL_cdf_t* file, off_t *rdoff)
{
	/* Read length of record dimension (numrecs) */
	if (file->fmt < 5) {
		/* NON_NEG = <non-negative INT> = <32-bit signed integer, Bigendian, two's complement> */
		return (cdf_non_neg_t) get_uint32_t(file,rdoff);
	} else {
		/* NON_NEG = <non-negative INT64> = <64-bit signed integer, Bigendian, two's complement> */
		return (cdf_non_neg_t) get_uint64_t(file,rdoff);
	}
}

/*-------------------------------------------------------------------------
 * Function:    get_offset
 *
 * Purpose:     Parse file header to get the offset value for a variable
 *
 * Return:      cdf_offset_t value.
 *-------------------------------------------------------------------------
 */
static cdf_offset_t
get_offset(H5VL_cdf_t* file, off_t *rdoff)
{
	/* Read length of record dimension (numrecs) */
	if (file->fmt < 2) {
		/* NON_NEG = <non-negative INT> = <32-bit signed integer, Bigendian, two's complement> */
		return (cdf_offset_t) get_uint32_t(file,rdoff);
	} else {
		/* NON_NEG = <non-negative INT64> = <64-bit signed integer, Bigendian, two's complement> */
		return (cdf_offset_t) get_uint64_t(file,rdoff);
	}
}

/*-------------------------------------------------------------------------
 * Function:    cdf_name_t
 *
 * Purpose:     Parse file header to get a cdf_name_t object
 *
 * Return:      Pointer to single cdf_name_t object.
 *-------------------------------------------------------------------------
 */
static cdf_name_t*
get_cdf_name_t(H5VL_cdf_t* file, off_t *rdoff) {

	cdf_name_t *name;
	name = (cdf_name_t *)calloc(1, sizeof(cdf_name_t));

	/* Populate nelems */
	name->nelems = get_non_neg(file, rdoff);
#ifdef ENABLE_CDF_VERBOSE
	printf("nelems check = %llu\n", name->nelems);
#endif

	/* Populate padded byte count */
	name->npadded = (cdf_non_neg_t)( (name->nelems/4)*4 + ((name->nelems%4>0)?4:0) );
#ifdef ENABLE_CDF_VERBOSE
	printf("padded byte count = %d\n", (int)name->npadded);
#endif

	/* Populate string of characters for the name */
	name->string = (char *) malloc((name->nelems) * sizeof(char));
#ifdef H5_VOL_HAVE_PARALLEL
	cdf_cache_read(file, rdoff, name->string, name->nelems);
#else
	pread(file->fd, name->string, name->nelems, *rdoff);
#endif
	*rdoff = *rdoff + name->npadded;

#ifdef ENABLE_CDF_VERBOSE
	printf("namestr = ");
	for (int c=0; c<name->nelems; c++){
		printf("%c", name->string[c]);
	}
	printf("\n");
#endif

	return name;
}

/*-------------------------------------------------------------------------
 * Function:    get_cdf_dims_t
 *
 * Purpose:     Parse file header to generate list of cdf_dim_t objects.
 *
 * Return:      Pointer to cdf_dim_t list.
 *-------------------------------------------------------------------------
 */
static cdf_dim_t*
get_cdf_dims_t(H5VL_cdf_t* file, off_t *rdoff) {

	cdf_dim_t *dims;
	dims = (cdf_dim_t *) calloc(file->ndims, sizeof(cdf_dim_t));

	/* Read through each of file->ndims dimension entries */
	for(int idim=0;idim<file->ndims;idim++){
		/* Get dimension name */
		dims[idim].name = get_cdf_name_t(file, rdoff);
		/* Get dimension length */
		dims[idim].length = get_non_neg(file, rdoff);
	}

	return dims;
}

/*-------------------------------------------------------------------------
 * Function:    get_nc_type_size
 *
 * Purpose:     Given the nc_type - return the number of bytes
 *
 * Return:      int (number of bytes for each element)
 *-------------------------------------------------------------------------
 */
static int
get_nc_type_size(cdf_nc_type_t nc_type) {

	int nc_type_size=0;

	switch(nc_type) {
		case 	NC_BYTE:	//=	\x00 \x00 \x00 \x01	//8-bit signed integers
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_BYTE\n");
#endif
			nc_type_size = 1;
			break;
		case 	NC_CHAR:		//=	\x00 \x00 \x00 \x02	//text characters
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_CHAR\n");
#endif
			nc_type_size = 1;
			break;
		case 	NC_SHORT:		//=	\x00 \x00 \x00 \x03	//16-bit signed integers
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_SHORT\n");
#endif
			nc_type_size = 2;
			break;
		case 	NC_INT:			//=	\x00 \x00 \x00 \x04	//32-bit signed integers
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_INT\n");
#endif
			nc_type_size = 4;
			break;
		case 	NC_FLOAT:		//=	\x00 \x00 \x00 \x05	//IEEE single precision floats
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_FLOAT\n");
#endif
			nc_type_size = 4;
			break;
		case 	NC_DOUBLE:	//=	\x00 \x00 \x00 \x06	//IEEE double precision floats
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_DOUBLE\n");
#endif
			nc_type_size = 8;
			break;
		case 	NC_UBYTE:		//=	\x00 \x00 \x00 \x07	//unsigned 1 byte integer
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_UBYTE\n");
#endif
			nc_type_size = 1;
			break;
		case 	NC_USHORT:	//=	\x00 \x00 \x00 \x08	//unsigned 2-byte integer
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_USHORT\n");
#endif
			nc_type_size = 2;
			break;
		case 	NC_UINT:		//=	\x00 \x00 \x00 \x09	//unsigned 4-byte integer
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_UINT\n");
#endif
			nc_type_size = 4;
			break;
		case 	NC_INT64:		//=	\x00 \x00 \x00 \x0A	//signed 8-byte integer
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_INT64\n");
#endif
			nc_type_size = 8;
			break;
		case 	NC_UINT64:	//=	\x00 \x00 \x00 \x0B	//unsigned 8-byte integer
#ifdef ENABLE_CDF_VERBOSE
			printf("NC_TYPE is NC_UINT64\n");
#endif
			nc_type_size = 8;
			break;
		default :
			printf("ERROR!! %d not a recognized nc_type\n", nc_type);
	}
	return nc_type_size;
}

/*-------------------------------------------------------------------------
 * Function:    get_att_values
 *
 * Purpose:     Parse file header to generate "value" buffer.
 *
 * Return:      Pointer to void "value" buffer.
 *-------------------------------------------------------------------------
 */
static void*
get_att_values(H5VL_cdf_t* file, off_t *rdoff, cdf_att_t *atts, int iatt) {

	void *values;

	atts[iatt].nc_type_size = get_nc_type_size(atts[iatt].nc_type);
	values = (void *) calloc(atts[iatt].nvals, atts[iatt].nc_type_size);

	/* Populate padded byte count */
	cdf_non_neg_t valsize = atts[iatt].nc_type_size * atts[iatt].nvals;

	cdf_non_neg_t npadded = (cdf_non_neg_t)( (valsize/4)*4 + ((valsize%4>0)?4:0) );
#ifdef ENABLE_CDF_VERBOSE
	printf("[%d] padded byte count = %llu\n", file->rank, npadded);
#endif

	/* Populate values from bytes in file */
#ifdef H5_VOL_HAVE_PARALLEL
	cdf_cache_read(file, rdoff, values, valsize);
#else
	pread(file->fd, values, valsize, *rdoff);
#endif
	*rdoff = *rdoff + npadded;


//#ifdef ENABLE_CDF_VERBOSE
	/* Check Values */
	//char *tmp = (char *) calloc(valsize, sizeof(char));
	//memcpy(tmp, values, valsize);
	//printf("att values = ");
	//for (int c=0; c<valsize; c++){
	//	printf("%c", tmp[c]);
	//}
	//printf("\n");
	//free(tmp);
//#endif

	return values;
}

/*-------------------------------------------------------------------------
 * Function:    get_cdf_atts_t
 *
 * Purpose:     Parse file header to generate list of cdf_att_t objects.
 *
 * Return:      Pointer to cdf_att_t list.
 *-------------------------------------------------------------------------
 */
static cdf_att_t*
get_cdf_atts_t(H5VL_cdf_t* file, off_t *rdoff, cdf_non_neg_t natts) {

	cdf_att_t *atts;
	atts = (cdf_att_t *) calloc(natts, sizeof(cdf_att_t));

	/* Read through each of file->natts dimension entries */
	for(int iatt=0;iatt<natts;iatt++){

		/* Get attribute name */
		atts[iatt].name = get_cdf_name_t(file, rdoff);

		/* Define a void pointer to parent file object */
		atts[iatt].file = (void *)file;

		/* Get attribute type */
		atts[iatt].nc_type = get_uint32_t(file,rdoff);
#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] atts[%d].nc_type = %d\n", file->rank, iatt, atts[iatt].nc_type);
#endif

		/* Get number of attribute values */
		atts[iatt].nvals = get_non_neg(file, rdoff);
#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] atts[%d].nvals = %llu\n", file->rank, iatt, atts[iatt].nvals);
#endif

		/* Read the attribute values */
		atts[iatt].values = get_att_values(file, rdoff, atts, iatt);

	}
	return atts;

}

/*-------------------------------------------------------------------------
 * Function:    get_cdf_vars_t
 *
 * Purpose:     Parse file header to generate list of cdf_var_t objects.
 *
 * Return:      Pointer to cdf_var_t list.
 *-------------------------------------------------------------------------
 */
static cdf_var_t*
get_cdf_vars_t(H5VL_cdf_t* file, off_t *rdoff) {

	cdf_var_t *vars;
	vars = (cdf_var_t *) calloc(file->nvars, sizeof(cdf_var_t));

	/* Read through each of file->nvars variable entries */
	for(int ivar=0;ivar<file->nvars;ivar++){

		/* Specification:
		 * var	=	name nelems [dimid ...] vatt_list nc_type vsize begin
		 */

		/* Get variable name */
		vars[ivar].name = get_cdf_name_t(file, rdoff);

		/* Get rank of this var (dimensionality) */
		vars[ivar].dimrank = get_non_neg(file, rdoff);
#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] vars[%d].dimrank = %llu\n", file->rank, ivar, vars[ivar].dimrank);
#endif
		if (vars[ivar].dimrank > 0) {
			vars[ivar].dimids = (cdf_non_neg_t *) calloc(vars[ivar].dimrank, sizeof(cdf_non_neg_t));
			vars[ivar].dimlens = (hsize_t *) calloc(vars[ivar].dimrank, sizeof(hsize_t));

			/* Get dimid list */
			for (int id=0; id<vars[ivar].dimrank; id++) {
				vars[ivar].dimids[id] = get_non_neg(file, rdoff);
				vars[ivar].dimlens[id] = file->dims[ vars[ivar].dimids[id] ].length;
#ifdef ENABLE_CDF_VERBOSE
				printf("[%d] vars[%d].dimids[%d] = %llu\n", file->rank, ivar, id, vars[ivar].dimids[id]);
				printf("[%d] vars[%d].dimlens[%d] = %llu\n", file->rank, ivar, id, vars[ivar].dimlens[id]);
#endif
			}

		}

		/* Populate vatt_list (attribute list for this var) */
		get_header_list_item(file, rdoff, &vars[ivar]);

		/* Get variable nc_type */
		vars[ivar].nc_type = get_uint32_t(file,rdoff);
		vars[ivar].nc_type_size = get_nc_type_size(vars[ivar].nc_type);
#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] vars[%d].nc_type = %d (size = %d)\n", file->rank, ivar, vars[ivar].nc_type, vars[ivar].nc_type_size);
#endif

		/* Get vsize */
		vars[ivar].vsize = get_non_neg(file, rdoff);
#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] vars[%d].vsize = %llu\n", file->rank, ivar, vars[ivar].vsize);
#endif

		/* Get file offset (begin) */
		vars[ivar].offset = get_offset(file, rdoff);
#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] vars[%d].offset = %llu\n", file->rank, ivar, vars[ivar].offset);
#endif

	}

	return vars;
}

/*-------------------------------------------------------------------------
 * Function:    get_header_list_item
 *
 * Purpose:     Parse a list item (dim_list, att_list, var_list)
 *              from the header. If var==NULL, the list item is attached to
 *              the parent H5VL_cdf_t (file) object.
 *
 * Return:      Void
 *-------------------------------------------------------------------------
 */
static void
get_header_list_item(H5VL_cdf_t* file, off_t *rdoff, cdf_var_t* var) {

	/* First check if first 4 bytes are 0 */
	uint32_t nc_item = get_uint32_t(file, rdoff);
	if ( nc_item > (uint32_t)0 ){


		/* "NC_DIMENSION" */
		if ( nc_item == (uint32_t)10 ) {

			/* Get nelems of dim_list (file->ndims) */
			file->ndims = (uint64_t) get_non_neg(file, rdoff);

			/* Populate list of dim objects */
			file->dims = get_cdf_dims_t(file, rdoff);

		/* "NC_VARIABLE" */
	} else if ( nc_item == (uint32_t)11 ) {

			/* Get nelems of var_list (file->nvars) */
			file->nvars = (uint64_t) get_non_neg(file, rdoff);
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] file->nvars = %llu\n", file->rank, file->nvars);
#endif

			/* Populate list of var objects */
			file->vars = get_cdf_vars_t(file, rdoff);

		/* "NC_ATTRIBUTE" */
		} else if ( nc_item == (uint32_t)12 ) {

			/* File attribute */
			if (var==NULL) {

				/* Get nelems of att_list (file->natts) */
				file->natts = (uint64_t) get_non_neg(file, rdoff);
				//printf("file->natts = %llu\n", file->natts);

				/* Populate list of att objects */
				file->atts = get_cdf_atts_t(file, rdoff, file->natts);

			/* Variable attribute */
			} else {

				/* Get nelems of att_list (file->natts) */
				var->natts = (uint64_t) get_non_neg(file, rdoff);
				//printf("file->natts = %llu\n", file->natts);

				/* Populate list of att objects */
				var->atts = get_cdf_atts_t(file, rdoff, var->natts);

			}

		} else {
			printf("ERROR - Expecting NC_DIMENSION, NC_VARIABLE or NC_ATTRIBUTE not %d.\n", nc_item);
		}

	} else {
		/* dim_list is ABSENT - Move rdoff by 4|8 bytes */
		if (file->fmt < 5) {
			*rdoff = *rdoff + 4;
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] ABSENT. Moving forward 4 bytes.\n", file->rank);
#endif
		} else {
			*rdoff = *rdoff + 8;
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] ABSENT. Moving forward 8 bytes.\n", file->rank);
#endif
		}
	}
	return;
}

#ifdef H5_VOL_HAVE_PARALLEL
/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_new_obj_mpio
 *
 * Purpose:     Create a new cdf vol object with MPI-IO
 *
 * Return:      Success:    Pointer to the new H5VL_cdf_t object
 *              Failure:    NULL
 *-------------------------------------------------------------------------
 */
static H5VL_cdf_t *
H5VL_cdf_new_obj_mpio(MPI_File fh, hid_t acc_tpl)
{
	H5VL_cdf_t *new_obj;
	new_obj = (H5VL_cdf_t *)calloc(1, sizeof(H5VL_cdf_t));
	new_obj->fh = fh;
	new_obj->fmt = 0;
	new_obj->numrecs = 0;
	new_obj->ndims = 0;
	new_obj->natts = 0;
	if (H5Pget_driver(acc_tpl) == H5FD_MPIO) {
		H5Pget_fapl_mpio(acc_tpl, &new_obj->comm, &new_obj->info);
#ifdef ENABLE_CDF_VERBOSE
		printf("NOTE: H5Pget_driver(acc_tpl) == H5FD_MPIO.\n");
#endif
	} else {
		new_obj->comm = MPI_COMM_WORLD;
		new_obj->info = MPI_INFO_NULL;
	}
	MPI_Comm_rank(new_obj->comm, &new_obj->rank);
	MPI_File_get_size( fh, &new_obj->size);
	if (new_obj->size > 0) {
		new_obj->cache = (char *)calloc(H5VL_CDF_CACHE_SIZE, sizeof(char));
	}
	new_obj->cache_offset = -1;
	return new_obj;
} /* end H5VL__cdf_new_obj() */

#else /* H5_VOL_HAVE_PARALLEL */

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_new_obj
 *
 * Purpose:     Create a new cdf vol object (NO MPIO)
 *
 * Return:      Success:    Pointer to the new H5VL_cdf_t object
 *              Failure:    NULL
 *-------------------------------------------------------------------------
 */
static H5VL_cdf_t *
H5VL_cdf_new_obj(int fd)
{
	H5VL_cdf_t *new_obj;
	new_obj = (H5VL_cdf_t *)calloc(1, sizeof(H5VL_cdf_t));
	new_obj->fd = fd;
	new_obj->fmt = 0;
	new_obj->numrecs = 0;
	new_obj->ndims = 0;
	new_obj->natts = 0;
	new_obj->rank = 0;
	return new_obj;
} /* end H5VL__cdf_new_obj() */

#endif /* H5_VOL_HAVE_PARALLEL */

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_free_obj
 *
 * Purpose:     Release a CDF object
 *
 * Return:      Success:    0
 *              Failure:    -1
 *-------------------------------------------------------------------------
 */
static herr_t
H5VL_cdf_free_obj(H5VL_cdf_t *obj)
{
	/* Free name strings for each dimension */
	for (int d=0; d<obj->ndims; d++){
		free(obj->dims[d].name->string);
		free(obj->dims[d].name);
	}

	/* Free values for each attribute */
	for (int d=0; d<obj->natts; d++){
		free(obj->atts[d].name->string);
		free(obj->atts[d].name);
		free(obj->atts[d].values);
	}

	/* Free values for each variable */
	for (int d=0; d<obj->nvars; d++){
		if (obj->vars[d].dimrank > 0){
			free(obj->vars[d].dimids);
			free(obj->vars[d].dimlens);
		}
		free(obj->vars[d].name->string);
		free(obj->vars[d].name);
		for (int a=0; a<obj->vars[d].natts; a++){
			free(obj->vars[d].atts[a].name->string);
			free(obj->vars[d].atts[a].name);
			free(obj->vars[d].atts[a].values);
		}
		free(obj->vars[d].atts);
	}

	/* Free list of dimensions */
	if (obj->ndims > 0)
		free(obj->dims);

	/* Free list of attributes */
	if (obj->natts > 0)
		free(obj->atts);

	/* Free list of attributes */
	if (obj->nvars > 0)
		free(obj->vars);

#ifdef H5_VOL_HAVE_PARALLEL
	/* Free cache */
	if (obj->size > 0) {
		free(obj->cache);
	}
#endif

	/* Free file name */
	free(obj->fname);

	/* Free vol object */
	free(obj);

	return 0;
} /* end H5VL_cdf_free_obj() */

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_read_header
 *
 * Purpose:     Read CDF-formatted file header & populate H5VL_cdf_t struct
 *              for a given file. Assume formatting according to:
 *              http://cucis.ece.northwestern.edu/projects/PnetCDF/CDF-5.html
 *
 * Return:      Success:    0
 *              Failure:    -1
 *-------------------------------------------------------------------------
 */
static herr_t
cdf_read_header(H5VL_cdf_t* file)
{
	herr_t retval=0;
	off_t rdoff = 0;

	/* Read "magic" part of header */
	get_cdf_spec(file, &rdoff);

	/* Read length of record dimension (numrecs) */
	file->numrecs = get_non_neg(file, &rdoff);
	/* Note: Above read will be "STREAMING" rather than "NON_NEG" if all bytes are \xFF */

	/* Read dimension list (dim_list) */
	get_header_list_item(file, &rdoff, NULL);

	/* Read attribute list (att_list) */
	get_header_list_item(file, &rdoff, NULL);

	/* Read variable list (var_list) */
	get_header_list_item(file, &rdoff, NULL);

	return retval;
}

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_register
 *
 * Purpose:     Register the VOL connector and retrieve an ID
 *              for it.
 *
 * Return:      Success:    The ID for the cdf VOL connector
 *              Failure:    -1
 *-------------------------------------------------------------------------
 */
hid_t
H5VL_cdf_register(void)
{
    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    /* Singleton register the pass-through VOL connector ID */
    if(H5I_VOL != H5Iget_type(H5VL_CDF_g))
        H5VL_CDF_g = H5VLregister_connector(&H5VL_cdf_g, H5P_DEFAULT);

    return H5VL_CDF_g;
} /* end H5VL_cdf_register() */


/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_init
 *
 * Purpose:     Initialize this VOL connector.
 *
 * Return:      Success:    0
 *              Failure:    -1
 *-------------------------------------------------------------------------
 */
static herr_t
H5VL_cdf_init(hid_t vipl_id)
{

#ifdef ENABLE_CDF_LOGGING
    printf("------- CDF VOL INIT\n");
#endif

    /* Shut compiler up about unused parameter */
    vipl_id = vipl_id;

    return 0;
} /* end H5VL_cdf_init() */


/*---------------------------------------------------------------------------
 * Function:    H5VL_cdf_term
 *
 * Purpose:     Terminate this VOL connector, performing any necessary
 *              operations for the connector that release connector-wide
 *              resources (usually created / initialized with the 'init'
 *              callback).
 *
 * Return:      Success:    0
 *              Failure:    (Can't fail)
 *---------------------------------------------------------------------------
 */
static herr_t
H5VL_cdf_term(void)
{

#ifdef ENABLE_CDF_LOGGING
    printf("------- CDF VOL TERM\n");
#endif

    /* Reset VOL ID */
    H5VL_CDF_g = H5I_INVALID_HID;

    return 0;
} /* end H5VL_cdf_term() */


/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_file_open
 *
 * Purpose:     Opens a file/container created with this connector
 *
 * Return:      Success:    Pointer to a file object
 *              Failure:    NULL
 *-------------------------------------------------------------------------
 */
static void *
H5VL_cdf_file_open(const char *name, unsigned flags, hid_t fapl_id,
    hid_t dxpl_id, void **req)
{
		H5VL_cdf_t *file;
		int err;

#ifdef H5_VOL_HAVE_PARALLEL

		MPI_File fh = 0;
		MPI_Offset fsize;

		/* Using MPIO VFD */
		err = MPI_File_open(MPI_COMM_WORLD, (char *)name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
		file = H5VL_cdf_new_obj_mpio( fh, fapl_id );
#ifdef ENABLE_CDF_LOGGING
		printf("------- CDF VOL FILE Open (MPIO) (size=%llu)\n",file->size);
#endif

#else /* H5_VOL_HAVE_PARALLEL */

		/* Not using MPIO - Use POSIX API instead */
		file = H5VL_cdf_new_obj( open(name, O_RDONLY) );
#ifdef ENABLE_CDF_LOGGING
		printf("------- CDF VOL FILE Open (POSIX)\n");
#endif

#endif /* H5_VOL_HAVE_PARALLEL */

		/* Try reading/parsing the header */
		if ( cdf_read_header(file) ) {
			cdf_handle_error( "Failed to read CDF header on rank 0.", __LINE__ );
		}

#ifdef ENABLE_CDF_VERBOSE
#ifdef H5_VOL_HAVE_PARALLEL
		MPI_Barrier(MPI_COMM_WORLD);
		if (file->rank == 0) {
#endif
			printf("----\n");
			printf("====\n");
			printf("#### Header is parsed.\n");
			printf("====\n");
			printf("----\n");
#ifdef H5_VOL_HAVE_PARALLEL
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

		/* Store this open file in a list of pointers */
		file->fname = (char *) malloc(strlen(name) * sizeof(char));
		strcpy(file->fname, name);
		file_objects_open[n_file_objects_open] = file;
		n_file_objects_open++;

		return (void *)file;
} /* end H5VL_cdf_file_open() */


/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_file_close
 *
 * Purpose:     Closes a file.
 *
 * Return:      Success:    0
 *              Failure:    -1, file not closed.
 *-------------------------------------------------------------------------
 */
static herr_t
H5VL_cdf_file_close(void *file, hid_t dxpl_id, void **req)
{
    H5VL_cdf_t *o = (H5VL_cdf_t *)file;
    herr_t ret_value;

#ifdef ENABLE_CDF_LOGGING
    printf("------- CDF VOL FILE Close\n");
#endif

#ifdef H5_VOL_HAVE_PARALLEL
    /* MPI_File Close */
		ret_value = (herr_t) MPI_File_close(&o->fh);
#else
    /* POSIX Close */
    ret_value = close(o->fd);
#endif

    H5VL_cdf_free_obj( o );

    return ret_value;
} /* end H5VL_cdf_file_close() */





/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
static void *
H5VL_cdf_group_open(void *obj, const H5VL_loc_params_t *loc_params, const char *varName,
                    hid_t gapl_id, hid_t dxpl_id, void **req) {
	H5VL_cdf_t *file;
	cdf_var_t *var;
	if (H5I_GROUP == loc_params->obj_type) {
		printf("H5VL_cdf_group_open with: H5I_GROUP\n");
		var = (cdf_var_t *) obj;
	} else if (H5I_FILE == loc_params->obj_type) {
		printf("H5VL_cdf_group_open with: H5I_FILE\n");
		file = (H5VL_cdf_t *) obj;
	}
	/*
	if ((varName != NULL) && (1 == strlen(varName)) && ('/' == varName[0])) { // h5ls request
		cdf_var_t *dummy = (cdf_var_t *)calloc(1, sizeof(cdf_var_t));
		//dummy->name = NULL;
		return dummy;
	}
	cdf_var_t *var = (H5VL_adios_var_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_var_t));
	var->name = (char *)SAFE_CALLOC(strlen(varName) + 1, sizeof(char));
	sprintf(var->name, "%s", varName);
	var->fileReader = bpFile;
	var->dimInfo = (H5VL_adios_dim_t *)SAFE_CALLOC(1, sizeof(H5VL_adios_dim_t));
	struct adios_index_comp_struct_v1 *varInAdios = adios_get_var_byname(var->fileReader, varName);
	var->ntimestep = adios_get_var_nsteps(varInAdios);
	*/
	return var;
}

static herr_t
H5VL_cdf_group_get(void *obj, H5VL_group_get_t get_type, hid_t dxpl_id,
                     void **req, va_list arguments) {
	herr_t ret_value = 0; /* Return value */
	switch (get_type) {
	case H5VL_GROUP_GET_INFO: {
		H5VL_cdf_t *file = (H5VL_cdf_t *)obj;
		H5VL_loc_params_t loc_params = va_arg(arguments, H5VL_loc_params_t);
		va_arg(arguments, H5G_info_t*); /* NOT SURE WHY THIS IS NEEDED TO WORK? */
		H5G_info_t *grp_info = va_arg(arguments, H5G_info_t*);
		if (loc_params.type == H5VL_OBJECT_BY_SELF) { /* H5Gget_info */
			//printf("(B) Inside H5VL_cdf_group_get -- grp_info->nlinks=%llu\n",grp_info->nlinks);
			grp_info->nlinks = (hsize_t)(file->nvars);
			//printf("(C) Inside H5VL_cdf_group_get -- grp_info->nlinks=%llu\n",grp_info->nlinks);
		}
		break;

// /* Information struct for group (for H5Gget_info/H5Gget_info_by_name/H5Gget_info_by_idx) */
// typedef struct H5G_info_t {
//     H5G_storage_type_t  storage_type;   /* Type of storage for links in group */
//     hsize_t     nlinks;                 /* Number of links in group */
//     int64_t     max_corder;             /* Current max. creation order value for group */
//     hbool_t     mounted;                /* Whether group has a file mounted on it */
// } H5G_info_t;

	}
	default:
		printf("ERROR: cdf_group_get() only supports get info\n");
		ret_value = -1;
	}
	return ret_value;
}
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */





/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_attr_open
 *
 * Purpose      Open attribute (For now, this must be BY NAME)
*
 * Return:      Pinter to attribute object
 *
 *-------------------------------------------------------------------------
 */
void *
H5VL_cdf_attr_open(void *obj, const H5VL_loc_params_t *loc_params, const char *attr_name,
                   hid_t aapl_id, hid_t dxpl_id, void **req)
{
	void *ret_value = NULL; /* Return value */
	H5VL_cdf_t *file;
	cdf_var_t * var;
	cdf_att_t *atts;
	cdf_non_neg_t natts;
	int rank;

	//printf("[%d] <%d> <%d> <<%d>> Looking for attribute with name = %s\n",file->rank, loc_params->type, H5VL_OBJECT_BY_NAME, loc_params->obj_type, attr_name);

	// H5I_FILE        = 1,        /* type ID for File objects                     */
	// H5I_GROUP,                  /* type ID for Group objects                    */
	// H5I_DATATYPE,               /* type ID for Datatype objects                 */
	// H5I_DATASPACE,              /* type ID for Dataspace objects                */
	// H5I_DATASET,                /* type ID for Dataset objects                  */
	// H5I_ATTR,                   /* type ID for Attribute objects                */

	if (loc_params->obj_type == H5I_FILE) {

		file = (H5VL_cdf_t *)obj;
		natts = file->natts;
		atts = file->atts;
		rank = file->rank;

	} else if (loc_params->obj_type == H5I_DATASET) {

		var = (cdf_var_t *)obj;
		natts = var->natts;
		atts = var->atts;
		rank = ((H5VL_cdf_t *)var->file)->rank;

	} else {
		printf("ERROR -- H5VL_cdf_attr_open loc_params->obj_type not supported yet.\n");
	}

	//if(loc_params->type == H5VL_OBJECT_BY_NAME) {
		for (int iatt=0; iatt<natts; iatt++ ){
			cdf_att_t *attr = &(atts[iatt]);
			if ( !strcmp( attr->name->string, attr_name) ) {
#ifdef ENABLE_CDF_VERBOSE
				printf("[%d] <%s> <%s> ATTRIBUTE MATCH!\n",rank, attr->name->string, attr_name);
#endif
				ret_value = (void *)attr;
				break;
			}
		}
	//} else {
	//	printf("ERROR -- H5VL_cdf_attr_open option not supported yet.\n");
	//}

	/* TODO: Need to add more useful error handling here (if name is not found etc) */

	return ret_value;
}

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_attr_read
 *
 * Purpose      Read the attribute using memcpy from cdf_att_t object.
*
 * Return:      Success:    Non-negative
 *              Failure:    Negative
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5VL_cdf_attr_read(void *obj, hid_t dtype_id, void *buf, hid_t dxpl_id, void **req)
{
	herr_t ret_value = 0;
	cdf_att_t *attr = (cdf_att_t *)obj;

	/* We already read this when we opened the file. Just do a memory copy */
	memcpy(buf, attr->values, (attr->nvals * attr->nc_type_size));

	return ret_value;
}

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_attr_close
 *
 * Purpose      Close attribute (not really used in this VOL)
*
 * Return:      Success:    Non-negative
 *              Failure:    Negative
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5VL_cdf_attr_close(void *obj, hid_t dxpl_id, void **req) {
	herr_t ret_value=0;
	obj = NULL;
	/* Don't need to free anything (taken care of when file is closed) */
	return ret_value;
}

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_dataset_open
 *
 * Purpose:     Opens a dataset
 *
 * Return:      Success:    Pointer to the new dataset
 *              Failure:    NULL
 *
 *-------------------------------------------------------------------------
 */
void *
H5VL_cdf_dataset_open(void *obj, const H5VL_loc_params_t *loc_params,
                      const char *name, hid_t dapl_id, hid_t dxpl_id, void **req)
{
	void *ret_value = NULL;             /* Return value */

	/* Assume an HDF5 dataset is the same as a CDF variable */
	/* Loop through variables, to find the desired object */
	H5VL_cdf_t *file = (H5VL_cdf_t *)obj;
	for (int ivar=0; ivar<file->nvars; ivar++ ){
		cdf_var_t *var = &(file->vars[ivar]);
		if ( !strcmp( var->name->string, name) ) {
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] <%s> <%s> MATCH! vsize = %llu, nc_type_size = %d, offset = %llu\n",file->rank, var->name->string, name, var->vsize, var->nc_type_size, var->offset);
#endif
			var->file = (void *)file; /* Need to set the file pointer */
			ret_value = (void *)var;
			break;
		}
	}

	/* TODO: Need to add more useful error handling here (if group name is not found etc) */

	return ret_value;
} /* end H5VL_cdf_dataset_open() */

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_dataset_get
 *
 * Purpose      Get specific information about the dataset
*
 * Return:      Success:    Non-negative
 *              Failure:    Negative
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5VL_cdf_dataset_get(void *obj, H5VL_dataset_get_t get_type, hid_t dxpl_id,
                     void **req, va_list arguments) {
  cdf_var_t *var = (cdf_var_t *)obj;
  herr_t ret_value = 0; /* Return value */

	switch (get_type) {

	/* H5Dget_space */
	case H5VL_DATASET_GET_SPACE: {
		hid_t *ret_id = va_arg(arguments, hid_t *);
		if (var->dimrank == 0) {
			*ret_id = H5Screate(H5S_SCALAR);
		} else {
			//// ADIOS VOL Code:
			//// if (var->is_global_dim >= 0)
			//// not sure how the local variables will behave here
			//// right now retrieve the datablock specified todimInfo in createVar()
			//*ret_id = H5Screate_simple(var->dimInfo->ndim, var->dimInfo->dims, NULL);
			//H5VL_cdf_t *file = (H5VL_cdf_t *) (var->file);
			*ret_id = H5Screate_simple(var->dimrank, var->dimlens, NULL);
		}
		break;
	}
	//case H5VL_DATASET_GET_TYPE: {
		//// ADIOS VOL Code:
		//struct adios_index_comp_struct_v1 *varInAdios = adios_get_var_byname(var->fileReader, var->name);
		//*ret_id = toHDF5type(varInAdios);
		//break;
	//}
	//case H5VL_DATASET_GET_DCPL: {
		//// ADIOS VOL Code:
		//*ret_id = H5Pcreate(H5P_DATASET_CREATE);
		//break;
	//}
	case H5VL_DATASET_GET_STORAGE_SIZE: {
		hsize_t *ret = va_arg(arguments, hsize_t *);
		*ret = (hsize_t)(var->vsize);
		//printf("H5VL_cdf_dataset_get -- var->vsize = %llu (*ret=%llu)\n",var->vsize,*ret);
		break;
	}
	default:
		printf("ERROR -- H5VL_DATASET_GET option not supported yet.\n");
		return -1;
	}

	return ret_value;
}

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_dataset_close
 *
 * Purpose      Close dataset (not really used in this VOL)
*
 * Return:      Success:    Non-negative
 *              Failure:    Negative
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5VL_cdf_dataset_close(void *dset, hid_t dxpl_id, void **req) {
	herr_t ret_value=0;
	dset = NULL;
	/* Don't need to free anything (taken care of when file is closed) */
	return ret_value;
}

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_dataset_read
 *
 * Purpose:     Reads data from dataset through the VOL
*
 * Return:      Success:    Non-negative
 *              Failure:    Negative
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5VL_cdf_dataset_read(void *dset, hid_t mem_type_id, hid_t mem_space_id,
                      hid_t file_space_id, hid_t plist_id, void *buf, void **req)
{

	herr_t ret_value=0;
	cdf_var_t *var = (cdf_var_t *)dset; /* Cast dset to cdf_var_t */
	H5VL_cdf_t *file = (H5VL_cdf_t *)(var->file); /* Get pointer to file object */
	char *charbuf = (char *)buf; /* Cast buffer to character for byte arithmetic */
	size_t dataTypeSize;
	H5S_sel_type memSelType, h5selType; /* Selection types in memory and file */
	hsize_t npoints, npoints_mem, nblocks_file, nblocks_mem;
	hsize_t *blockinfo_file;
	herr_t h5_status;
	int ndims, i, n;
	H5_posMap *sourceSelOrderInC = NULL;
	H5_posMap *targetSelOrderInC = NULL;
#ifdef H5_VOL_HAVE_PARALLEL
	H5FD_mpio_xfer_t xfer_mode = H5FD_MPIO_INDEPENDENT;
	MPI_Status mpi_status;
	MPI_Datatype structtype;
#endif

#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] <%s> In H5VL_cdf_dataset_read\n",file->rank, var->name->string);
#endif

	dataTypeSize = H5Tget_size(mem_type_id);
	memSelType = H5S_SEL_ALL;
	if (mem_space_id != 0) {
		memSelType = H5Sget_select_type(mem_space_id);
	}

	h5selType = H5S_SEL_ALL;
	if (file_space_id != 0) {
		h5selType = H5Sget_select_type(file_space_id);
	}

	if ((memSelType == H5S_SEL_NONE) || (h5selType == H5S_SEL_NONE)) {
		/* Nothing was selected, do nothing */
		return 0;
	}

#ifdef H5_VOL_HAVE_PARALLEL
	H5Pget_dxpl_mpio(plist_id, &xfer_mode);
#ifdef ENABLE_CDF_VERBOSE
	printf("[%d] <%s> xfer_mode %d\n", file->rank, var->name->string, xfer_mode);
#endif
#endif

	if (h5selType == H5S_SEL_ALL) {

#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] <%s> Reading %llu values for file offset %llu\n", file->rank, var->name->string, var->vsize, var->offset);
		printf("[%d] <%s> nc_type_size = %d\n",file->rank,var->name->string,var->nc_type_size);
#endif
#ifdef H5_VOL_HAVE_PARALLEL
		/* Use MPI-IO to read entire variable */
		if (xfer_mode == H5FD_MPIO_COLLECTIVE)
			MPI_File_read_at_all( file->fh, var->offset, charbuf, (var->vsize) , MPI_BYTE, &mpi_status );
		else
			MPI_File_read_at( file->fh, var->offset, charbuf, (var->vsize) , MPI_BYTE, &mpi_status );
#else
		/* Use POSIX to read entire variable */
		pread(file->fd, charbuf, var->vsize, (off_t)var->offset);
#endif

		/* Swap bytes for each value (need to add rigorous test for "when" to do this) */
		for (i=0; i<(var->vsize); i+=var->nc_type_size){
			bytestr_rev(&charbuf[i],var->nc_type_size);
		}

	} else if (h5selType == H5S_SEL_HYPERSLABS) {
		/* Generate arrays of flattened positions for each point in the selection.
		 * The H5_posMap type will have an index (posCompact), and a position in
		 * file or memory space (posInSource).  The position is the element index
		 * for a flattened representation of the selection.
		 */
		GetSelOrder(file_space_id, h5selType, &sourceSelOrderInC);
		GetSelOrder(mem_space_id, memSelType, &targetSelOrderInC);
		npoints = H5Sget_select_npoints(file_space_id);
		npoints_mem = H5Sget_select_npoints(mem_space_id);

#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] <%s> Hyperslab selection has %llu points.\n", file->rank, var->name->string, npoints);
		printf("[%d] <%s> Variable size is %llu.\n", file->rank, var->name->string, var->vsize);
#endif

		if ((var->vsize) == (npoints * var->nc_type_size)) {

			/* The hyperslab is actually the entire variable, read everything */
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] <%s> Reading %llu values for file offset %llu\n",file->rank,var->name->string,var->vsize, var->offset);
			printf("[%d] <%s> nc_type_size = %d\n",file->rank,var->name->string,var->nc_type_size);
#endif
#ifdef H5_VOL_HAVE_PARALLEL
			/* Use MPI-IO to read entire variable */
			MPI_File_read_at( file->fh, var->offset, charbuf, (var->vsize) , MPI_BYTE, &mpi_status );
#else
			/* Use POSIX to read entire variable */
			pread(file->fd, charbuf, var->vsize, (off_t)var->offset);
#endif
			/* Swap bytes for each value (need to add rigorous test for "when" to do this) */
			for (i=0; i<(var->vsize); i+=var->nc_type_size){
				bytestr_rev(&charbuf[i],var->nc_type_size);
			}

		} else {

			/* Hyperslab is a proper subset of elements - create block list */
			nblocks_file = H5Sget_select_hyper_nblocks(file_space_id);
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] <%s> nblocks_file = %llu\n", file->rank, var->name->string, nblocks_file);
#endif
			ndims = (var->dimrank);
			blockinfo_file = (hsize_t *) malloc( sizeof(hsize_t) * 2 * ndims * nblocks_file );
			h5_status = H5Sget_select_hyper_blocklist(file_space_id, (hsize_t)0, nblocks_file, blockinfo_file);

			/* Lets create a new buffer to actually read in the data */
			char *output_data = (char *) malloc( (var->nc_type_size) * npoints );

			/* Calculate the maximum number of contiguous chunks, and the chunk size */
			uint64_t max_chunks = nblocks_file;
			uint64_t chunk_len_s = 0;
			for (i = ndims-1; i >= 0; i--) {
				if (i < ndims-1) {
					max_chunks *= (blockinfo_file[0 + ndims + i] - blockinfo_file[0 + i] + 1);
				} else {
					chunk_len_s = (blockinfo_file[0 + ndims + i] - blockinfo_file[0 + i] + 1);
				}
			}

			/* Allocate arrays to store the offset of each chunk
			 * (Assume the selection will result in all chunks being the same length,
			 * with the magnitude being the size in the 0th dimension)
			 */
#ifdef H5_VOL_HAVE_PARALLEL
			MPI_Aint *chunk_off = (MPI_Aint *) malloc (sizeof(MPI_Aint) * max_chunks);
			MPI_Datatype *chunk_typ = (MPI_Datatype *) malloc (sizeof(MPI_Datatype) * max_chunks);
#else
			off_t *chunk_off = (off_t *) malloc (sizeof(off_t) * max_chunks);
#endif
			int *chunk_len = (int *) malloc (sizeof(int) * max_chunks);

			/* Populate chunk_off */
			int last_ind = sourceSelOrderInC[0].posInSource;
			int chunk_ind = 0;  /* number of contiguous chunks (current index) */
			int chunk_cnt = 1;  /* num elements in current chunk */
			chunk_off[chunk_ind] = last_ind * (var->nc_type_size);
			chunk_len[chunk_ind] = npoints * (var->nc_type_size);
#ifdef H5_VOL_HAVE_PARALLEL
			chunk_typ[chunk_ind] = MPI_BYTE;
#endif
			for (i = 1; i<npoints; i++) {
				int this_ind = sourceSelOrderInC[i].posInSource;
				if (this_ind != (last_ind+1)) {
					chunk_len[chunk_ind] = chunk_cnt * (var->nc_type_size);
					chunk_ind++;
					chunk_off[chunk_ind] = this_ind * (var->nc_type_size);
					chunk_len[chunk_ind] = (npoints - i) * (var->nc_type_size);
#ifdef H5_VOL_HAVE_PARALLEL
					chunk_typ[chunk_ind] = MPI_BYTE;
#endif
					chunk_cnt = 0;
				}
				chunk_cnt++;
				last_ind = this_ind;
			}
			max_chunks = (chunk_ind+1);

#ifdef H5_VOL_HAVE_PARALLEL
			MPI_Type_create_struct( (int)max_chunks, chunk_len, chunk_off, chunk_typ, &structtype );
			MPI_Type_commit( &structtype );
			MPI_File_set_view( file->fh, var->offset, MPI_BYTE, structtype, "native", MPI_INFO_NULL );
			if (xfer_mode == H5FD_MPIO_COLLECTIVE)
				MPI_File_read_all( file->fh, output_data, (npoints * var->nc_type_size), MPI_BYTE, &mpi_status );
			else
				MPI_File_read( file->fh, output_data, (npoints * var->nc_type_size), MPI_BYTE, &mpi_status );
			MPI_Type_free(&structtype);
#else
			/* POSIX-based read (Just read one chunk at a time) */
			uint64_t mem_ind = 0;
			for (i = 0; i<max_chunks; i++) {
				off_t total_offset = ((off_t)var->offset) + chunk_off[i];
				pread(file->fd, &output_data[mem_ind], chunk_len[i], total_offset);
				mem_ind += chunk_len[i];
			}
#endif
			/* Swap bytes for each value (need to add rigorous test for "when" to do this) */
			for (i=0; i<(npoints * var->nc_type_size); i+=var->nc_type_size){
				bytestr_rev(&output_data[i],var->nc_type_size);
			}

			/* Assign data to proper memory-space selection */
			assignToMemSpace(sourceSelOrderInC, targetSelOrderInC, npoints, var->nc_type_size, output_data, charbuf);
			//DEBUG (direct copy): memcpy(charbuf, output_data, (npoints * var->nc_type_size));

			free(chunk_len);
			free(chunk_off);
#ifdef H5_VOL_HAVE_PARALLEL
			free(chunk_typ);
#endif
			free(output_data);
			free(blockinfo_file);
		}
		free(sourceSelOrderInC);
		free(targetSelOrderInC);

	} else {
		printf("ERROR!!! Only H5S_SEL_ALL available for now.\n");
	}


	return ret_value;
}

/*-------------------------------------------------------------------------
 * Function:    cdf_vol_var_get_attname
 *
 * Purpose      Populates "attname" with attribute name at specified
 *              index (iatts)
 *
 * Return:      Success:    Non-negative
 *              Failure:    Negative
 *-------------------------------------------------------------------------
 */
herr_t
cdf_vol_var_get_attname(const char *file_name, const char *var_name, char *att_name, int iatts) {
	herr_t ret_value=-1;
	int i,j;
	for (i=0; i<n_file_objects_open; i++) {
		H5VL_cdf_t *obj = file_objects_open[i];
		if ( !strcmp( obj->fname, file_name) ) {
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] <%s> <%s> FILE MATCH!\n",obj->rank, obj->fname, file_name);
#endif
			for (j=0; j<obj->nvars; j++) {
				cdf_var_t *var = &(obj->vars[j]);
				if ( !strcmp( var->name->string, var_name) ) {
#ifdef ENABLE_CDF_VERBOSE
					printf("[%d] <%s> <%s> VAR MATCH!\n",obj->rank, var->name->string, var_name);
#endif
					cdf_att_t attr = var->atts[iatts];
					char *name = attr.name->string;
					strcpy(att_name, name);
					ret_value = 0;
					break;
				}
			}
			break;
		}
	}
	return ret_value;
}

/*-------------------------------------------------------------------------
 * Function:    cdf_vol_var_get_natts
 *
 * Purpose      Populates "natts" with number of attributes.
 *
 * Return:      Success:    Non-negative
 *              Failure:    Negative
 *-------------------------------------------------------------------------
 */
herr_t
cdf_vol_var_get_natts(const char *file_name, const char *var_name, int *natts) {
	herr_t ret_value=-1;
	int i,j;
	for (i=0; i<n_file_objects_open; i++) {
		H5VL_cdf_t *obj = file_objects_open[i];
		if ( !strcmp( obj->fname, file_name) ) {
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] <%s> <%s> FILE MATCH!\n",obj->rank, obj->fname, file_name);
#endif
			for (j=0; j<obj->nvars; j++) {
				cdf_var_t *var = &(obj->vars[j]);
				if ( !strcmp( var->name->string, var_name) ) {
#ifdef ENABLE_CDF_VERBOSE
					printf("[%d] <%s> <%s> VAR MATCH! (var->natts: %llu)\n",obj->rank, var->name->string, var_name, var->natts);
#endif
					*natts = var->natts;
					ret_value = 0;
					break;
				}
			}
			break;
		}
	}
	return ret_value;
}

 /*-------------------------------------------------------------------------
  * Function:    cdf_vol_file_get_iname
  *
  * Purpose      Populates "iname" with name of variable or attribute at given
	*              index ("item").
  *
  * Return:      Success:    Non-negative
  *              Failure:    Negative
  *-------------------------------------------------------------------------
  */
herr_t
cdf_vol_file_get_iname(const char *file_name, char *iname, int item, cdf_vol_get_t type) {
	herr_t ret_value=-1;
	char *name;
	int i,j;
	for (i=0; i<n_file_objects_open; i++) {
		H5VL_cdf_t *obj = file_objects_open[i];
		if ( !strcmp( obj->fname, file_name) ) {
			if (type == CDF_VOL_VAR) {
				cdf_var_t var =  obj->vars[item];
				name = var.name->string;
			} else {
				cdf_att_t attr =  obj->atts[item];
				name = attr.name->string;
			}
			//iname = (char *) malloc( strlen(name) * sizeof(char) );
			strcpy(iname, name);
			ret_value = 0;
			break;
		}
	}
	return ret_value;
}

/*-------------------------------------------------------------------------
 * Function:    cdf_vol_file_get_nitems
 *
 * Purpose      Populates "nitems" with number of variables of attributes,
 *              depending on "type" parameter.
 *
 * Return:      Success:    Non-negative
 *              Failure:    Negative
 *-------------------------------------------------------------------------
 */
herr_t
cdf_vol_file_get_nitems(const char *file_name, int *nitems, cdf_vol_get_t type) {
	herr_t ret_value=-1;
	int i,j;
	for (i=0; i<n_file_objects_open; i++) {
		H5VL_cdf_t *obj = file_objects_open[i];
		if ( !strcmp( obj->fname, file_name) ) {
#ifdef ENABLE_CDF_VERBOSE
			printf("[%d] <%s> <%s> NAME MATCH!\n",obj->rank, obj->fname, file_name);
#endif
			if (type == CDF_VOL_VAR)
				*nitems = obj->nvars;
			else
				*nitems = obj->natts;
			ret_value = 0;
			break;
		}
	}
	return ret_value;
}
