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
#define ENABLE_CDF_LOGGING
//#define ENABLE_CDF_VERBOSE

/* Hack for missing va_copy() in old Visual Studio editions
 * (from H5win2_defs.h - used on VS2012 and earlier)
 */
#if defined(_WIN32) && defined(_MSC_VER) && (_MSC_VER < 1800)
#define va_copy(D,S)      ((D) = (S))
#endif

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
	cdf_non_neg_t natts; /* number of attributes for this var */
	cdf_att_t *atts; /* list of cdf_att_t objects for this var */
} cdf_var_t;

/* The CDF VOL info object */
typedef struct H5VL_cdf_t {
	uint8_t fmt; /* CDF File Spec (1,2,5) */
	cdf_non_neg_t numrecs; /* length of record dimension */
	cdf_non_neg_t ndims; /* number of dimensions */
	cdf_non_neg_t natts; /* number of attributes */
	cdf_non_neg_t nvars; /* number of variables */
	cdf_dim_t *dims; /* list of cdf_dim_t objects */
	cdf_att_t *atts; /* list of cdf_att_t objects */
	cdf_var_t *vars; /* list of cdf_var_t objects */
#ifdef H5_HAVE_PARALLEL
	MPI_File fh; /* MPI File Handle */
	MPI_Comm comm; /* MPI Communicator */
	MPI_Comm info; /* MPI Info */
	MPI_Comm rank; /* MPI Rank */
	MPI_Offset size; /* Total Size of read-only file */
	char *cache; /* In-memory cache for parsing header (size = H5VL_CDF_CACHE_SIZE) */
	MPI_Offset cache_offset; /* What is the file offset for first byte of cache */
#else
	int fd; /* posix file object */
#endif
} H5VL_cdf_t;


/********************* */
/* Function prototypes */
/********************* */

/* Helper routines */
static void cdf_handle_error(char* errmsg, int lineno);
#ifdef H5_HAVE_PARALLEL
static H5VL_cdf_t *H5VL_cdf_new_obj_mpio(MPI_File fh, hid_t acc_tpl);
#else
static H5VL_cdf_t *H5VL_cdf_new_obj(int fd);
#endif
static herr_t H5VL_cdf_free_obj(H5VL_cdf_t *obj);
static void bytestr_rev(char *p, int size);
static uint32_t bytestr_to_uint32_t(char *bytestr);
static uint64_t bytestr_to_uint64_t(char *bytestr);
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
        NULL, //H5VL_cdf_attr_open,                         /* open */
        NULL, //H5VL_cdf_attr_read,                         /* read */
        NULL, //H5VL_cdf_attr_write,                        /* write */
        NULL, //H5VL_cdf_attr_get,                          /* get */
        NULL, //H5VL_cdf_attr_specific,                     /* specific */
        NULL, //H5VL_cdf_attr_optional,                     /* optional */
        NULL, //H5VL_cdf_attr_close                         /* close */
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
        NULL, //H5VL_cdf_group_open,                        /* open */
        NULL, //H5VL_cdf_group_get,                         /* get */
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

#ifdef H5_HAVE_PARALLEL
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
			printf("Initializing cache to size = %llu\n", rdsize);
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
bytestr_to_uint32_t(char *bytestr)
{
	/* Little Endian Case */
	//uint32_t myInt1 = (uint32_t)bytestr[0] + ((uint32_t)bytestr[1] << 8) + ((uint32_t)bytestr[2] << 16) + ((uint32_t)bytestr[3] << 24);
	/* Big Endian Case */
	uint32_t myInt1 = ((uint32_t)bytestr[0] << 24) + ((uint32_t)bytestr[1] << 16) + ((uint32_t)bytestr[2] << 8) + (uint32_t)bytestr[3];
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
bytestr_to_uint64_t(char *bytestr)
{
	/* Big Endian Case */
	uint64_t myInt1 = ((uint64_t)bytestr[0] << 56) + ((uint64_t)bytestr[1] << 48) + ((uint64_t)bytestr[2] << 40) + ((uint64_t)bytestr[3] << 32) +  ((uint64_t)bytestr[4] << 24) + ((uint64_t)bytestr[5] << 16) + ((uint64_t)bytestr[6] << 8) + (uint64_t)bytestr[7];
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
#ifdef H5_HAVE_PARALLEL
		cdf_cache_read(file, rdoff, &bytestr[0], 4);
#else
		pread(file->fd, &bytestr[0], 4, *rdoff);
#endif
	*rdoff = *rdoff + 4;
#ifdef ENABLE_CDF_VERBOSE
	printf(" bytes: %d %d %d %d\n", (uint8_t)bytestr[0], (uint8_t)bytestr[1], (uint8_t)bytestr[2], (uint8_t)bytestr[3]);
#endif
	return (uint64_t)bytestr_to_uint32_t( bytestr );
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
#ifdef H5_HAVE_PARALLEL
		cdf_cache_read(file, rdoff, &bytestr[0], 8);
#else
		pread(file->fd, &bytestr[0], 8, *rdoff);
#endif
	*rdoff = *rdoff + 8;
	return bytestr_to_uint64_t( bytestr );
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
#ifdef H5_HAVE_PARALLEL
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
#ifdef H5_HAVE_PARALLEL
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
#ifdef H5_HAVE_PARALLEL
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

			/* Get dimid list */
			for (int id=0; id<vars[ivar].dimrank; id++) {
				vars[ivar].dimids[id] = get_non_neg(file, rdoff);
#ifdef ENABLE_CDF_VERBOSE
				printf("[%d] vars[%d].dimids[%d] = %llu\n", file->rank, ivar, id, vars[ivar].dimids[id]);
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

#ifdef H5_HAVE_PARALLEL
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

#else /* H5_HAVE_PARALLEL */

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
	return new_obj;
} /* end H5VL__cdf_new_obj() */

#endif /* H5_HAVE_PARALLEL */

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
		if (obj->vars[d].dimrank > 0)
			free(obj->vars[d].dimids);
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

#ifdef H5_HAVE_PARALLEL
	/* Free cache */
	if (obj->size > 0) {
		free(obj->cache);
	}
#endif

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
		H5VL_cdf_info_t *info;
		H5VL_cdf_t *file;
		int err;

#ifdef H5_HAVE_PARALLEL

		MPI_File fh = 0;
		MPI_Offset fsize;

		/* Using MPIO VFD */
		err = MPI_File_open(MPI_COMM_WORLD, (char *)name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
		file = H5VL_cdf_new_obj_mpio( fh, fapl_id );
#ifdef ENABLE_CDF_LOGGING
		printf("------- CDF VOL FILE Open (MPIO) (size=%llu)\n",file->size);
#endif

#else /* H5_HAVE_PARALLEL */

		/* Not using MPIO - Use POSIX API instead */
		file = H5VL_cdf_new_obj( open(name, O_RDONLY) );
#ifdef ENABLE_CDF_LOGGING
		printf("------- CDF VOL FILE Open (POSIX)\n");
#endif

#endif /* H5_HAVE_PARALLEL */

		/* Try reading/parsing the header */
		if ( cdf_read_header(file) ) {
			cdf_handle_error( "Failed to read CDF header on rank 0.", __LINE__ );
		}

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

#ifdef H5_HAVE_PARALLEL
    /* MPI_File Close */
		ret_value = (herr_t) MPI_File_close(&o->fh);
#else
    /* POSIX Close */
    ret_value = close(o->fd);
#endif

    H5VL_cdf_free_obj( o );

    return ret_value;
} /* end H5VL_cdf_file_close() */


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

	/* Need to add more useful error handling here (if group name is not found etc) */

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

	// /* CDF variable object (structure) */
	// typedef struct cdf_var_t {
	// 	cdf_name_t *name;
	// 	void *file; /* "void" Pointer back to parent file object */
	// 	cdf_non_neg_t offset; /* Offset byte in file for this variable (begin) */
	// 	cdf_non_neg_t vsize; /* vsize */
	// 	cdf_nc_type_t nc_type; /* var data type */
	// 	int nc_type_size;
	// 	cdf_non_neg_t dimrank; /* Dimensionality (rank) of this variable */
	// 	cdf_non_neg_t *dimids; /* Dimension ID (index into dim_list) for variable
	//                           * shape. We say this is a "record variable" if and only
	//                           * if the first dimension is the record dimension.  */
	// 	cdf_non_neg_t natts; /* number of attributes for this var */
	// 	cdf_att_t *atts; /* list of cdf_att_t objects for this var */
	// } cdf_var_t;

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
		}
		break;
	}
	case H5VL_DATASET_GET_TYPE: {
		//// ADIOS VOL Code:
		//struct adios_index_comp_struct_v1 *varInAdios = adios_get_var_byname(var->fileReader, var->name);
		//*ret_id = toHDF5type(varInAdios);
		break;
	}
	case H5VL_DATASET_GET_DCPL: {
		//// ADIOS VOL Code:
		//*ret_id = H5Pcreate(H5P_DATASET_CREATE);
		break;
	}
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
 * Purpose:	Reads data from dataset through the VOL
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
	MPI_Status status;
	int i;

#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] <%s> In H5VL_cdf_dataset_read\n",file->rank, var->name->string);
#endif

	size_t dataTypeSize = H5Tget_size(mem_type_id);

	H5S_sel_type memSelType = H5S_SEL_ALL;
	if (mem_space_id != 0) {
		memSelType = H5Sget_select_type(mem_space_id);
	}

	H5S_sel_type h5selType = H5S_SEL_ALL;
	if (file_space_id != 0)
		h5selType = H5Sget_select_type(file_space_id);

	if ((memSelType == H5S_SEL_NONE) || (h5selType == H5S_SEL_NONE)) {
		/* Nothing was selected, do nothing */
		return 0;
	}

	if (h5selType == H5S_SEL_ALL) {

#ifdef ENABLE_CDF_VERBOSE
		printf("[%d] <%s> Reading %llu values for offset %llu\n",file->rank,var->name->string,var->vsize, var->offset);
		printf("[%d] <%s> nc_type_size = %d\n",file->rank,var->name->string,var->nc_type_size);
#endif
		MPI_File_read_at( file->fh, var->offset, (char *)buf, (var->vsize) , MPI_BYTE, &status );

		/* Swap bytes for each value (need to add rigorous test for "when" to do this) */
		for (i=0; i<(var->vsize); i+=var->nc_type_size){
			bytestr_rev(&buf[i],var->nc_type_size);
		}

	} else {
		printf("ERROR!!! Only H5S_SEL_ALL available for now.\n");
	}


	return ret_value;
}
