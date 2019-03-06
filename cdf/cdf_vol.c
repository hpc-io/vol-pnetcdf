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
/* #define ENABLE_CDF_LOGGING */

/* Hack for missing va_copy() in old Visual Studio editions
 * (from H5win2_defs.h - used on VS2012 and earlier)
 */
#if defined(_WIN32) && defined(_MSC_VER) && (_MSC_VER < 1800)
#define va_copy(D,S)      ((D) = (S))
#endif

/************/
/* Typedefs */
/************/

/* The CDF VOL info object */
typedef struct H5VL_cdf_t {
	int  fd; /* posix file object */
	uint8_t fmt; /* CDF File Spec (1,2,5) */
	uint64_t numrecs; /* length of record dimension */
	uint64_t ndims; /* number of dimensions */
} H5VL_cdf_t;

/********************* */
/* Function prototypes */
/********************* */

/* Helper routines */
static H5VL_cdf_t *H5VL_cdf_new_obj(int fd);
static herr_t H5VL_cdf_free_obj(H5VL_cdf_t *obj);
static uint32_t bytestr_to_uint32_t(char *bytestr);
static uint64_t bytestr_to_uint64_t(char *bytestr);
static uint32_t get_uint32_t(H5VL_cdf_t* file, off_t *rdoff);
static uint64_t get_uint64_t(H5VL_cdf_t* file, off_t *rdoff);
static uint64_t get_non_neg(H5VL_cdf_t* file, off_t *rdoff);
static uint8_t get_cdf_spec(H5VL_cdf_t* file, off_t *rdoff);

/* "Management" callbacks */
static herr_t H5VL_cdf_init(hid_t vipl_id);
static herr_t H5VL_cdf_term(void);

/* Attribute callbacks */

/* Dataset callbacks */

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
        NULL, //H5VL_cdf_dataset_open,                      /* open */
        NULL, //H5VL_cdf_dataset_read,                      /* read */
        NULL, //H5VL_cdf_dataset_write,                     /* write */
        NULL, //H5VL_cdf_dataset_get,                       /* get */
        NULL, //H5VL_cdf_dataset_specific,                  /* specific */
        NULL, //H5VL_cdf_dataset_optional,                  /* optional */
        NULL, //H5VL_cdf_dataset_close                      /* close */
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


static uint32_t
bytestr_to_uint32_t(char *bytestr)
{
	/* Note: Little Endian Case */
	//uint32_t myInt1 = (uint32_t)bytestr[0] + ((uint32_t)bytestr[1] << 8) + ((uint32_t)bytestr[2] << 16) + ((uint32_t)bytestr[3] << 24);
	/* Big Endian Case */
	uint32_t myInt1 = ((uint32_t)bytestr[0] << 24) + ((uint32_t)bytestr[1] << 16) + ((uint32_t)bytestr[2] << 8) + (uint32_t)bytestr[3];
	/* Print Result */
	//printf("uint32_t value = %d (bytes = %d %d %d %d) \n", myInt1,(uint8_t)bytestr[0],(uint8_t)bytestr[1],(uint8_t)bytestr[2],(uint8_t)bytestr[3] );
	return myInt1;
}


static uint64_t
bytestr_to_uint64_t(char *bytestr)
{
	/* Big Endian Case */
	uint64_t myInt1 = ((uint64_t)bytestr[0] << 56) + ((uint64_t)bytestr[1] << 48) + ((uint64_t)bytestr[2] << 40) + ((uint64_t)bytestr[3] << 32) +  ((uint64_t)bytestr[4] << 24) + ((uint64_t)bytestr[5] << 16) + ((uint64_t)bytestr[6] << 8) + (uint64_t)bytestr[7];
	/* Print Result */
	//printf("uint32_t value = %llu (bytes = %d %d %d %d %d %d %d %d) \n", myInt1,(uint8_t)bytestr[0],(uint8_t)bytestr[1],(uint8_t)bytestr[2],(uint8_t)bytestr[3],(uint8_t)bytestr[4],(uint8_t)bytestr[5],(uint8_t)bytestr[6],(uint8_t)bytestr[7] );
	return myInt1;
}

static uint32_t
get_uint32_t(H5VL_cdf_t* file, off_t *rdoff)
{
	char bytestr[4];
	pread(file->fd, &bytestr[0], 4, *rdoff); *rdoff = *rdoff + 4;
	return (uint64_t)bytestr_to_uint32_t( bytestr );
}

static uint64_t
get_uint64_t(H5VL_cdf_t* file, off_t *rdoff)
{
	char bytestr[8];
	pread(file->fd, &bytestr[0], 8, *rdoff); *rdoff = *rdoff + 8;
	return bytestr_to_uint64_t( bytestr );
}

static uint8_t
get_cdf_spec(H5VL_cdf_t* file, off_t *rdoff)
{
	char bytestr[4];
	pread(file->fd, &bytestr[0], 4, *rdoff); *rdoff = *rdoff + 4;
	file->fmt = (uint8_t) bytestr[3];
	printf("------- (FILE FORMAT = %c%c%c-%d) \n", bytestr[0], bytestr[1], bytestr[2], (uint8_t) bytestr[3]);
	return (uint8_t) bytestr[3];
}

static uint64_t
get_non_neg(H5VL_cdf_t* file, off_t *rdoff)
{
	/* Read length of record dimension (numrecs) */
	if (file->fmt < 5) {
		/* NON_NEG = <non-negative INT> = <32-bit signed integer, Bigendian, two's complement> */
		return (uint64_t) get_uint32_t(file,rdoff);
	} else {
		/* NON_NEG = <non-negative INT64> = <64-bit signed integer, Bigendian, two's complement> */
		return get_uint64_t(file,rdoff);
	}
}

/*-------------------------------------------------------------------------
 * Function:    H5VL_cdf_new_obj
 *
 * Purpose:     Create a new cdf vol object
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

    return new_obj;
} /* end H5VL__cdf_new_obj() */


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
    free(obj);
    return 0;
} /* end H5VL_cdf_free_obj() */


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

    printf("------- CDF VOL INIT\n");

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

    printf("------- CDF VOL TERM\n");

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
		off_t rdoff = 0;
		char bytes_sp[4];
		char bytes_dp[8];

		printf("------- CDF VOL FILE Open\n");

		/* For now - Just using POSIX to prototype/test basic VOL */
		file = H5VL_cdf_new_obj( open(name, O_RDONLY) );

		/* Read "magic" part of header */
		get_cdf_spec(file, &rdoff);

		/* Read length of record dimension (numrecs) */
		file->numrecs = (uint64_t) get_non_neg(file, &rdoff);
		/* Note: Above read will be "STREAMING" rather than "NON_NEG" if all bytes are \xFF */

		/* Read dimension list (dim_list) */

		/* First check if first 4 bytes are 0 */
		uint32_t nc_dim = get_uint32_t(file, &rdoff);
		if ( nc_dim > (uint32_t)0 ){

			/* Make Sure these 4 bytes are "NC_DIMENSION" */
			if ( nc_dim != (uint32_t)10 ) {
				printf("ERROR - Expecting NC_DIMENSION (%d).\n", nc_dim);
			} else {
				//printf("GOOD - NC_DIMENSION (%d).\n", nc_dim);
			}

			/* Get nelms of dim_list (file->ndims) */
			file->ndims = (uint64_t) get_non_neg(file, &rdoff);
			printf("ndims = %llu\n", file->ndims);

		} else {
			/* dim_list is ABSENT - Move rdoff by 4|8 bytes */
			if (file->fmt < 5)
				rdoff = rdoff + 4;
			else
				rdoff = rdoff + 8;
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

    printf("------- CDF VOL FILE Close\n");

    /* For now - Just using POSIX to prototype/test basic VOL */
    ret_value = close(o->fd);
    H5VL_cdf_free_obj( o );

    return ret_value;
} /* end H5VL_cdf_file_close() */


