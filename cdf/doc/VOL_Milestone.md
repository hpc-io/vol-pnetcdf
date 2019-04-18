# HDF5 VOL Connector for CDF Files


>**Richard Zamora, Venkatram Vishwanath & Paul Coffman**

>*Argonne National Laboratory*

>April 18th 2019


--

### Abstract

*The HDF5 technology suite offers a rich data model and flexible API for high-performance I/O.  The library interface can be used to greatly simplify the task of performing parallel I/O on large multi-dimensional arrays.  For example, hyper-slab selections make it relatively straightforward to access non-contiguous data in parallel.  Since the API is particularly well-suited for multi-dimensional data, many of it's features are also well-suited for accessing data stored in the array-based [Common Data Format (CDF)](https://cdf.gsfc.nasa.gov/) [1], which is popular in the climate-modeling community.  For this work, we have implemented a new Virtual Object Layer (VOL) connector, enabling HDF5 to read files written by CDF-based libraries, such as [PnetCDF](http://cucis.ece.northwestern.edu/projects/PnetCDF/) [2]. The CDF VOL prototype allows users to read data from CDF-formatted files, without the need for additional CDF-based libraries,  using the HDF5 API.  In this document, we will discuss the details and performance of the CDF VOL connector prototype.*

--

### Introduction

The [Common Data Format (CDF)](https://cdf.gsfc.nasa.gov/) [1] is a self-describing data format for the storage of both scalar and multidimensional data.  Like HDF5 [4], CDF aims to provide users with both performance and inter-platform portability. [Parallel netCDF (PnetCDF)](http://cucis.ece.northwestern.edu/projects/PnetCDF/) [2] is a parallel I/O library for accessing files stored in CDF-1, 2, and 5 formats.  It was originally introduced in 2003 as an alternative to the existing [Network Common Data Form (NetCDF)](https://www.unidata.ucar.edu/software/netcdf/) library from Unidata [3], which only allowed users to read and write CDF-formatted data in serial at the time.  Since the introduction of versions 4.0 and 4.1 of Unidata's library, netCDF has supported parallel I/O through both HDF5 and PnetCDF backends, respectively. This means that netCDF-4 is able to read and write files in either HDF5 or CDF format.

The goal of this work is to enable read-only interoperability between the HDF5 library and CDF-formatted data.  To this end, we have implemented a Virtual Object Layer (VOL) connector, allowing user applications to read CDF data using the HDF5 API.  In this document, we will present the implementation details of the CDF VOL prototype, and discuss opportunities to improve its functionality and performance in future iterations.  We also provide instructions for building and running a simple benchmark example.

### The Virtual Object Layer in HDF5

To date, the most recent [development branch of the HDF5 library repository](https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5/browse) includes a Virtual Object Layer (VOL) API [7]. As illustrated in **Fig 1**, the VOL is essentially a new software layer that allows all HDF5 API calls to be redirected through a specific VOL *Connector* (connectors are labeled as "plugins" in the figure). The goal is to provide an application with the HDF5 data model and API, but allow different underlying storage mechanisms. For all work discussed in this document, we used the VOL interface available in the development branch of the HDF5 library (cloned March 1st, 2019; commit b23079de3af).  

![](vol-image.png)

**Fig 1** *Illustration of the overall HDF5 library architecture with the VOL interface. All API calls are directed to a specific VOL connector (or "plugin" in the figure).*


### CDF Schema Mapping

The CDF file format is relatively well-suited for the HDF5 VOL interface, because both data models were designed with multidimensional arrays in mind. In CDF, these arrays are called *variables*, while in HDF5 they are called *datasets*. What makes the mapping even more convenient is the fact that the CDF data model is essentially a subset of HDF5.  The CDF data schema consists of two fundamental components (variables and attributes), while HDF5 consists of three (groups, datasets, and attributes).  Since CDF does not support any notion of hierarchical groups, the most straightforward strategy is to treat CDF variables as HDF5 datasets, and to treat all CDF files as if there is only one group (the **root** group). In CDF, attributes can either be global (attached to the file), or variable-local (attached to a variable). In HDF5 attributes can be attached to both 

![](cdf-layout.png)

**Fig AAA** *General layout of a CDF-formatted file. The formatting calls for three general sections: (1) A file header containing all necessary metadata, (2) The non-record variable data, and (3) the record variable data.*

### Implementation Details

The CDF VOL-Connector prototype discussed here is available in the public [`pnetcdf_vol`](https://xgitlab.cels.anl.gov/ExaHDF5/pnetcdf_vol) repository, hosted on xGitLab [10]. Most of the VOL-connector source code can be found in `/cdf/cdf_vol.c`, where the `H5VL_cdf_g` connector (of type `H5VL_class_t`) is defined as follows:

```
/* CDF VOL connector class struct */
const H5VL_class_t H5VL_cdf_g = {
    H5VL_CDF_VERSION,                      /* version      */
    (H5VL_class_value_t)H5VL_CDF_VALUE,    /* value        */
    H5VL_CDF_NAME,                         /* name         */
    0,                                     /* capability flags */
    H5VL_cdf_init,                         /* initialize   */
    H5VL_cdf_term,                         /* terminate    */
    sizeof(H5VL_cdf_info_t),               /* info size    */
...
    {                                           /* attribute_cls */
...
        H5VL_cdf_attr_open,                         /* open */
        H5VL_cdf_attr_read,                         /* read */
...
        H5VL_cdf_attr_close                         /* close */
    },
    {                                           /* dataset_cls */
...
        H5VL_cdf_dataset_open,                      /* open */
        H5VL_cdf_dataset_read,                      /* read */
...
        H5VL_cdf_dataset_get,                       /* get */
...
        H5VL_cdf_dataset_close                      /* close */
    },
...
    {                                           /* file_cls */
...
        H5VL_cdf_file_open,                         /* open */
...
        H5VL_cdf_file_close                         /* close */
    },
    {                                           /* group_cls */
...
        H5VL_cdf_group_open,                        /* open */
        H5VL_cdf_group_get,                         /* get */
...
    },
...
};

```

In the code snippet shown above, all lines with `...` correspond to API callback functions that are currently undefined for the CDF connector (meaning the vales are set to `NULL`).  Since the purpose of the CDF prototype is to provide a minimal COL implementation for reading CDF-formatted files, a majority of the available callback functions are non-critical (and some are even non-applicable).

In addition to the thirteen callback functions included within the `H5VL_cdf_g` definition, we also add a number of helper/utility functions, as well as four additional application-accessible API functions (`cdf_vol_var_get_attname`, `cdf_vol_var_get_natts`, `cdf_vol_file_get_nitems`, and `cdf_vol_file_get_iname`).  The additional API functions should be viewed as a *temporary* interface that should eventually be replaced by a more-complete `H5VL_cdf_g` implementation.

#### Parsing the File Header

The typical procedure for using the CDF VOL connector to read a CDF-formatted file within a user application (once the connector is registered and added to the file-access property list) starts with the opening of the file (using `H5Dopen`).  When passed through the HDF5 VOL, the `H5Dopen` call is redirected to the `H5VL_cdf_file_open` function, which ultimately reads through the file header to populate a `H5VL_cdf_t` structure for the open file. Throughout the CDV VOL implementation, the following structures are used to organize the metadata that is read-in from the file header:

- `H5VL_cdf_t`: Files
- `cdf_var_t`: Variables
- `cdf_att_t`: Attributes
- `cdf_dim_t`: Dimensions
- `cdf_name_t`: Names
- `cdf_offset_t`: Offsets
- `cdf_non_neg_t`: Non-negative values

Once the header is completely parsed, the `H5VL_cdf_t` structure includes all the necessary information for reading variables and/or attributes from the file. The structure is defines as:

```
/* The CDF VOL "file" object */
typedef struct H5VL_cdf_t {
	char *fname; /* Store the file name */
	uint8_t fmt; /* CDF File Spec (1,2,5) */
	cdf_non_neg_t numrecs; /* length of record dimension */
	cdf_non_neg_t record_stride; /* Stride between records of same variable */
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
```

The creation/population of the `H5VL_cdf_t` object is performed by the `H5VL_cdf_new_obj` (or `H5VL_cdf_new_obj_mpio`) function, which uses a collection of utility functions to parse through any file header that is in compliance with the CDF-1, CDF-2, or CDF-5 specification [9].  In order to provide a VOL-connector implementation with both parallel and serial support, the header can be parsed in parallel (using MPI-IO) or in serial (using POSIX).  For the current implementation, the MPI-based implementation does not attempt to actually discretize the header.  Instead, `H5VL_CDF_CACHE_SIZE` bytes are read in by the 0th rank and broadcast to other ranks for redundant/parallel processing. this procedure is followed iteratively, until the entire header is parsed.

#### Reading Multidimensional Variables

![](record-layout.png)

**Fig BBB** *Layout of non-record variables (top) and record variables (bottom) in a CDF-formatted file. For non-record variables, the entire variable is stored contiguously in the file, using row-major mapping of array elements.  For record variables, each record (containing one or more variables) is stored contiguously.*

### Performance



![](theta-results.png)

**Fig CCC** *Performance results for the CDF VOL prototype, using the `vol_test` benchmark on ALCF Theta. The results correspond to the average timing of 10 trials, with each trial using: 2MB per process (2 1MB variables), 1-8 nodes, 32 ppn, and 52 Lustre stripes (8MB stripe size).*

### Building & Running the CDF VOL Benchmark

--

#### Missing Features (TODO)

- Custom API calls were added to allow a user application to query the count and names of attributes and/or variables (`cdf_vol_var_get_attname`, `cdf_vol_var_get_natts`, `cdf_vol_file_get_nitems`, and `cdf_vol_file_get_iname`).  In the future, the same functionality should be implemented within the official VOL interface.
- The `H5VL_DATASET_GET_TYPE` case needs to be added to the `H5VL_cdf_dataset_get` implementation. Otherwise, the user must know the datatype of the variable a priori.
- Reading of hyper-slab selections can be further optimized for cases with simple/regular patterns.  For example, it may not be necessary to explicitly map the file-space selection to the memory-space selection if both cases are contiguous (and in the same order).
- Need to test (and possibly debug) the code at scale (128+ nodes). Preliminary tests at scale (128+ nodes) show problematic behavior.

#### References

[1] CDF (NASA):
https://cdf.gsfc.nasa.gov/

[2] PnetCDF (Argonne & Northwestern):
http://cucis.ece.northwestern.edu/projects/PnetCDF/

[3] NetCDF (Unidata):
https://www.unidata.ucar.edu/software/netcdf/

[4] HDF5:
https://www.hdfgroup.org/solutions/hdf5/

[5] VOL Presentation by HDF Group:
https://slideplayer.com/slide/7607091/

[6] VOL User Guide:
https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5doc/browse/RFCs/HDF5/VOL/user_guide

[7] HDF5 Source Code. 
https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5/browse

[8] ADIOS VOL Example:
https://bitbucket.org/berkeleylab/exahdf5/src/master/vol_plugins/swift/

[9] CDF Specification:
http://cucis.ece.northwestern.edu/projects/PnetCDF/CDF-5.html

[10] CDF VOL Repository:
https://xgitlab.cels.anl.gov/ExaHDF5/pnetcdf_vol