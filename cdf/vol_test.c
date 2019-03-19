#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "hdf5.h"
#include "cdf_vol.h"
#define FILE "testfile.nc"
#include <pnetcdf.h>

/*
 *    PNETCDF code for writing the file is taken from the github repository:
 *    https://github.com/Parallel-NetCDF/PnetCDF/blob/master/examples/C/collective_write.c
 */

/* Set SIMPLE_TEST=0 for multi-dimensional tests */
#define SIMPLE_TEST 0
#define NDIMS       3
#define NUM_VARS   10
#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}
static int verbose=1;

static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d of %s: %s\n", lineno, __FILE__, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

static
void print_info(MPI_Info *info_used)
{
    int     flag;
    char    info_cb_nodes[64], info_cb_buffer_size[64];
    char    info_striping_factor[64], info_striping_unit[64];

    strcpy(info_cb_nodes,        "undefined");
    strcpy(info_cb_buffer_size,  "undefined");
    strcpy(info_striping_factor, "undefined");
    strcpy(info_striping_unit,   "undefined");

    MPI_Info_get(*info_used, "cb_nodes", 64, info_cb_nodes, &flag);
    MPI_Info_get(*info_used, "cb_buffer_size", 64, info_cb_buffer_size, &flag);
    MPI_Info_get(*info_used, "striping_factor", 64, info_striping_factor, &flag);
    MPI_Info_get(*info_used, "striping_unit", 64, info_striping_unit, &flag);

    printf("MPI hint: cb_nodes        = %s\n", info_cb_nodes);
    printf("MPI hint: cb_buffer_size  = %s\n", info_cb_buffer_size);
    printf("MPI hint: striping_factor = %s\n", info_striping_factor);
    printf("MPI hint: striping_unit   = %s\n", info_striping_unit);
}

int write_cdf_col(MPI_Comm comm, char *filename, int cmode, int len) {

	char str[512];
	int i, j, rank, nprocs, ncid, bufsize, err, nerrs=0;
	int *buf[NUM_VARS], psizes[NDIMS], dimids[NDIMS], varids[NUM_VARS];
	double write_timing, max_write_timing, write_bw;
	MPI_Offset gsizes[NDIMS], starts[NDIMS], counts[NDIMS];
	MPI_Offset write_size, sum_write_size;
	MPI_Info info_used;

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &nprocs);

	for (i=0; i<NDIMS; i++)
	    psizes[i] = 0;

	MPI_Dims_create(nprocs, NDIMS, psizes);
	starts[0] = rank % psizes[0];
	starts[1] = (rank / psizes[1]) % psizes[1];
	starts[2] = (rank / (psizes[0] * psizes[1])) % psizes[2];
	//printf("[%d] starts = %llu %llu %llu\n", rank, starts[0], starts[1], starts[2]);
	//printf("[%d] psizes = %d %d %d\n", rank, psizes[0], psizes[1], psizes[2]);

	bufsize = 1;
	for (i=0; i<NDIMS; i++) {
	    gsizes[i] = (MPI_Offset)len * psizes[i];
	    starts[i] *= len;
	    counts[i]  = len;
	    bufsize   *= len;
	}

	/* allocate buffer and initialize with non-zero numbers */
	for (i=0; i<NUM_VARS; i++) {
	    buf[i] = (int *) malloc(bufsize * sizeof(int));
	    for (j=0; j<bufsize; j++) buf[i][j] = rank * i + 123 + j;
	}

	MPI_Barrier(comm);
	write_timing = MPI_Wtime();

	/* create the file */
	cmode |= NC_CLOBBER;
	err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid);
	if (err != NC_NOERR) {
	    printf("Error at %s:%d ncmpi_create() file %s (%s)\n",
	           __FILE__,__LINE__,filename,ncmpi_strerror(err));
	    MPI_Abort(comm, -1);
	    exit(1);
	}

	/* define dimensions */
	for (i=0; i<NDIMS; i++) {
	    sprintf(str, "%c", 'x'+i);
	    err = ncmpi_def_dim(ncid, str, gsizes[i], &dimids[i]); ERR
	}

	/* define variables */
	for (i=0; i<NUM_VARS; i++) {
	    sprintf(str, "var%d", i);
	    err = ncmpi_def_var(ncid, str, NC_INT, NDIMS, dimids, &varids[i]); ERR
	}

	/* exit the define mode */
	err = ncmpi_enddef(ncid); ERR

	/* get all the hints used */
	err = ncmpi_inq_file_info(ncid, &info_used); ERR

	/* write one variable at a time */
	for (i=0; i<NUM_VARS; i++) {
	    err = ncmpi_put_vara_int_all(ncid, varids[i], starts, counts, buf[i]);
	    ERR
	}

	/* close the file */
	err = ncmpi_close(ncid); ERR

	write_timing = MPI_Wtime() - write_timing;

	write_size = bufsize * NUM_VARS * sizeof(int);
	for (i=0; i<NUM_VARS; i++) free(buf[i]);

	MPI_Reduce(&write_size, &sum_write_size, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
	MPI_Reduce(&write_timing, &max_write_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

	if (rank == 0 && verbose) {
	    float subarray_size = (float)bufsize*sizeof(int)/1048576.0;
	    print_info(&info_used);
	    printf("Local array size %d x %d x %d integers, size = %.2f MB\n",len,len,len,subarray_size);
	    sum_write_size /= 1048576.0;
	    printf("Global array size %lld x %lld x %lld integers, write size = %.2f GB\n",
	           gsizes[0], gsizes[1], gsizes[2], sum_write_size/1024.0);

	    write_bw = sum_write_size/max_write_timing;
	    printf(" procs    Global array size  exec(sec)  write(MB/s)\n");
	    printf("-------  ------------------  ---------  -----------\n");
	    printf(" %4d    %4lld x %4lld x %4lld %8.2f  %10.2f\n\n", nprocs,
	           gsizes[0], gsizes[1], gsizes[2], max_write_timing, write_bw);
	}
	MPI_Info_free(&info_used);

return nerrs;

}

int write_cdf( ) {
    int ret, ncfile, nprocs, rank, dimid1, varid1, varid2, ndims=1;
    char filename[256], buf[13] = "Hello World\n";
		MPI_Offset start, count=1;
    int data;
		double data_dbl;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    strcpy(filename, FILE);

    ret = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_dim(ncfile, "d1", nprocs, &dimid1);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_var(ncfile, "v1", NC_INT, ndims, &dimid1, &varid1);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_var(ncfile, "v2", NC_DOUBLE, ndims, &dimid1, &varid2);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_put_att_text(ncfile, NC_GLOBAL, "string", 13, buf);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    /* All processors defined the dimensions, attributes, and variables,
     * but here in ncmpi_enddef is the one place where metadata I/O
     * happens.  Behind the scenes, rank 0 takes the information and writes
     * the netcdf header.  All processes communicate to ensure they have
     * the same (cached) view of the dataset */

    ret = ncmpi_enddef(ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    start=rank, count=1, data=(rank+10);
		data_dbl=((double)rank + 10.0);

    /* in this simple example every process writes its rank to two 1d variables */
    ret = ncmpi_put_vara_int_all(ncfile, varid1, &start, &count, &data);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_put_vara_double_all(ncfile, varid2, &start, &count, &data_dbl);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_close(ncfile);
		if (ret != NC_NOERR) handle_error(ret, __LINE__);

    return 0;
}

int read_cdf( ) {

    int i, j, rank, nprocs, ret;
    int ncfile, ndims, nvars, ngatts, unlimited;
    int var_ndims, var_natts;;
    MPI_Offset *dim_sizes, var_size;
    MPI_Offset *start, *count;
    char filename[256], varname[NC_MAX_NAME+1];
    int *dimids=NULL;
    nc_type type;
    int *data=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    strcpy(filename, FILE);

    ret = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    /* reader knows nothing about dataset, but we can interrogate with query
     * routines: ncmpi_inq tells us how many of each kind of "thing"
     * (dimension, variable, attribute) we will find in the file  */

    /* no communication needed after ncmpi_open: all processors have a cached
     * view of the metadata once ncmpi_open returns */

    ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    /* we do not really need the name of the dimension or the variable for
     * reading in this example.  we could, in a different example, take the
     * name of a variable on the command line and read just that one */

    dim_sizes = (MPI_Offset*) calloc(ndims, sizeof(MPI_Offset));
    /* netcdf dimension identifiers are allocated sequentially starting
     * at zero; same for variable identifiers */
    for(i=0; i<ndims; i++)  {
        ret = ncmpi_inq_dimlen(ncfile, i, &(dim_sizes[i]) );
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
    }

    for(i=0; i<nvars; i++) {
        /* obtain the number of dimensions of variable i, so we can allocate
         * the dimids array */
        ret = ncmpi_inq_varndims(ncfile, i, &var_ndims);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);
        dimids = (int*) malloc(var_ndims * sizeof(int));

        /* much less coordination in this case compared to rank 0 doing all
         * the i/o: everyone already has the necessary information */
        ret = ncmpi_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                &var_natts);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        start = (MPI_Offset*) calloc(var_ndims, sizeof(MPI_Offset));
        count = (MPI_Offset*) calloc(var_ndims, sizeof(MPI_Offset));

        /* we will simply decompose along one dimension.  Generally the
         * application has some algorithm for domain decomposition.  Note
         * that data decomposition can have an impact on i/o performance.
         * Often it's best just to do what is natural for the application,
         * but something to consider if performance is not what was
         * expected/desired */

        start[0] = (dim_sizes[dimids[0]]/nprocs)*rank;
        count[0] = (dim_sizes[dimids[0]]/nprocs);
        var_size = count[0];

        for (j=1; j<var_ndims; j++) {
            start[j] = 0;
            count[j] = dim_sizes[dimids[j]];
            var_size *= count[j];
        }

        switch(type) {
            case NC_INT:
                data = (int*) calloc(var_size, sizeof(int));
                ret = ncmpi_get_vara_all(ncfile, i, start, count, data,
                        var_size, MPI_INT);
                free(data);
                if (ret != NC_NOERR) handle_error(ret, __LINE__);
                break;
            default:
                /* we can do this for all the known netcdf types but this
                 * example is already getting too long  */
                fprintf(stderr, "unsupported NetCDF type \n");
        }

        free(start);
        free(count);
        free(dimids);
    }

    ret = ncmpi_close(ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);
    free(dim_sizes);

    return 0;
}

int main(int argc, char **argv) {

	hid_t file_id, int_dataset_id, dbl_dataset_id, dataspace_id;  /* identifiers */
	hsize_t dims[2];
	herr_t status;
	char connector_name[25];
	hid_t acc_tpl;
	hid_t fapl;
	hid_t vol_id;
	char name[25];
	ssize_t len;
	hsize_t bytecnt;
	int *dset_data_int;
	double *dset_data_dbl;
	int nprocs, rank, i;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	//if(argc !=2){
	//	printf("Input: connector name, e.g., cdf\n");
	//	return 0;
	//}
	//strcpy(connector_name,argv[1]);
	strcpy(connector_name,"cdf");

	if(SIMPLE_TEST) {

		/* Write File with pNetCDF */
		write_cdf();

		/* Read File with pNetCDF */
		//read_cdf();

	} else {

		/* Write multi-dimensional dataset */
		int dlen = 8;
		int cmode = 0;
		write_cdf_col(MPI_COMM_WORLD, FILE, cmode, dlen);

	}

	/* Open File etc */
	fapl = H5Pcreate (H5P_FILE_ACCESS);
	vol_id = H5VLregister_connector (&H5VL_cdf_g, fapl);
	assert(vol_id > 0);
	assert(H5VLis_connector_registered(connector_name) == 1);
	acc_tpl = H5Pcreate (H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(acc_tpl, MPI_COMM_WORLD, MPI_INFO_NULL);
	H5Pset_vol(acc_tpl, vol_id, &fapl);
	file_id = H5Fopen(FILE, H5F_ACC_RDWR, acc_tpl);

	if(SIMPLE_TEST) {

		/* Read simple integer dataset */
		int_dataset_id = H5Dopen(file_id, "v1", H5P_DEFAULT);
		bytecnt = H5Dget_storage_size(int_dataset_id);
		dset_data_int = (int *) malloc( bytecnt );
		status = H5Dread(int_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data_int);
		printf("dset_data_int[%d] = %d\n", rank, dset_data_int[rank]);
		H5Dclose(int_dataset_id);
		free(dset_data_int);

		/* Read simple double-precision dataset */
		dbl_dataset_id = H5Dopen(file_id, "v2", H5P_DEFAULT);
		bytecnt = H5Dget_storage_size(dbl_dataset_id);
		dset_data_dbl = (double *) malloc( bytecnt );
		status = H5Dread(dbl_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data_dbl);
		printf("dset_data_dbl[%d] = %f\n", rank, dset_data_dbl[rank]);
		H5Dclose(dbl_dataset_id);
		free(dset_data_dbl);

	} else {

		/* Read 3D integer dataset (H5S_ALL) */
		int_dataset_id = H5Dopen(file_id, "var0", H5P_DEFAULT);
		bytecnt = H5Dget_storage_size(int_dataset_id);
		dset_data_int = (int *) malloc( bytecnt );
		status = H5Dread(int_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data_int);
		printf("dset_data_int[%d] = %d\n", rank, dset_data_int[rank]);
		H5Dclose(int_dataset_id);
		free(dset_data_int);




		/* Read simple hyperslab selection from var1 */
		int_dataset_id = H5Dopen(file_id, "var1", H5P_DEFAULT);
		dset_data_int = (int *) malloc( H5Dget_storage_size(int_dataset_id) );
		hid_t dataspace = H5Dget_space (int_dataset_id);    /* dataspace handle */
		int dimrank      = H5Sget_simple_extent_ndims (dataspace);
		hsize_t *dims_out = (hsize_t *) malloc( dimrank * sizeof(hsize_t) );
		int status_n  = H5Sget_simple_extent_dims (dataspace, dims_out, NULL);
		printf("\nDim-Rank: %d -- Dimensions: %lu x %lu x %lu\n", dimrank, (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]), (unsigned long)(dims_out[2]));
		int psizes[NDIMS];
		MPI_Offset starts[NDIMS];
		for (i=0; i<NDIMS; i++)
			psizes[i] = 0;
		MPI_Dims_create(nprocs, NDIMS, psizes);
		starts[0] = rank % psizes[0];
		starts[1] = (rank / psizes[1]) % psizes[1];
		starts[2] = (rank / (psizes[0] * psizes[1])) % psizes[2];
		printf("[%d] starts = %llu %llu %llu\n", rank, starts[0], starts[1], starts[2]);
		printf("[%d] psizes = %d %d %d\n", rank, psizes[0], psizes[1], psizes[2]);
		/*
		 * Define hyperslab in the dataset.
		 */
		hsize_t     stride[NDIMS];
		hsize_t     block[NDIMS];
		hsize_t     count[NDIMS];              /* size of the hyperslab in the file */
		hsize_t     offset[NDIMS];             /* hyperslab offset in the file */
		for (i=0; i<NDIMS; i++) {
			block[i] = dims_out[i] / psizes[i];
			stride[i] = block[i];
			count[i] = 1;
			offset[i] = starts[i] * block[i];
		}
		herr_t status_h5 = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, stride, count, block);
		printf("[%d] offset = %llu %llu %llu\n", rank, offset[0], offset[1], offset[2]);
		printf("[%d] stride = %llu %llu %llu\n", rank, stride[0], stride[1], stride[2]);
		printf("[%d] count = %llu %llu %llu\n", rank, count[0], count[1], count[2]);
		printf("[%d] block = %llu %llu %llu\n", rank, block[0], block[1], block[2]);

		/*
		 * Define the memory dataspace.
		 */
		hsize_t dimsm[NDIMS];
		for (i=0; i<NDIMS; i++) {
			dimsm[i] = count[i];
		}
		hid_t memspace = H5Screate_simple (NDIMS, dimsm, NULL);

		/*
		 * Define memory hyperslab.
		 */
		hsize_t     stride_out[NDIMS];
		hsize_t     block_out[NDIMS];
		hsize_t     count_out[NDIMS];              /* size of the hyperslab in the file */
		hsize_t     offset_out[NDIMS];             /* hyperslab offset in the file */
		int hypersize=1;
		for (i=0; i<NDIMS; i++) {
			block_out[i] = block[i];
			stride_out[i] = block[i];
			count_out[i] = 1;
			offset_out[i] = 0;
			hypersize *= block[i];
		}
		status_h5 = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, stride_out, count_out, block_out);

		/*
		 * Read data from hyperslab in the file into the hyperslab in
		 * memory and display.
		 */
		int* data_out = (int *) malloc( hypersize * sizeof(int) );
		printf("[%d] data_out size = %d integers.\n", rank, hypersize);
		status_h5 = H5Dread (int_dataset_id, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, data_out);
		printf("data_out[0] = %d\n", data_out[0]);
		H5Dclose(int_dataset_id);
		free(data_out);
		free(dims_out);
		free(dset_data_int);
	}

	/* Close File etc */
	H5Fclose(file_id);
	H5Pclose(acc_tpl);
	H5Pclose(fapl);
	H5VLterminate(vol_id);
	H5VLunregister_connector(vol_id);
	assert(H5VLis_connector_registered(connector_name) == 0);

	MPI_Finalize();
	return 0;
}
