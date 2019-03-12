#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "hdf5.h"
#include "cdf_vol.h"
#define FILE "testfile.nc"
#include <pnetcdf.h>


static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d of %s: %s\n", lineno, __FILE__, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}


int write_cdf( ) {
    int ret, ncfile, nprocs, rank, dimid, varid1, varid2, ndims=1;
    char filename[256], buf[13] = "Hello World\n";
		MPI_Offset start, count=1;
    int data;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    strcpy(filename, FILE);

    ret = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_dim(ncfile, "d1", nprocs, &dimid);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_var(ncfile, "v1", NC_INT, ndims, &dimid, &varid1);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_var(ncfile, "v2", NC_INT, ndims, &dimid, &varid2);
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

    start=rank, count=1, data=rank;

    /* in this simple example every process writes its rank to two 1d variables */
    ret = ncmpi_put_vara_int_all(ncfile, varid1, &start, &count, &data);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_put_vara_int_all(ncfile, varid2, &start, &count, &data);
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

	hid_t file_id, dataset_id, dataspace_id;  /* identifiers */
	hsize_t dims[2];
	herr_t status;
	char connector_name[25];
	hid_t acc_tpl;
	hid_t fapl;
	hid_t vol_id;
	char name[25];
	ssize_t len;
	int dset_data[1024];

	MPI_Init(&argc, &argv);

	if(argc !=2){
		printf("Input: connector name, e.g., cdf\n");
		return 0;
	}
	strcpy(connector_name,argv[1]);


	/* Write File with pNetCDF */
	write_cdf();

	/* Read File with pNetCDF */
	read_cdf();

	fapl = H5Pcreate (H5P_FILE_ACCESS);

	vol_id = H5VLregister_connector (&H5VL_cdf_g, fapl);
	assert(vol_id > 0);
	assert(H5VLis_connector_registered(connector_name) == 1);

	acc_tpl = H5Pcreate (H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(acc_tpl, MPI_COMM_WORLD, MPI_INFO_NULL);
	H5Pset_vol(acc_tpl, vol_id, &fapl);

	file_id = H5Fopen(FILE, H5F_ACC_RDWR, acc_tpl);

	dataset_id = H5Dopen(file_id, "v1", H5P_DEFAULT);

	status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dset_data[0]);
	printf("dset_data[0] = %d\n",dset_data[0]);

	H5Fclose(file_id);
	H5Pclose(acc_tpl);
	H5Pclose(fapl);
	H5VLterminate(vol_id);
	H5VLunregister_connector(vol_id);
	assert(H5VLis_connector_registered(connector_name) == 0);

	MPI_Finalize();
	return 0;
}
