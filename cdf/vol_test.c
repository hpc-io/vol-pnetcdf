#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "hdf5.h"
#include "cdf_vol.h"
#define FILE "testfile.h5"
#include <pnetcdf.h>


static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d of %s: %s\n", lineno, __FILE__, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}


int write_cdf( ) {
    int ret, ncid=0, nprocs, rank, dimid, varid1=0, varid2=0, ndims=1;
    char filename[256], buf[13] = "Hello World\n";
    int *data=NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    strcpy(filename, "testfile.nc");

    if (rank == 0) {
        ret = ncmpi_create(MPI_COMM_SELF, filename,
                           NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncid);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_def_dim(ncid, "d1", nprocs, &dimid);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_def_var(ncid, "v1", NC_INT, ndims, &dimid, &varid1);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_def_var(ncid, "v2", NC_INT, ndims, &dimid, &varid2);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_put_att_text(ncid, NC_GLOBAL, "string", 13, buf);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_enddef(ncid);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        /* first reason this approach is not scalable:  need to allocate
        * enough memory to hold data from all processors */
        data = (int*)calloc(nprocs, sizeof(int));
    }

    /* second reason this approach is not scalable: sending to rank 0
     * introduces a serialization point, even if using an optimized
     * collective routine */
    MPI_Gather(&rank, 1, MPI_INT, data, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        /* and lastly, the third reason this approach is not scalable: I/O
         * happens from a single processor.  This approach can be ok if the
         * amount of data is quite small, but almost always the underlying
         * MPI-IO library can do a better job */
        MPI_Offset start[1], count[1];
        start[0]=0, count[0]=nprocs;
        ret = ncmpi_put_vara_int_all(ncid, varid1, start, count, data);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_put_vara_int_all(ncid, varid2, start, count, data);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        ret = ncmpi_close(ncid);
        if (ret != NC_NOERR) handle_error(ret, __LINE__);

        free(data);
    }

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


	MPI_Init(&argc, &argv);

	if(argc !=2){
		printf("Input: connector name, e.g., cdf\n");
		return 0;
	}
	strcpy(connector_name,argv[1]);


	/* Write File with pNetCDF */
	write_cdf();

	fapl = H5Pcreate (H5P_FILE_ACCESS);

	vol_id = H5VLregister_connector (&H5VL_cdf_g, fapl);
	assert(vol_id > 0);
	assert(H5VLis_connector_registered(connector_name) == 1);

	acc_tpl = H5Pcreate (H5P_FILE_ACCESS);
	H5Pset_vol(acc_tpl, vol_id, &fapl);

	file_id = H5Fopen(FILE, H5F_ACC_RDWR, acc_tpl);

	H5Fclose(file_id);
	H5Pclose(acc_tpl);
	H5Pclose(fapl);
	H5VLterminate(vol_id);
	H5VLunregister_connector(vol_id);
	assert(H5VLis_connector_registered(connector_name) == 0);

	MPI_Finalize();
	return 0;
}
