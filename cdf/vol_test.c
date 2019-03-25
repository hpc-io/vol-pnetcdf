#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "hdf5.h"
#include "cdf_vol.h"
#define FILE "testfile.nc"
#include <pnetcdf.h>

/*    HDF5: Read-Only CDF VOL Connector Test/Benchmark.
 *
 *    Author: R.J. Zamora
 *
 *    Use:  make clean; make
 *          rm *.h5 *.nc; mpirun -n 4 ./vol_test.ex --dimlen 128 --col
 *
 *    NOTE: PNETCDF code for writing the file is from the github repository:
 *    https://github.com/Parallel-NetCDF/PnetCDF/blob/master/examples/C/collective_write.c
 */

#define NDIMS          3
#define DIMLEN         8
#define NUM_VARS      10

static void handle_error(int status, int lineno)
{
	if(status!=NC_NOERR) {
		fprintf(stderr, "Error at line %d of %s: %s\n", lineno, __FILE__, ncmpi_strerror(status));
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

static void print_info(MPI_Info *info_used)
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
	char attrbuf[13] = "Hello World\n";
	char varattrbuf[13] = "Hello Var01\n";

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
		for (j=0; j<bufsize; j++) buf[i][j] = (rank+1)*10; // * i + 123 + j;
	}

	MPI_Barrier(comm);
	write_timing = MPI_Wtime();

	/* create the file */
	cmode |= NC_CLOBBER;
	err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid);
	if (err != NC_NOERR) {
		printf("Error at %s:%d ncmpi_create() file %s (%s)\n", __FILE__,__LINE__,filename,ncmpi_strerror(err));
		MPI_Abort(comm, -1);
		exit(1);
	}

	/* define dimensions */
	for (i=0; i<NDIMS; i++) {
		sprintf(str, "%c", 'x'+i);
		err = ncmpi_def_dim(ncid, str, gsizes[i], &dimids[i]);
		handle_error(err, __LINE__);
	}

	/* define variables */
	for (i=0; i<NUM_VARS; i++) {
		sprintf(str, "var%d", i);
		err = ncmpi_def_var(ncid, str, NC_INT, NDIMS, dimids, &varids[i]);
		handle_error(err, __LINE__);
	}

	/* Write attribute (Attached to file) */
	err = ncmpi_put_att_text(ncid, NC_GLOBAL, "string", 13, attrbuf);
	handle_error(err, __LINE__);

	/* Write attribute (Attached to variable #1) */
	err = ncmpi_put_att_text(ncid, varids[1], "varstring", 13, varattrbuf);
	handle_error(err, __LINE__);

	/* exit the define mode */
	err = ncmpi_enddef(ncid); handle_error(err, __LINE__);

	/* get all the hints used */
	err = ncmpi_inq_file_info(ncid, &info_used); handle_error(err, __LINE__);

	/* write one variable at a time */
	for (i=0; i<NUM_VARS; i++) {
		err = ncmpi_put_vara_int_all(ncid, varids[i], starts, counts, buf[i]);
		handle_error(err, __LINE__);
	}

	/* close the file */
	err = ncmpi_close(ncid);
	handle_error(err, __LINE__);

	write_timing = MPI_Wtime() - write_timing;

	write_size = bufsize * NUM_VARS * sizeof(int);
	for (i=0; i<NUM_VARS; i++) free(buf[i]);

	MPI_Reduce(&write_size, &sum_write_size, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
	MPI_Reduce(&write_timing, &max_write_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

	if (rank == 0) {
		float subarray_size = (float)bufsize*sizeof(int)/1048576.0;
		print_info(&info_used);
		printf("Local array size %d x %d x %d integers, size = %.2f MB\n",len,len,len,subarray_size);
		sum_write_size /= 1048576.0;
		printf("Global array size %lld x %lld x %lld integers, write size = %.2f GB\n", gsizes[0], gsizes[1], gsizes[2], sum_write_size/1024.0);
		write_bw = sum_write_size/max_write_timing;
		printf(" procs    Global array size  exec(sec)  write(MB/s)\n");
		printf("-------  ------------------  ---------  -----------\n");
		printf(" %4d    %4lld x %4lld x %4lld %8.2f  %10.2f\n\n", nprocs, gsizes[0], gsizes[1], gsizes[2], max_write_timing, write_bw);
	}
	MPI_Info_free(&info_used);

	return nerrs;
}

int main(int argc, char **argv) {

	hid_t file_id, int_dataset_id, dbl_dataset_id, dataspace_id, dxpl_plist_id;  /* identifiers */
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
	int dlen = DIMLEN;
	int cmode = 0;
	int use_collective = 0;
	double open_time=0.0, read_time_all=0.0, read_time_hyper=0.0;
	double max_open_time=0.0, max_read_time_all=0.0, max_read_time_hyper=0.0;
	long read_size_all=1.0, read_size_hyper=1.0;
	long sum_read_size_all=1.0, sum_read_size_hyper=1.0;
	double size_rpt;

	/* Parse Inpu Args */
	for (i=1;i<argc;i++) {
		if (strcmp(argv[i],"--dimlen") == 0) {
			i++; dlen = atoi(argv[i]);
		} else if (strcmp(argv[i],"--cmode") == 0) {
			i++; cmode = atoi(argv[i]);
		} else if (strcmp(argv[i],"--col") == 0) {
			use_collective = 1;
		} else {
			printf("ERROR - unrecognized parameter: %s.  Exitting.\n",argv[i]);
			exit(-1);
		}
	}
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	strcpy(connector_name,"cdf");

	/* Write multi-dimensional dataset */
	write_cdf_col(MPI_COMM_WORLD, FILE, cmode, dlen);

	/* Open File etc */
	fapl = H5Pcreate (H5P_FILE_ACCESS);
	vol_id = H5VLregister_connector (&H5VL_cdf_g, fapl);
	assert(vol_id > 0);
	assert(H5VLis_connector_registered(connector_name) == 1);
	acc_tpl = H5Pcreate (H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(acc_tpl, MPI_COMM_WORLD, MPI_INFO_NULL);
	H5Pset_vol(acc_tpl, vol_id, &fapl);
	open_time = MPI_Wtime();
	file_id = H5Fopen(FILE, H5F_ACC_RDWR, acc_tpl);
	open_time = MPI_Wtime()-open_time;
	MPI_Reduce(&open_time, &max_open_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	/* Set transfer property list */
	dxpl_plist_id = H5Pcreate(H5P_DATASET_XFER);
	if (use_collective)
		H5Pset_dxpl_mpio(dxpl_plist_id, H5FD_MPIO_COLLECTIVE);
	else
		H5Pset_dxpl_mpio(dxpl_plist_id, H5FD_MPIO_INDEPENDENT);

	/* Read 3D integer dataset (H5S_ALL) */
	int_dataset_id = H5Dopen(file_id, "var0", H5P_DEFAULT);
	bytecnt = H5Dget_storage_size(int_dataset_id);
	dset_data_int = (int *) malloc( bytecnt );
	read_time_all = MPI_Wtime();
	status = H5Dread(int_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, dxpl_plist_id, dset_data_int);
	read_time_all = MPI_Wtime() - read_time_all;
	read_size_all = bytecnt;
	MPI_Reduce(&read_size_all, &sum_read_size_all, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&read_time_all, &max_read_time_all, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	//printf("[%d] data_out size = %llu integers.\n", rank, bytecnt/sizeof(int) );
	//for (i=0; i<dlen*NDIMS; i++) {
	//	printf("[%d] dset_data_int[%d] = %d\n", rank, i, dset_data_int[i]);
	//}
	H5Dclose(int_dataset_id);
	free(dset_data_int);

	/* Read simple hyperslab selection from var1 */
	int_dataset_id = H5Dopen(file_id, "var1", H5P_DEFAULT);
	dset_data_int = (int *) malloc( H5Dget_storage_size(int_dataset_id) );
	hid_t dataspace = H5Dget_space (int_dataset_id);    /* dataspace handle */
	int dimrank      = H5Sget_simple_extent_ndims (dataspace);
	hsize_t *dims_out = (hsize_t *) malloc( dimrank * sizeof(hsize_t) );
	int status_n  = H5Sget_simple_extent_dims (dataspace, dims_out, NULL);
	//printf("\nDim-Rank: %d -- Dimensions: %lu x %lu x %lu\n", dimrank, (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]), (unsigned long)(dims_out[2]));
	int psizes[NDIMS];
	MPI_Offset starts[NDIMS];
	for (i=0; i<NDIMS; i++)
		psizes[i] = 0;
	MPI_Dims_create(nprocs, NDIMS, psizes);
	starts[0] = rank % psizes[0];
	starts[1] = (rank / psizes[1]) % psizes[1];
	starts[2] = (rank / (psizes[0] * psizes[1])) % psizes[2];
	//printf("[%d] starts = %llu %llu %llu\n", rank, starts[0], starts[1], starts[2]);
	//printf("[%d] psizes = %d %d %d\n", rank, psizes[0], psizes[1], psizes[2]);

	/* Define hyperslab in the dataset. */
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
	//printf("[%d] offset = %llu %llu %llu\n", rank, offset[0], offset[1], offset[2]);
	//printf("[%d] stride = %llu %llu %llu\n", rank, stride[0], stride[1], stride[2]);
	//printf("[%d] count = %llu %llu %llu\n", rank, count[0], count[1], count[2]);
	//printf("[%d] block = %llu %llu %llu\n", rank, block[0], block[1], block[2]);

	/* Define the memory dataspace. */
	hsize_t dimsm[NDIMS];
	for (i=0; i<NDIMS; i++) {
		dimsm[i] = block[i];
	}
	hid_t memspace = H5Screate_simple (NDIMS, dimsm, NULL);

	/* Define memory hyperslab. */
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
	//printf("[%d] data_out size = %d integers.\n", rank, hypersize);
	read_time_hyper = MPI_Wtime();
	status_h5 = H5Dread (int_dataset_id, H5T_NATIVE_INT, memspace, dataspace, dxpl_plist_id, data_out);
	read_time_hyper = MPI_Wtime() - read_time_hyper;
	read_size_hyper = hypersize * sizeof(int);
	MPI_Reduce(&read_size_hyper, &sum_read_size_hyper, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&read_time_hyper, &max_read_time_hyper, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	for (i=0; i<hypersize; i++) {
		if (data_out[i] != ((rank+1)*10)) {
			printf("ERROR!!! ~~~ [%d] data_out[%d] = %d\n", rank, i, data_out[i]);
		}
	}
	free(dset_data_int);
	free(data_out);
	free(dims_out);

	/* Read Attribute */
	hid_t attr1, attr2;
	char attr_data[13];
	char var_attr_data[13];
	attr1 = H5Aopen_name (file_id, "string");
	H5Aread(attr1, H5T_NATIVE_CHAR, attr_data);
	attr2 = H5Aopen_name (int_dataset_id, "varstring");
	H5Aread(attr2, H5T_NATIVE_CHAR, var_attr_data);
	H5Dclose(int_dataset_id);

	/* Summarize the results */
	if (rank == 0) {
		printf(" File Attribute = %s", attr_data);
		printf(" Variable Attribute = %s", var_attr_data);
		printf(" FILE_OPEN:   time[s]=%f\n", max_open_time);
		size_rpt = ((double)sum_read_size_all) / 1048576.0;
		printf(" READ_ALL:    size[MB]=%f  time[s]=%f  bandwidth[MB/s]=%f\n", size_rpt, max_read_time_all, size_rpt/max_read_time_all);
		size_rpt = ((double)sum_read_size_hyper) / 1048576.0;
		printf(" READ_HYPER:  size[MB]=%f  time[s]=%f  bandwidth[MB/s]=%f\n\n", size_rpt, max_read_time_hyper, size_rpt/max_read_time_hyper);
	}

	/* Close File etc */
	H5Aclose(attr1);
	H5Aclose(attr2);
	H5Fclose(file_id);
	H5Pclose(acc_tpl);
	H5Pclose(dxpl_plist_id);
	H5Pclose(fapl);
	H5VLterminate(vol_id);
	H5VLunregister_connector(vol_id);
	assert(H5VLis_connector_registered(connector_name) == 0);

	MPI_Finalize();
	return 0;
}
