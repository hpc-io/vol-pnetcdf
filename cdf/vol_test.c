#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "hdf5.h"
#include "cdf_vol.h"
#define FILE "testfile.h5"

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

	if(argc !=2){
		printf("Input: connector name, e.g., cdf\n");
		return 0;
	}
	strcpy(connector_name,argv[1]);


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

	return 0;
}
