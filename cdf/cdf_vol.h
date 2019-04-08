/*
 * Purpose:	The public header file for the cdf VOL connector.
 */

#ifndef _H5VLcdf_H
#define _H5VLcdf_H

/* Identifier for the cdf VOL connector */
#define H5VL_CDF	(H5VL_cdf_register())

/* Characteristics of the cdf VOL connector */
#define H5VL_CDF_NAME        "cdf"
#define H5VL_CDF_VALUE       611           /* VOL connector ID */
#define H5VL_CDF_VERSION     0
#define H5VL_CDF_CACHE_SIZE  1024
#define MAX_FILES_OPEN       1024
#define LITTLE_ENDIAN_CDFVL  1            /* Define Endianess of Machine */

#define H5VL_CFD_MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef enum {
	CDF_VOL_ATT,
	CDF_VOL_VAR
} cdf_vol_get_t;

extern const H5VL_class_t H5VL_cdf_g;

/* Extra API */
herr_t cdf_vol_var_get_attname(const char *file_name, const char *var_name, char *att_name, int iatts);
herr_t cdf_vol_var_get_natts(const char *file_name, const char *var_name, int *natts);
herr_t cdf_vol_file_get_nitems(const char *file_name, int *nitems, cdf_vol_get_t type);
herr_t cdf_vol_file_get_iname(const char *file_name, char *iname, int item, cdf_vol_get_t type);

#ifdef __cplusplus
extern "C" {
#endif

H5_DLL hid_t H5VL_cdf_register(void);

#ifdef __cplusplus
}
#endif

#endif /* _H5VLcdf_H */
