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

#define H5VL_CFD_MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

extern const H5VL_class_t H5VL_cdf_g;
extern hid_t native_connector_id;

/* Pass-through VOL connector info */
typedef struct H5VL_cdf_info_t {
    hid_t vol_id; /* VOL ID */
} H5VL_cdf_info_t;


#ifdef __cplusplus
extern "C" {
#endif

H5_DLL hid_t H5VL_cdf_register(void);

#ifdef __cplusplus
}
#endif

#endif /* _H5VLcdf_H */
