#include "adios_vol.h"
#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
// typedef int hsize_t;

int cmpfuncWobble(const void *a, const void *b) {
  H5_posMap *wa = (H5_posMap *)a;
  H5_posMap *wb = (H5_posMap *)b;

  return (wa->posInSource - wb->posInSource);
}

void traverseBlock(int ndim, int currDim, hsize_t *start, hsize_t *count,
                   hsize_t *m, hsize_t previous, H5_posMap **result) {
  hsize_t k, n;
  if (ndim == currDim + 1) {
    for (k = 0; k < count[currDim]; k++) {
      hsize_t pos = start[currDim] + k;
      // printf("... pos = %llu\n", previous + pos * m[currDim]);

      (**result).posInSource = previous + pos * m[currDim];
      hsize_t nextCompactPos = (**result).posCompact + 1;
      (*result)++;
      (*result)->posCompact = nextCompactPos;
    }
    return;
  }

  for (k = start[currDim]; k < start[currDim] + count[currDim]; k++) {
    hsize_t base = k * m[currDim] + previous;
    traverseBlock(ndim, currDim + 1, start, count, m, base, result);
  }
}

void getMultiplier(int ndim, hsize_t *dims, hsize_t *result) {
  int n, k;
  for (n = 0; n < ndim; n++) {
    result[n] = 1;
    for (k = n + 1; k < ndim; k++) {
      result[n] *= dims[k];
    }
  }
}

void GetSelOrder(hid_t space_id, H5S_sel_type space_type, H5_posMap **result) {
  if (space_type == H5S_SEL_ALL) {
    return;
  }

  int ndim = H5Sget_simple_extent_ndims(space_id);
  hsize_t dims[ndim];
  H5Sget_simple_extent_dims(space_id, dims, NULL);

  hsize_t mm[ndim];
  getMultiplier(ndim, dims, mm);

  if (space_type == H5S_SEL_HYPERSLABS) {
    hsize_t nblocks = H5Sget_select_hyper_nblocks(space_id);
    hsize_t total = H5Sget_select_npoints(space_id);
    *result = (H5_posMap *)malloc(sizeof(H5_posMap) * (total + 1));
    (*result)->posCompact = 0;

    // NOTE: a simple bounding box is not treated as one block
    // rather each element is a block, thus nblocks = total
    // for a simple bounding box.  not sure how to work around it
    // if (nblocks == total) {
    // this call crashes
    //  H5S_get_select_bounds(space_id, boundstart, boundends);
    //}
    hsize_t *blockinfo =
        (hsize_t *)SAFE_MALLOC(sizeof(hsize_t) * 2 * ndim * nblocks);
    herr_t status =
        H5Sget_select_hyper_blocklist(space_id, (hsize_t)0, nblocks, blockinfo);

    H5_posMap *startPointer = *result;

    int i = 0, n = 0;

    for (n = 0; n < nblocks; n++) {
      uint64_t blockSize = 1;
      uint64_t h_start[ndim], h_count[ndim];
      for (i = 0; i < ndim; i++) {
        int pos = 2 * ndim * n;
        h_start[i] = blockinfo[pos + i];
        h_count[i] = blockinfo[pos + ndim + i] - h_start[i] + 1;
        blockSize *= h_count[i];
      }

      // printf(" got block: start[%llu, %llu], count[%llu, %llu]\n",
      // h_start[0], h_start[1], h_count[0], h_count[1]);
      traverseBlock(ndim, 0, h_start, h_count, mm, 0, result);
    }

    // now reorder the pos
    qsort(startPointer, total, sizeof(H5_posMap), cmpfuncWobble);

    *result = startPointer;
    SAFE_FREE(blockinfo);
  }

  // H5S_SEL_POINTS
}

void assignToMemSpace(H5_posMap *sourceSelOrderInC,
                      H5_posMap *targetSelOrderInC, hsize_t total,
                      size_t dataTypeSize, char *adiosData, char *buf) {
  size_t k;
  hsize_t n;

  if (sourceSelOrderInC == NULL) { // fileSpace=ALL memSpace=partial
    for (n = 0; n < total; n++) {
      hsize_t dd = n;
      hsize_t ss = targetSelOrderInC[n].posInSource;
      for (k = 0; k < dataTypeSize; k++) {
        buf[ss * dataTypeSize + k] = adiosData[dd * dataTypeSize + k];
      }
    }
    return;
  }
  if (targetSelOrderInC == NULL) {
    for (n = 0; n < total; n++) {
      // hsize_t ss = targetSelOrderInC[n].posInSource;
      hsize_t dd = sourceSelOrderInC[n].posCompact;
      hsize_t ss = n;
      for (k = 0; k < dataTypeSize; k++) {
        buf[ss * dataTypeSize + k] = adiosData[dd * dataTypeSize + k];
      }
    }
    return;
  }

  for (n = 0; n < total; n++) {
    hsize_t ss = targetSelOrderInC[n].posInSource;
    hsize_t dd = sourceSelOrderInC[n].posCompact;

    // printf("dd[%llu, %llu]\n", sourceSelOrderInC[n].posCompact,
    // sourceSelOrderInC[n].posInSource);  printf("n=%llu, total=%llu, ss=%llu
    // ts=%lu\n", n, total, ss*dataTypeSize, dataTypeSize);
    for (k = 0; k < dataTypeSize; k++) {
      buf[ss * dataTypeSize + k] = adiosData[dd * dataTypeSize + k];
    }
  }
}
