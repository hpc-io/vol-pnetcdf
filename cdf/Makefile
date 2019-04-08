FNAME=vol_test

VOLNAME=cdf_vol
EXE=${FNAME}.ex

# Theta
HDF5_ROOT=/home/zamora/hdf5_root_dir
BLD_DIR=theta-develop
MPICC=cc -std=c99

default: ${EXE}

PNETCDF_INSTALL_DIR=${PARALLEL_NETCDF_DIR}
HDF5_INSTALL_DIR=${HDF5_ROOT}/library/install/${BLD_DIR}
HDF5_SRC_DIR=${HDF5_ROOT}/xgitlabrepos/hdf5/src

${VOLNAME}.o: ${VOLNAME}.c ${VOLNAME}.h
	${MPICC} -c -g -O3 -I${PNETCDF_INSTALL_DIR}/include -I${HDF5_INSTALL_DIR}/include ${VOLNAME}.c -o ${VOLNAME}.o
	#${MPICC} -c -g -O3 -I${PNETCDF_INSTALL_DIR}/include -I${HDF5_INSTALL_DIR}/include -I${HDF5_SRC_DIR} ${VOLNAME}.c -o ${VOLNAME}.o

${FNAME}.o: ${FNAME}.c ${VOLNAME}.o
	${MPICC} -c -g -O3 -I${PNETCDF_INSTALL_DIR}/include -I${HDF5_INSTALL_DIR}/include ${FNAME}.c -o ${FNAME}.o
	#${MPICC} -c -g -O3 -I${PNETCDF_INSTALL_DIR}/include -I${HDF5_INSTALL_DIR}/include -I${HDF5_SRC_DIR} ${FNAME}.c -o ${FNAME}.o

${EXE}: ${FNAME}.o ${VOLNAME}.o
	${MPICC} ${FNAME}.o ${VOLNAME}.o -o ${EXE} -L${HDF5_INSTALL_DIR}/lib -lhdf5 -L${PNETCDF_INSTALL_DIR}/lib -lpnetcdf

install:
	mv -f ${EXE} bin/${EXE}

clean:
	rm -f ${FNAME}.o
	rm -f ${VOLNAME}.o
	rm -f bin/${EXE}
