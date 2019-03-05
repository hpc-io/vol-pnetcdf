FNAME=vol_test

VOLNAME=cdf_vol
EXE=${FNAME}.ex

# Mac
HDF5_ROOT=/Users/rzamora/IO/CCIO
BLD_DIR=ccio
MPICC=/Users/rzamora/IO/CCIO/library/install/ccio/bin/h5pcc
#MPICC=mpicc

default: ${EXE}

HDF5_INSTALL_DIR=${HDF5_ROOT}/library/install/${BLD_DIR}

${VOLNAME}.o: ${VOLNAME}.c ${VOLNAME}.h
	${MPICC} -c -g -O3 -I${HDF5_INSTALL_DIR}/include ${VOLNAME}.c -o ${VOLNAME}.o

${FNAME}.o: ${FNAME}.c ${VOLNAME}.o
	${MPICC} -c -g -O3 -I${HDF5_INSTALL_DIR}/include ${FNAME}.c -o ${FNAME}.o

${EXE}: ${FNAME}.o ${VOLNAME}.o
	${MPICC} ${FNAME}.o ${VOLNAME}.o -o ${EXE} -L${HDF5_INSTALL_DIR}/lib -lhdf5

clean:
	rm -f ${FNAME}.o
	rm -f ${VOLNAME}.o
	rm -f ${EXE}
	rm -f *.h5
