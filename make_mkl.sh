workDir=/data/gaoxj/work/xQTLmapping/2_xQTLMapping
GSLROOT=/data/gaoxj/lib/gsl-2.7
gccLibDir=/lib64

source /data/gaoxj/lib/intel/setvars.sh --force

export C_INCLUDE_PATH=$C_INCLUDE_PATH:${GSLROOT}/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${GSLROOT}/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${GSLROOT}/lib
export LIBRARY_PATH=$LIBRARY_PATH:${GSLROOT}/lib

g++ -std=c++11 -O2\
 -DMKL_ILP64  -m64  -I"${MKLROOT}/include"\
 -o $workDir/fastQTLmapping.o -c $workDir/fastQTLmapping.cpp

g++ -std=c++11 -fopenmp\
 -static\
 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${gccLibDir}/libpthread.a ${gccLibDir}/libm.a ${gccLibDir}/libdl.a ${GSLROOT}/lib/libgsl.a ${GSLROOT}/lib/libgslcblas.a $workDir/fastQTLmapping.o -Wl,--end-group\
 -o $workDir/fastQTLmapping 

