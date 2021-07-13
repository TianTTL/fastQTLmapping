workDir=/thinker/storage/org/liufanGroup/gaoxj/work/xQTLMapping/2_xQTLMapping
GSLROOT=/thinker/storage/org/liufanGroup/public_lib/gsl-2.6

export C_INCLUDE_PATH=$C_INCLUDE_PATH:${GSLROOT}/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${GSLROOT}/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${GSLROOT}/lib
export LIBRARY_PATH=$LIBRARY_PATH:${GSLROOT}/lib

g++ -std=c++11 -O2 -fopenmp \
-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl \
-lgsl \
-DMKL_ILP64 -m64 -I"${MKLROOT}/include" \
-o $workDir/fastQTLmapping $workDir/fastQTLmapping.cpp
