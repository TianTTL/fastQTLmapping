# workDir=fastQTLmapping_install_path
# GSLROOT=GSL_install_directory
# MKLROOT=MKL_install_directory

g++ -std=c++11 -fopenmp \
 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl \
 -L${GSLROOT}/lib -lgsl -lgslcblas \
 -DMKL_ILP64  -m64  -I"${MKLROOT}/include" \
 -o $workDir/fastQTLmapping \
 -c $workDir/fastQTLmapping.cpp

