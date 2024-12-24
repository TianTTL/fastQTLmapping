CC = gcc

#PLEASE specify the GSL and MKL path here 
GSLROOT = /usr/local/
MKLROOT = /opt/intel/oneapi/mkl/latest

HOSTFILES := fastQTLmapping.cpp
OBJS := fastQTLmapping.o
BIN := fastQTLmapping

override CXXFLAGS += -std=c++11 -O3 -fopenmp
override INCDIR += -DMKL_ILP64 -m64 -I${MKLROOT}/include -I${GSLROOT}/include
override LDFLAGS += -m64 -Wl,--start-group ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lstdc++ -lpthread -lm -ldl ${GSLROOT}/lib/libgsl.a ${GSLROOT}/lib/libgslcblas.a

all: fastQTLmapping.h clipp.h
	@echo "starting compiling fastQTLmapping"
	@$(CC) $(CXXFLAGS) $(INCDIR) -c -o $(OBJS) $(HOSTFILES) 
	@echo "starting linking fastQTLmapping"
	@$(CC) $(CXXFLAGS) -o $(BIN) $(OBJS) $(LDFLAGS)
	@echo "compile done"
	@rm -f $(OBJS)
