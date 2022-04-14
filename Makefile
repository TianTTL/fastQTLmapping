CC = g++
GSLROOT = /thinker/storage/org/liufanGroup/public_lib/gsl-2.6
MKLROOT = /thinker/storage/org/liufanGroup/public_software/intel/oneapi/mkl/2021.3.0

HOSTFILES := fastQTLmapping.cpp
OBJS := fastQTLmapping.o
BIN := fastQTLmapping

override CXXFLAGS += -std=c++11 -O3 -fopenmp 
override INCDIR += -DMKL_ILP64 -m64 -I${MKLROOT}/include
override LDFLAGS += -L${MKLROOT}/lib/intel64 -L${GSLROOT}/lib
override LIBS += -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lgsl -lgslcblas

all: fastQTLmapping.h clipp.h
	@echo "starting compiling fastQTLmapping"
	@$(CC) $(CXXFLAGS) $(INCDIR) -c -o $(OBJS) $(HOSTFILES) 
	@$(CC) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $(BIN) $(OBJS)
	@echo "compile done"
	@rm -f $(OBJS)
