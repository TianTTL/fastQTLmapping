CC = g++

#PLEASE specify the GSL and MKL path here 
GSLROOT = /usr/local/
MKLROOT = /opt/intel/oneapi/mkl/latest

HOSTFILES := fastQTLmapping.cpp
OBJS := fastQTLmapping.o
BIN := fastQTLmapping

override CXXFLAGS += -std=c++11 -O3 -fopenmp 
override INCDIR += -DMKL_ILP64 -m64 -I${MKLROOT}/include -I${GSLROOT}/include
override LDFLAGS += -L${MKLROOT}/lib/intel64 -L${GSLROOT}/lib
override LIBS += -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lgsl -lgslcblas

all: fastQTLmapping.h clipp.h
	@echo "starting compiling fastQTLmapping"
	@$(CC) $(CXXFLAGS) $(INCDIR) -c -o $(OBJS) $(HOSTFILES) 
	@$(CC) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $(BIN) $(OBJS)
	@echo "compile done"
	@rm -f $(OBJS)
