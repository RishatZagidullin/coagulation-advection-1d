GCC = g++


FFTW_INSTALL_PATH = ${FFTW}#ZHORES
FFTWCFLAGS=-I${FFTW_INSTALL_PATH}/include
FFTWLDFLAGS=-L${FFTW_INSTALL_PATH}/lib
FFTW_LIB = -lfftw3

LIB  = -lgomp -lm
GCC_LIB = -lblas -llapacke
COPT = -DFFTW -fopenmp -O3
GCC_COPT = -pthread -std=c++11
NVCC_COPT = -Xcompiler -std=c++11
VPATH = ../:../coagulation:

SOURCES = blas_double.cpp blas_int.cpp cross.cpp matrix.cpp maxvol.cpp parallel_cross_omp.cpp tensor.cpp tensor_train.cpp convolution.cpp wrappers.cpp main.cpp
OBJECTS = $(SOURCES:.cpp=_cpp.o)

all: solver.exe

solver.exe: $(OBJECTS)
	$(GCC) $^ -o $@ $(LIB) $(GCC_LIB) $(FFTWLDFLAGS) $(FFTW_LIB) $(GCC_COPT) $(COPT)

%_cpp.o: %.cpp
	$(GCC) -c $< -o $@ $(COPT) $(FFTWCFLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe
