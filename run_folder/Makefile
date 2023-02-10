GCC = g++

FFTW_LIB = -lfftw3
LIB  = -lgomp -lm
GCC_LIB = -lblas -llapacke
COPT = -DFFTW -fopenmp -O3
GCC_COPT = -pthread -std=c++11
VPATH = ../:../coagulation:

SOURCES = blas_double.cpp blas_int.cpp cross.cpp matrix.cpp maxvol.cpp parallel_cross_omp.cpp tensor.cpp tensor_train.cpp convolution.cpp wrappers.cpp main.cpp

OBJECTS = $(SOURCES:%.cpp=$(OBJDIR)/%.o)

OBJDIR = build

all: solver.exe

solver.exe: $(OBJECTS)
	$(GCC) $^ -o $@ $(LIB) $(GCC_LIB) $(FFTW_LIB) $(GCC_COPT) $(COPT)

$(OBJDIR)/%.o: %.cpp
	$(GCC) -c $< -o $@ $(GCC_COPT) $(COPT) 

clean:
	rm -rf $(OBJDIR)/*.o *.exe
