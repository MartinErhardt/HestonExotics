CXX = g++ 
LDD = g++
CXXFLAGS= -I src/inc -Wpedantic -Wall -Wextra -std=c++17 -oFast -fopenmp -fprofile-generate -flto -march=native -ffloat-store -ffast-math -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -fcx-fortran-rules #-fsingle-precision-constant
CLANGFLAGS = -Wpedantic -Wall -Werror -Wextra -std=c++17 -oFast -fopenmp -march=native -ffast-math -fno-rounding-math -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -I src/inc -DEIGEN_USE_MKL_ALL
LDFLAGS = -lgcov --coverage  -lcurl -llevmar -lfftw3 -fopenmp
CLANGLDFLAGS = -lstdc++  -lcurl -llevmar -lfftw3 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
LIBS=simdjson.cpp
LIB_OBJS=simdjson.o
GETLIBS=curl -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/simdjson.h -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/simdjson.cpp
SRCS = $(shell find -name '*.cpp')
OBJS = $(addsuffix .o,$(basename $(SRCS)))
OBJS_ICC = $(addsuffix _icc.o,$(basename $(SRCS)))
ALL_OBJS=$(shell find -name '*.o')
LIB_OBJS_ICC=simdjson_icc.o
#patch < ~/Projekt/HestonExotics/levmar_patch.diff

HestonExotics: $(OBJS) $(LIB_OBJS) 
	$(LDD) $(LDFLAGS) -o HestonExotics $^
	rm $(shell find -name '*.gcda')
simdjson.o:
	$(GETLIBS)
	$(CXX) $(CXXFLAGS) -c -o $(LIB_OBJS) $(LIBS) 
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

%_icc.o: %.cpp
	dpcpp $(CLANGFLAGS) -c -o $@ $^
simdjson_icc.o:
	$(GETLIBS)
	dpcpp $(CLANGFLAGS) -c -o $(LIB_OBJS_ICC) $(LIBS) 
icc: $(OBJS_ICC) $(LIB_OBJS_ICC) 
	dpcpp $(CLANGLDFLAGS) -o HestonExotics_icc $^
clean:
	rm $(ALL_OBJS)
doc: doxygen_conf
	doxygen doxygen_conf
doxygen_conf:
	doxygen -g doxygen_conf
	patch -p0 < doxygen_conf.patch
.PHONY: clean doc
