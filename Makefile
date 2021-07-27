SRCS = $(shell find -name '*.cpp')
OBJS = $(addprefix ./bin/,$(addsuffix .o,$(basename $(SRCS))))
SJSON_SRC=simdjson.cpp simdjson.h
SJSON_OBJ=bin/simdjson.o
SHISHUA_INC=shishua/shishua.h
#SHISHUA_OBJ=bin/shishua.o
AS241_SRC=src/as241.f90
AS241_OBJ=bin/as241.o
FC=gfortran

#patch < ~/Projekt/HestonExotics/levmar_patch.diff
ifeq ($(CXX), dpcpp)
	FFLAGS= -std=f2008ts -fdefault-real-8
	CXXFLAGS = -Wpedantic -Wall -Wextra -std=c++17 -oFast -fopenmp -march=native -DFASTMATH -ffast-math -fno-rounding-math -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -I src/inc -I shishua -DEIGEN_USE_MKL_ALL -DFP_SIZE=8
	LDFLAGS = -lstdc++  -lcurl -llevmar -lfftw3 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
else
	CXX = g++
	ifeq ($(DBG), true)
	    FFLAGS= -std=f2008ts -fdefault-real-8
	    CXXFLAGS= -I src/inc -I shishua -Wpedantic -Wall -Werror -Wextra -std=c++17 -g -fopenmp  -march=native
	else
	    FFLAGS= -std=f2008ts -fdefault-real-8 -flto
	    CXXFLAGS= -I src/inc -I shishua -Wpedantic -Wall -Wextra -std=c++17  -fopenmp -oFast  -march=native -ffloat-store -DFASTMATH -ffast-math -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -fcx-fortran-rules -DFP_SIZE=8 -flto #-fsingle-precision-constant -fprofile-generate;
	endif
	LDFLAGS = -lgfortran -lcurl  -llevmar -lfftw3 -fopenmp -flto #-lgcov --coverage 
endif

all: | bin_dirs HestonExotics
bin_dirs:
	mkdir -p bin
	mkdir -p bin/src
HestonExotics: $(SJSON_OBJ) $(SHISHUA_INC) $(AS241_OBJ) $(OBJS) 
	$(CXX) -o HestonExotics  $(OBJS)  $(AS241_OBJ) $(LDFLAGS)
	rm -f *.gcda 2> /dev/null
$(SJSON_SRC):%:
	curl -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/$@
$(SJSON_OBJ):%:$(SJSON_SRC)
	$(CXX) $(CXXFLAGS) -c -o $@ $(basename $(notdir $@)).cpp
$(SHISHUA_INC):%:
	git clone https://github.com/espadrine/shishua.git
	patch shishua/shishua-avx2.h shishua_inline_patch.diff
#$(AS241_SRC):%:
#	curl -o $@ http://lib.stat.cmu.edu/apstat/241 
$(AS241_OBJ):%:$(AS241_SRC)
	$(FC) $(FFLAGS) -c $^ -o $@
	
#$(SHISHUA_OBJ):%:$(SHISHUA_SRC)
#	$(CXX) $(CXXFLAGS) -DHEADER='"shishua.h"' -c -o $@ $(SHISHUA_SRC)
bin/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^
clean:
	rm bin/src/*
fclean:
	rm -f $(SJSON_OBJ) $(AS241_OBJ)
	rm bin/src/*
doc: doxygen_conf
	doxygen doxygen_conf
doxygen_conf:
	doxygen -g doxygen_conf
	patch -p0 < doxygen_conf.patch
.PHONY: clean doc bin_dirs HestonExotics
