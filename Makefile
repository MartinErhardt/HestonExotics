SRCS = $(shell find -name '*.cpp')
OBJS = $(addprefix ./bin/,$(addsuffix .o,$(basename $(SRCS))))
SJSON_SRC=simdjson.cpp simdjson.h
SJSON_OBJ=bin/simdjson.o
SHISHUA_INC=shishua/shishua.h
#SHISHUA_OBJ=bin/shishua.o
AS241_SRC=src/as241.f90
AS241_OBJ=bin/as241.o
FC=gfortran
INCS= -I src/inc -I shishua -I levmar-2.6
MACROS= -DFP_SIZE=8

#patch < ~/Projekt/HestonExotics/levmar_patch.diff
ifeq ($(CXX), dpcpp)
	FFLAGS= -std=f2008ts -fdefault-real-8
	CXXFLAGS = $(INCS) -Wpedantic -Wall -Wextra -std=c++17 -oFast -fopenmp -march=native -DFASTMATH -ffast-math -fno-rounding-math -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -DEIGEN_USE_MKL_ALL $(MACROS)
	LDFLAGS = -lstdc++  -lcurl -lsqlite3 -llevmar -I levmar-2.6 -lfftw3 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -L./levmar-2.6 -llevmar
else
	CXX = g++
	ifeq ($(DBG), true)
	    FFLAGS= -std=f2008ts -fdefault-real-8
	    CXXFLAGS= $(INCS) -Wpedantic -Wall -Wextra -std=c++17 -g -fopenmp  -march=native $(MACROS)
	else
	    FFLAGS=  -std=f2008ts -fdefault-real-8 -flto
	    CXXFLAGS= $(INCS) -Wpedantic -Wall -Wextra -std=c++17  -fopenmp -oFast  -march=native -ffloat-store -DFASTMATH -ffast-math -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -fcx-fortran-rules $(MACROS) -flto #-fsingle-precision-constant -fprofile-generate;
	endif
	LDFLAGS = -lgfortran -lcurl -lsqlite3 -lfftw3 -fopenmp -flto -lblas -llapack -L./levmar-2.6 -llevmar  #-lgcov --coverage 
endif

all: | bin_dirs hexo
bin_dirs:
	mkdir -p bin
	mkdir -p bin/src
hexo: $(SJSON_OBJ) $(SHISHUA_INC) $(AS241_OBJ) levmar-2.6 $(OBJS) 
	$(CXX) -o hexo  $(OBJS)  $(AS241_OBJ) $(LDFLAGS)
	rm -f *.gcda 2> /dev/null
$(SJSON_SRC):%:
	curl -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/$@
$(SJSON_OBJ):%:$(SJSON_SRC)
	$(CXX) $(CXXFLAGS) -c -o $@ $(basename $(notdir $@)).cpp
$(SHISHUA_INC):%:
	git clone https://github.com/espadrine/shishua.git
	patch -p0 < shishua_inline_patch.diff
#$(AS241_SRC):%:
#	curl -o $@ http://lib.stat.cmu.edu/apstat/241 
$(AS241_OBJ):%:$(AS241_SRC)
	$(FC) $(FFLAGS) -c $^ -o $@
levmar-2.6:
	curl -O http://users.ics.forth.gr/~lourakis/levmar/levmar-2.6.tgz
	tar -xzf levmar-2.6.tgz
	patch -p0 < levmar_patch.diff
	cd levmar-2.6 && make
	rm levmar-2.6.tgz
#$(SHISHUA_OBJ):%:$(SHISHUA_SRC)
#	$(CXX) $(CXXFLAGS) -DHEADER='"shishua.h"' -c -o $@ $(SHISHUA_SRC)
bin/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^
clean:
	rm bin/src/*
fclean:
	rm -f $(SJSON_OBJ) $(AS241_OBJ)
	rm bin/src/*
loc:
	find . -regextype posix-extended -regex "./src/.*(.h|.cpp|.f90)" | xargs wc -l
doc: doxygen_conf
	doxygen doxygen_conf
doxygen_conf:
	doxygen -g doxygen_conf
	patch -p0 < doxygen_conf.diff
.PHONY: clean doc bin_dirs hexo loc
