LIBS=simdjson.cpp
GETLIBS=curl -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/simdjson.h -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/simdjson.cpp
SRCS = $(shell find -name '*.cpp')
OBJS = $(addprefix ./bin/,$(addsuffix .o,$(basename $(SRCS))))
OBJS_ICC = $(addprefix ./bin/,$(addsuffix _icc.o,$(basename $(SRCS))))
ALL_OBJS=$(addprefix ./bin/,$(shell find -name '*.o'))
#patch < ~/Projekt/HestonExotics/levmar_patch.diff
ifeq ($(CXX), dpcpp)
	CXXFLAGS = -Wpedantic -Wall -Werror -Wextra -std=c++17 -oFast -fopenmp -march=native -ffast-math -fno-rounding-math -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -I src/inc -DEIGEN_USE_MKL_ALL
	LDFLAGS = -lstdc++  -lcurl -llevmar -lfftw3 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
else
	CXX = g++
	CXXFLAGS= -I src/inc -Wpedantic -Wall -Werror -Wextra -std=c++17 -fopenmp -oFast -fprofile-generate -flto -march=native -ffloat-store -ffast-math -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -fcx-fortran-rules #-fsingle-precision-constant
	CXXDBGFLAGS= -I src/inc -Wpedantic -Wall -Werror -Wextra -std=c++17 -g -fopenmp
	LDFLAGS = -lcurl -llevmar -lfftw3 -fopenmp -lgcov --coverage 
endif
HestonExotics: $(OBJS) $(LIBS)
	$(CXX) $(LDFLAGS) -o HestonExotics $(OBJS)
	rm -f *.gcda 2> /dev/null
simdjson.cpp:
	$(GETLIBS)
bin/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^
clean:
	rm bin/src/*
fclean:
	rm -f bin/simdjson.o
	rm bin/src/*
doc: doxygen_conf
	doxygen doxygen_conf
doxygen_conf:
	doxygen -g doxygen_conf
	patch -p0 < doxygen_conf.patch
.PHONY: clean doc
