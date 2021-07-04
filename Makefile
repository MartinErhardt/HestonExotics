CXX = g++ 
LDD = g++ 
CXXFLAGS= -I src/inc -Wpedantic -Wall -Wextra -oFast -fopenmp -fprofile-generate -flto -march=native -ffloat-store -ffast-math -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -frounding-math -fsingle-precision-constant -fcx-fortran-rules 
CLANGFLAGS = -Wpedantic -Wall -Werror -oFast -fopenmp -march=native -ffast-math -fno-rounding-math -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math -I src/inc
LDFLAGS = -lgcov --coverage  -lcurl
CLANGLDFLAGS = -lstdc++  -lcurl
LIBS=simdjson.cpp
LIB_OBJS=simdjson.o
GETLIBS=curl -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/simdjson.h -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/simdjson.cpp
SRCS = $(shell find -name '*.cpp')
OBJS = $(addsuffix .o,$(basename $(SRCS)))
OBJS_ICC = $(addsuffix _icc.o,$(basename $(SRCS)))
ALL_OBJS=$(shell find -name '*.o')
LIB_OBJS_ICC=simdjson_icc.o
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
