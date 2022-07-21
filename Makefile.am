bin_PROGRAMS = hexo
hexo_SOURCES = src/inc/AsianContract.h\
  src/inc/DB.h\
  src/inc/HDistribution.h\
  src/inc/Main.h\
  src/inc/RNG.h\
  src/inc/SWIFT.h\
  src/inc/UI.h\
  src/inc/VanillaContract.h\
  src/inc/BSM.h\
  src/inc/HCalibration.h\
  src/inc/HSimulation.h\
  src/inc/my_math.h\
  src/inc/SDE.h\
  src/inc/Types.h\
  src/inc/UnitTest.h\
  src/inc/WebAPI.h\
  src/HSimulation.tpp\
  src/BSM.cpp\
  src/DB.cpp\
  src/HCalibration.cpp\
  src/HDistribution.cpp\
  src/Main.cpp\
  src/RNG.cpp\
  src/SWIFT.cpp\
  src/UI.cpp\
  src/UnitTest.cpp\
  src/as241.f90\
  src/WebAPI.cpp\
  simdjson.cpp
SJSON_SRC=simdjson.cpp simdjson.h
SJSON_OBJ=bin/simdjson.o
SHISHUA_INC=shishua/shishua.h
if DEBUG
  AM_CPPFLAGS= -I src/inc -I shishua -I levmar-2.6 -DFP_SIZE=8
  AM_CXXFLAGS= -Wpedantic -Wall -Wextra -Werror -g -fsanitize=undefined -fsanitize=address
  AM_LDFLAGS = -lgfortran -lcurl -lsqlite3 -lfftw3 -fopenmp -flto -lblas -llapack -g -fsanitize=undefined -fsanitize=address
  LDADD = levmar-2.6/liblevmar.a
else
  AM_CPPFLAGS= -I src/inc -I shishua -I levmar-2.6 -DFASTMATH -DNDEBUG -DFP_SIZE=8
  AM_CXXFLAGS= -ffast-math\
    -fno-rounding-math\
    -fno-signaling-nans\
    -fcx-limited-range\
    -fno-math-errno\
    -funsafe-math-optimizations\
    -fassociative-math\
    -freciprocal-math\
    -ffinite-math-only\
    -fno-signed-zeros\
    -fno-trapping-math
  AM_LDFLAGS = -lgfortran -lcurl -lsqlite3 -lfftw3 -fopenmp -flto -lblas -llapack #-lgcov --coverage
  LDADD = levmar-2.6/liblevmar.a
endif
EXTRA_DIST = *.diff 

prereq: | bin_dirs $(SJSON_SRC) $(SHISHUA_INC) levmar-2.6 
bin_dirs:
	@echo ${hexo_SOURCES}
	mkdir -p bin
	mkdir -p bin/src
$(SJSON_SRC):%:
	curl -O https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/$@
$(SHISHUA_INC):%:
	git clone https://github.com/espadrine/shishua.git
	patch -p0 < shishua_inline_patch.diff
$(AS241_OBJ):%:$(AS241_SRC)
	$(FC) $(FFLAGS) -c $^ -o $@
levmar-2.6:
	curl -O http://users.ics.forth.gr/~lourakis/levmar/levmar-2.6.tgz
	tar -xzf levmar-2.6.tgz
	patch -p0 < levmar_patch.diff
	cd levmar-2.6 && make
	rm levmar-2.6.tgz
loc:
	find . -regextype posix-extended -regex "./src/.*(.h|.cpp|.tpp|.f90)" | xargs wc -l
doc: doxygen_conf
	doxygen doxygen_conf
doxygen_conf:
	doxygen -g doxygen_conf
	patch -p0 < doxygen_conf.diff
#clean-hexo:
#	@rm hexo_OBJECTS
#clean-local:
#	@rm -rf bin/*
#	@rm config.status configure config.log
#	@rm Makefile
#	@rm -r autom4te.cache/
#	@rm aclocal.m4
#	@rm install-sh missing Makefile.in
.PHONY: clean-local doc bin_dirs loc
