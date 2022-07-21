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
if NO_CURL
  LCURL=curl/lib/libcurl.a
  CURL_INC=-I curl/include
endif
#if NO_SSL
#  LSSL=ssl/lib/libssl.a
#  LCRYPTO=ssl/lib/libcrypto.a
#endif
if NO_SQLITE
  LSQLITE = sqlite/lib/libsqlite3.a
  SQLITE_INC=-I ./sqlite/include
endif
if NO_FFTW
  LFFTW = fftw/lib/libfftw3.a
  FFTW_INC=-I ./fftw/include
endif
DEP_INC=${SQLITE_INC} ${FFTW_INC} ${CURL_INC}
DEP_LOCAL=${LFFTW} ${LSQLITE} ${LCURL} ${LSSL}
if DEBUG
  AM_CPPFLAGS= -I src/inc -I shishua -I levmar-2.6 ${DEP_INC} -DFP_SIZE=8
  AM_CXXFLAGS= -Wpedantic -Wall -Wextra -Werror -g -fsanitize=undefined -fsanitize=address
  AM_LDFLAGS = -lgfortran  -fopenmp -flto -lblas -llapack -g -fsanitize=undefined -fsanitize=address
  LDADD = levmar-2.6/liblevmar.a ${DEP_LOCAL}
else
  AM_CPPFLAGS= -I src/inc -I shishua -I levmar-2.6 ${DEP_INC} -DFASTMATH -DNDEBUG -DFP_SIZE=8
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
  AM_LDFLAGS = -lgfortran -fopenmp -flto -lblas -llapack -lssl -lcrypto#-lgcov --coverage
  LDADD = levmar-2.6/liblevmar.a ${DEP_LOCAL}
endif

EXTRA_DIST = *.diff 

prereq: | bin_dirs $(SJSON_SRC) $(SHISHUA_INC) levmar-2.6 ${DEP_LOCAL}
bin_dirs:
	@echo ${hexo_SOURCES}
	mkdir -p bin
	mkdir -p bin/src
$(SJSON_SRC):%:
	wget -c https://raw.githubusercontent.com/simdjson/simdjson/master/singleheader/$@
$(SHISHUA_INC):%:
	git clone https://github.com/espadrine/shishua.git
	patch -p0 < shishua_inline_patch.diff
$(AS241_OBJ):%:$(AS241_SRC)
	$(FC) $(FFLAGS) -c $^ -o $@
levmar-2.6:
	wget -c http://users.ics.forth.gr/~lourakis/levmar/levmar-2.6.tgz
	tar -xzf levmar-2.6.tgz
	patch -p0 < levmar_patch.diff
	cd levmar-2.6 && make
	rm levmar-2.6.tgz
SQLITE_VER=3390100
${LSQLITE}:
	wget -c https://www.sqlite.org/2022/sqlite-autoconf-${SQLITE_VER}.tar.gz -4
	tar -xzf sqlite-autoconf-${SQLITE_VER}.tar.gz
	mkdir -p sqlite
	cd sqlite-autoconf-${SQLITE_VER} &&\
	  ./configure --prefix=${PWD}/sqlite &&\
	  make &&\
	  make install &&\
	  cd ..
FFTW_VER=3.3.10
${LFFTW}:
	wget -c ftp://ftp.fftw.org/pub/fftw/fftw-${FFTW_VER}.tar.gz
	tar -xzf fftw-${FFTW_VER}.tar.gz
	mkdir -p fftw
	cd fftw-${FFTW_VER} &&\
	  ./configure --prefix=${PWD}/fftw &&\
	  make &&\
	  make install &&\
	  cd ..
#SSL_VER=1.1.1
#${LSSL} ${LCR1YPTO}:
#	wget -c https://www.openssl.org/source/openssl-${SSL_VER}q.tar.gz -4
#	tar -xzf openssl-${SSL_VER}q.tar.gz
#	mkdir -p ssl
#	cd openssl-${SSL_VER}q &&\
#	  ./config --prefix=${PWD}/ssl --no-comp --openssldir=${PWD}/ssl &&\
#	  make depend &&\
#	  make &&\
#	  make install &&\
#	  cd ..
CURL_VER=7.84.0
${LCURL}:	
	wget -c https://curl.se/download/curl-${CURL_VER}.tar.gz -4
	tar -xzf curl-${CURL_VER}.tar.gz
	mkdir -p curl
	cd curl-${CURL_VER} &&\
	  ./configure --with-ssl=${PWD}/ssl\
	    --without-zlib\
	    --without-nghttp2\
	    --without-libidn2\
	    --without-brotli\
	    --without-zstd\
	    --disable-thread\
	    --disable-shared\
	    --prefix=${PWD}/curl &&\
	  make &&\
	  make install &&\
	  cd ..
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
