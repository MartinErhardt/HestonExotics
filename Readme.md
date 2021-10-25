
![GCC](https://img.shields.io/static/v1?logo=github&label=GCC11&message=passing&color=Blue)
![ICC](https://img.shields.io/static/v1?logo=github&label=ICC&message=passing&color=Blue)
![MinGW64](https://img.shields.io/static/v1?logo=github&label=MinGW64&message=passing&color=Blue)
![License](https://img.shields.io/static/v1?label=License&message=MPLv2&color=blue)
![docs](https://img.shields.io/static/v1?label=docs&message=doxygen&color=green)

hexo
===============================================
Tool for pricing exotic options in the Heston model
### Installation
The project can be compiled and installed using (GNU-)make.
```
$ make [-j4] [CXX=dpcpp|g++]
/** compile and install project */
$ make doc
/** compile documentation and install into doc/html/index.html*/
```
### Known issues
- RNG generation is 3x faster when using Intel compiler(further profiling/disassembly required)
- standard out is locked by curl when downloading
- if compiled in MinGW64 standard
### Dependencies
- blas
- curl
- eigen3
- f2c
- fftw3
- lapack
- libcurl
- sqlite3
