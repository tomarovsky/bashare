module load gcc/gcc-11.1
export PATH=${PATH}:/mnt/lustre/tools/texinfo/texinfo-6.7/bin

# compile time varts
CURL="/home/tools/curl/curl-7.84.0/"
PCRE="/home/tools/pcre/pcre-8.39/"
XZ="/home/tools/xz/xz-5.2.3/"
ZLIB="/home/tools/zlib/zlib-1.2.8/"
BZIP="/home/tools/bzip2/bzip2-1.0.6/"

IFLAGS="-I${CURL}/include -I${PCRE}/include -I${XZ}/include -I${BZIP}/include -I${ZLIB}/include -I/mnt/lustre/tools/gcc/gcc-11.1/include"
LFLAGS="-L${CURL}/lib -L${PCRE}/lib -L${XZ}/lib -L${BZIP}/lib -L${ZLIB}/lib -L/mnt/lustre/tools/gcc/gcc-11.1/lib64"
PKGTIFF=/home/tools/libtiff/libtiff-4.1.0/lib/pkgconfig

export CFLAGS=${IFLAGS}
export CXXFLAGS=${IFLAGS}
export LDFLAGS=${LFLAGS}
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${PKGTIFF}

export LD_LIBRARY_PATH="/home/tools/boost/boost_1_76_0_py38/lib:/home/tools/ncl/ncl-2.1.18/lib/ncl:/home/tools/ncl/ncl-2.1.18/lib"
export LIBRARY_PATH=${LIBRARY_PATH}:/mnt/lustre/tools/xz/xz-5.2.3/lib/
export PATH=/mnt/lustre/tools/curl/curl-7.84.0/bin:/home/tools/cmake/cmake-3.21.0/bin:${PATH}
