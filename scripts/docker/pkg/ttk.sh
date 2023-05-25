#! /bin/bash
set -e

require-pkgs \
    build-essential         \
    cmake                   \
    curl                    \
    libboost-system-dev     \
    libcgns-dev             \
    libeigen3-dev           \
    libexpat1-dev           \
    libfreetype6-dev        \
    libhdf5-dev             \
    libjpeg-dev             \
    libjsoncpp-dev          \
    liblz4-dev              \
    liblzma-dev             \
    libnetcdf-cxx-legacy-dev\
    libnetcdf-dev           \
    libogg-dev              \
    libpng-dev              \
    libprotobuf-dev         \
    libpugixml-dev          \
    libsqlite3-dev          \
    libgraphviz-dev	        \
    libtheora-dev           \
    libtiff-dev             \
    libxml2-dev             \
    ninja-build             \
    protobuf-compiler       \
    python3-dev             \
    python3-numpy-dev       \
    wget                    \
    zlib1g-dev              \
    git



/sbin/ldconfig

if [ -n "${DEV}" ]; then
        #echo "DEVELOPER MODE"
        exit
fi

# get source code
git clone https://github.com/m-s-will/ttk.git
cd ttk
git checkout mpi_container

# actually compile
cmake-default \
    -DTTK_BUILD_DOCUMENTATION=OFF \
    -DTTK_BUILD_PARAVIEW_PLUGINS=ON \
    -DTTK_BUILD_STANDALONE_APPS=OFF \
    -DTTK_BUILD_VTK_WRAPPERS=ON \
    -DTTK_BUILD_VTK_PYTHON_MODULE=OFF \
    -DTTK_ENABLE_DOUBLE_TEMPLATING=OFF \
    -DTTK_ENABLE_CPU_OPTIMIZATION=OFF \
    -DTTK_ENABLE_OPENMP=ON \
    -DTTK_ENABLE_KAMIKAZE=ON \
    -DTTK_ENABLE_MPI=ON \
    -DVTK_MODULE_ENABLE_ttkExTreeM=NO \
    -DVTK_MODULE_ENABLE_ttkPairExtrema=NO \
    -DVTK_MODULE_ENABLE_ttkGradientGraph=NO \
    ..

# call Ninja manually to ignore duplicate targets
# cmake --build .

# ninja -w dupbuild=warn install
# cmake --install .

# popd
