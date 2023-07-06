ompi=4.1.5
ompi_url=4.1
ompi_prefix=openmpi-$ompi

require-pkgs \
	wget			\
	build-essential	\

wget https://download.open-mpi.org/release/open-mpi/v$ompi_url/$ompi_prefix.tar.gz
tar xvzf $ompi_prefix.tar.gz
cd $ompi_prefix
./configure
make -j
make install
make clean
cd ..
rm -rf $ompi_prefix