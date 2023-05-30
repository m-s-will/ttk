mpich=4.0.2
mpich_prefix=mpich-$mpich

require-pkgs \
	wget			\
	build-essential	\

wget https://www.mpich.org/static/downloads/$mpich/$mpich_prefix.tar.gz
tar xvzf $mpich_prefix.tar.gz
cd $mpich_prefix
./configure --disable-fortran
make -j
make install
make clean
cd ..
rm -rf $mpich_prefix