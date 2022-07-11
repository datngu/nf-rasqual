## BUILD DOCKER AUTOMATICLY

docker build -t rasqual:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag rasqual:v0.0.0 ndatth/rasqual:v0.0.0
docker push ndatth/rasqual:v0.0.0




## testing
#docker run -it --rm -h datn --name ras ubuntu:18.04
docker run -it -v /home/datn/github/nf-rasqual:/host/ -h datn --name ras ubuntu:18.04

# install 
apt-get update && apt-get install -y \
    make \
    gcc \
    wget \
    git \
    libz-dev


wget http://www.netlib.org/clapack/clapack.tgz clapack-3.2.1.tgz
tar zxvf clapack-3.2.1.tgz
cd CLAPACK-3.2.1
mv make.inc.example make.inc
make
# link file
ln -s lapack_LINUX.a liblapack.a
ln -s tmglib_LINUX.a libtmglib.a
ln -s blas_LINUX.a libblas.a

cd /
# gsl
wget https://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz gsl-2.5.tar.gz

tar zxvf gsl-2.5.tar.gz
cd gsl-2.5
./configure --prefix=$PWD
make
make install

cd /

git clone https://github.com/kauralasoo/rasqual.git /rasqual

## build rasqual

RASQUALDIR=/rasqual
cd $RASQUALDIR/src
# Not run!  Please export your environment.
export CFLAGS="-I/CLAPACK-3.2.1/INCLUDE -I/CLAPACK-3.2.1/F2CLIBS -I/gsl-2.5"
export LDFLAGS="-L/CLAPACK-3.2.1 -L/CLAPACK-3.2.1/F2CLIBS -L/gsl-2.5/lib"
export LD_LIBRARY_PATH=/gsl-2.5/lib:$LD_LIBRARY_PATH
make
make install