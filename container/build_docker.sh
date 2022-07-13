## BUILD DOCKER AUTOMATICLY

docker build -t rasqual:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag rasqual:v0.0.0 ndatth/rasqual:v0.0.0
docker push ndatth/rasqual:v0.0.0


## testing
#docker run -it --rm -h datn --name ras ubuntu:18.04
docker run -it -v /home/datn/github/nf-rasqual:/host/ -h datn --name ras ubuntu:18.04


##################
#FROM nfcore/base:2.1
docker run -it -h datn --name ras nfcore/base:2.1
#LABEL authors="Dat T Nguyen - ndat<at>utexas.edu" \
      description="Docker image containing all requirements for running RASQUAL" 


#RUN 
apt-get update && apt-get install -y \
    make \
    gcc \
    wget \
    tar \
    libz-dev

# crtl p q

# copy libs
docker cp lib_src/clapack-3.2.1.tgz ras:clapack-3.2.1.tgz
docker cp lib_src/gsl-2.5.tar.gz ras:/gsl-2.5.tar.gz
docker cp lib_src/rasqual_kauralasoo.tgz ras:/rasqual.tgz

docker attach ras

# CLAPACK
#RUN 
tar zxvf clapack-3.2.1.tgz
cd /CLAPACK-3.2.1
mv make.inc.example make.inc
#RUN
ulimit -s 100000 
make
## link file
#RUN 
ln -s lapack_LINUX.a liblapack.a
#RUN 
ln -s tmglib_LINUX.a libtmglib.a
#RUN 
ln -s blas_LINUX.a libblas.a
cd /

# gsl
#RUN 
tar zxvf gsl-2.5.tar.gz
cd /gsl-2.5
#RUN 
./configure --prefix=$PWD
#RUN 
make
#RUN 
make install
cd /

# RASQUAL
export CFLAGS="-I/CLAPACK-3.2.1/INCLUDE -I/CLAPACK-3.2.1/F2CLIBS -I/gsl-2.5"
export LDFLAGS="-L/CLAPACK-3.2.1 -L/CLAPACK-3.2.1/F2CLIBS -L/gsl-2.5/lib"
export LD_LIBRARY_PATH="/gsl-2.5/lib":$LD_LIBRARY_PATH

export CFLAGS="-I/CLAPACK-3.2.1/INCLUDE -I/CLAPACK-3.2.1/F2CLIBS -I/gsl-2.5"
export LDFLAGS="-L/CLAPACK-3.2.1 -L/CLAPACK-3.2.1/F2CLIBS -L/gsl-2.5/lib"
export LD_LIBRARY_PATH=/gsl-2.5/lib:$LD_LIBRARY_PATH

tar zxvf rasqual.tgz
cd /rasqual/src
make
make install
cd /rasqual/src/ASVCF
make
make install

cd /
export PATH=/rasqual/src:$PATH
export PATH=/rasqual/src/ASVCF:$PATH

# environment for running supporting scripts
ADD environment.yml /
RUN conda install mamba -n base -c conda-forge -y
RUN mamba env create -f /environment.yml
ENV PATH /opt/conda/envs/rasqual/bin:$PATH


























##################### install 
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