FROM nfcore/base:2.1

LABEL authors="Dat T Nguyen - ndat<at>utexas.edu" \
      description="Docker image containing all requirements for running RASQUAL" 


RUN apt-get update && apt-get install -y \
    make \
    gcc \
    wget \
    tar \
    libz-dev \
    liblapack-dev \
    libgsl-dev

# copy libs
# ADD lib_src/clapack-3.2.1.tgz /clapack-3.2.1.tgz
# ADD lib_src/gsl-2.5.tar.gz /gsl-2.5.tar.gz
# ADD lib_src/rasqual_kauralasoo.tgz /rasqual.tgz
#RUN wget http://www.netlib.org/clapack/clapack.tgz -O clapack-3.2.1.tgz
#RUN wget https://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz -O gsl-2.5.tar.gz
RUN git clone https://github.com/kauralasoo/rasqual.git /rasqual

# CLAPACK
#RUN tar -xzf clapack-3.2.1.tgz
#WORKDIR /CLAPACK-3.2.1
#RUN ln -s make.inc.example make.inc
# fix memory errors in compiling.
#RUN ulimit -s 100000
#RUN make
## link file
#RUN ln -s lapack_LINUX.a liblapack.a
#RUN ln -s tmglib_LINUX.a libtmglib.a
#RUN ln -s blas_LINUX.a libblas.a
#WORKDIR /

# gsl
#RUN tar -xzf gsl-2.5.tar.gz
#WORKDIR /gsl-2.5
#RUN ./configure --prefix=$PWD
#RUN make
#RUN make install
#WORKDIR /

# RASQUAL
#ENV CFLAGS "-I/CLAPACK-3.2.1/INCLUDE -I/CLAPACK-3.2.1/F2CLIBS -I/gsl-2.5"
#ENV LDFLAGS "-L/CLAPACK-3.2.1 -L/CLAPACK-3.2.1/F2CLIBS -L/gsl-2.5/lib"
#ENV LD_LIBRARY_PATH /gsl-2.5/lib:$LD_LIBRARY_PATH
#export CFLAGS="-I/CLAPACK-3.2.1/INCLUDE -I/CLAPACK-3.2.1/F2CLIBS -I/gsl-2.5"
#export LDFLAGS="-L/CLAPACK-3.2.1 -L/CLAPACK-3.2.1/F2CLIBS -L/gsl-2.5/lib"
#export LD_LIBRARY_PATH=/gsl-2.5/lib:$LD_LIBRARY_PATH

#RUN tar -xzf rasqual.tgz
WORKDIR /rasqual/src
RUN make
RUN make install
WORKDIR /rasqual/src/ASVCF
RUN make
RUN make install

WORKDIR /
ENV PATH /rasqual/src:$PATH
ENV PATH /rasqual/src/ASVCF:$PATH

# environment for running supporting scripts
ADD environment.yml /
RUN conda install mamba -n base -c conda-forge -y
RUN mamba env create -f /environment.yml
ENV PATH /opt/conda/envs/rasqual/bin:$PATH