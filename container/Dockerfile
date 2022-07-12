FROM nfcore/base:2.1

LABEL authors="Dat T Nguyen - ndat<at>utexas.edu" \
      description="Docker image containing all requirements for running RASQUAL" 


RUN apt-get update && apt-get install -y \
    make \
    gcc \
    wget \
    tar \
    libz-dev

# copy libs
ADD lib_src/clapack-3.2.1.tgz /clapack-3.2.1.tgz
ADD lib_src/gsl-2.5.tar.gz /gsl-2.5.tar.gz
ADD lib_src/rasqual_kauralasoo.tgz /rasqual.tgz

# CLAPACK
#RUN tar zxvf clapack-3.2.1.tgz
WORKDIR /CLAPACK-3.2.1
COPY make.inc.example make.inc
RUN make
## link file
RUN ln -s lapack_LINUX.a liblapack.a
RUN ln -s tmglib_LINUX.a libtmglib.a
RUN ln -s blas_LINUX.a libblas.a
WORKDIR /

# gsl
RUN tar zxvf gsl-2.5.tar.gz
WORKDIR /gsl-2.5
RUN ./configure --prefix=$PWD
RUN make
RUN make install
WORKDIR /

# RASQUAL
ENV CFLAGS -I/CLAPACK-3.2.1/INCLUDE -I/CLAPACK-3.2.1/F2CLIBS -I/gsl-2.5:$CFLAGS
ENV LDFLAGS -L/CLAPACK-3.2.1 -L/CLAPACK-3.2.1/F2CLIBS -L/gsl-2.5/lib:$LDFLAGS
ENV LD_LIBRARY_PATH /gsl-2.5/lib:$LD_LIBRARY_PATH

RUN tar zxvf rasqual.tgz
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