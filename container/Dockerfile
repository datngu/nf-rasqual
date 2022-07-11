FROM ubuntu:18.04

MAINTAINER Dat T Nguyen "ndat<at>utexas.edu"
LABEL authors="Dat T Nguyen" \
      description="Docker image containing all requirements for running RASQUAL" 


RUN apt-get update && apt-get install -y \
    make \
    gcc \
    wget \
    tar \
    git \

RUN wget http://www.netlib.org/clapack/clapack.tgz /clapack.tgz
RUN tar zxvf clapack.tgz
cd 
mv make.inc.example make.inc
make




# RASQUAL
RUN git clone https://github.com/kauralasoo/rasqual.git /rasqual && \
     make -C /rasqual/src && \
     make -C /rasqual/src/ASVCF && \
     cp /rasqual/src/rasqual /usr/bin


# environment for running supporting scripts
ADD environment.yml /
RUN conda install mamba -n base -c conda-forge -y
RUN mamba env create -f /environment.yml
ENV PATH /opt/conda/envs/rasqual/bin:$PATH