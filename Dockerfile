FROM ubuntu:16.10

MAINTAINER Dat T Nguyen "ndat<at>utexas.edu"
LABEL authors="Dat T Nguyen" \
      description="Docker image containing all requirements for running RASQUAL" 

RUN apt-get update && apt-get install -y \
    software-properties-common \
    git

RUN add-apt-repository -y universe && apt-get install -y \
    f2c \
    gcc \
    libblas-dev \
    libcblas-dev \
    libclapack-dev \
    libgsl-dev \
    liblapack-dev \
    make \
    zlib1g-dev 

RUN git clone https://github.com/kauralasoo/rasqual.git /rasqual && \
     make -C /rasqual/src && \
     make -C /rasqual/src/ASVCF && \
     cp /rasqual/src/rasqual /usr/bin






#ADD environment.yml /

#RUN conda install mamba -n base -c conda-forge -y
#RUN mamba env create -f /environment.yml
#ENV PATH /opt/conda/envs/rasqual/bin:$PATH

# RASQUAL
#RUN git clone https://github.com/kauralasoo/rasqual.git /rasqual && \
#     make -C /rasqual/src && \
#     make -C /rasqual/src/ASVCF && \
#     cp /rasqual/src/rasqual /usr/bin
