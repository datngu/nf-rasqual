FROM nfcore/base:2.1

MAINTAINER Dat T Nguyen "ndat<at>utexas.edu"
LABEL authors="Dat T Nguyen" \
      description="Docker image containing all requirements for running RASQUAL" 

ADD environment.yml /

RUN conda install mamba -n base -c conda-forge -y
RUN mamba env create -f /environment.yml
ENV PATH /opt/conda/envs/rasqual/bin:$PATH

# RASQUAL
RUN git clone https://github.com/kauralasoo/rasqual.git /rasqual && \
     make -C /rasqual/src && \
     make -C /rasqual/src/ASVCF && \
     cp /rasqual/src/rasqual /usr/bin
