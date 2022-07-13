# ## BUILD DOCKER AUTOMATICLY - not working right now.
# docker build -t rasqual:v0.0.0 -f Dockerfile .
# # v0.0.0
# docker tag rasqual:v0.0.0 ndatth/rasqual:v0.0.0
# docker push ndatth/rasqual:v0.0.0

# it is noted that usesing Dockerfile will results in erorrs in CLAPACK installation - I don't know why and I don't know how to fix it
docker run -it -h datn --name ras nfcore/base:2.1

apt-get update && apt-get install -y \
    make \
    gcc \
    wget \
    tar \
    libz-dev \
    liblapack-dev \
    libgsl-dev



# CLAPACK install
wget http://www.netlib.org/clapack/clapack.tgz -O clapack-3.2.1.tgz

tar zxvf clapack-3.2.1.tgz
cd /CLAPACK-3.2.1
mv make.inc.example make.inc
# fix errors
ulimit -s 100000 
make
## link file
ln -s lapack_LINUX.a liblapack.a
ln -s tmglib_LINUX.a libtmglib.a
ln -s blas_LINUX.a libblas.a
cd /

# gsl
wget https://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz -O gsl-2.5.tar.gz

tar zxvf gsl-2.5.tar.gz
cd /gsl-2.5
./configure --prefix=$PWD
make
make install
cd /

## add libraries
export CFLAGS="-I/CLAPACK-3.2.1/INCLUDE -I/CLAPACK-3.2.1/F2CLIBS -I/gsl-2.5"
export LDFLAGS="-L/CLAPACK-3.2.1 -L/CLAPACK-3.2.1/F2CLIBS -L/gsl-2.5/lib"
export LD_LIBRARY_PATH="/gsl-2.5/lib":$LD_LIBRARY_PATH

# rasqual
git clone https://github.com/natsuhiko/rasqual.git /rasqual
cd /rasqual/src/
make
#make install
cd /rasqual/src/ASVCF/
chmod +x *
make

cd /
# export PATH=/rasqual/src:$PATH
# export PATH=/rasqual/src/ASVCF:$PATH

# environment for running supporting scripts
wget https://raw.githubusercontent.com/datngu/nf-rasqual/main/container/environment.yml -O /environment.yml
conda install mamba -n base -c conda-forge -y
mamba env create -f /environment.yml
conda activate rasqual
#export PATH=/opt/conda/envs/rasqual/bin:$PATH
ln -s /rasqual/src/rasqual /bin/rasqual
ln -s /rasqual/src/rasqual /bin/createASVCF.sh


######################################
# EXPORT IMAGE
######################################
docker commit --change "LABEL Docker image containing all requirements for running RASQUAL" ras ndatth/rasqual:v0.0.0
# docker tag rasqual:v0.0.0 ndatth/rasqual:v0.0.0
docker push ndatth/rasqual:v0.0.0

######################################
# TESTING docker
######################################

docker run --rm -v $PWD:/data ndatth/rasqual:v0.0.0 rasqual

docker run -it --rm -v $PWD:/data ndatth/rasqual:v0.0.0 createASVCF.sh
