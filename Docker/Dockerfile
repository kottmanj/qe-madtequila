FROM ubuntu:18.04
MAINTAINER JSK
LABEL "About"="Madness+Tequila in one container"

USER root

WORKDIR /app

# basics
RUN apt-get update && apt-get install -y \
        build-essential \
        cmake   \
        git \
        libgtk2.0-dev \
        pkg-config \
        libswscale-dev\
        liblapack-dev \
        libatlas-base-dev \
        libgomp1 \
        python3 \
        python3-pip \
        mpich \
        wget

RUN pip3 install mkl


# install madness
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
RUN apt-get update
RUN apt-get install -y intel-mkl-64bit-2020.1-102
RUN export MKLROOT=/opt/intel/compilers_and_libraries/linux/mkl

RUN wget https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2
RUN tar --bzip2 -xf boost_1_73_0.tar.bz2

RUN git clone https://github.com/dpilger26/NumCpp.git numcpp

RUN git clone https://github.com/kottmanj/madness.git madsrc
RUN mkdir madroot

WORKDIR /madroot/
RUN cmake -D ENABLE_MKL=ON -D MKL_ROOT_DIR=/opt/intel/compilers_and_libraries/linux/mkl -D CMAKE_CXX_FLAGS='-O3 -DNDEBUG -march=native -I/numcpp/include/ -I/boost_1_73_0/'  /madsrc/
RUN make

RUN ln -s /madroot/src/apps/pno/pno_integrals /usr/bin/pno_integrals
RUN export MAD_ROOT_DIR=/madroot/

# install tequila
RUN apt update -y
RUN apt install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa -y
RUN apt install python3.7 -y
RUN python3.7 --version
RUN python3.7 -m pip --version 
RUN python3.7 -m pip install --upgrade pip
RUN python3.7 -m pip --version
RUN rm /usr/bin/python
RUN ln -s /usr/bin/python3.7 /usr/bin/python
RUN python --version 
RUN python3.7 -m pip install git+https://github.com/aspuru-guzik-group/tequila.git

ENTRYPOINT bash
