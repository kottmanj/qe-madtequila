FROM ubuntu:18.04
MAINTAINER JSK
LABEL "About"="Madness+Tequila in one container"

USER root

WORKDIR /app

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
        mpich \
        wget

# install python3.7 in the same way as for the standard orquestra container
# so that we can use this container with the python3 runtime of orquestra

RUN apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get install -y python3.7 && \
    apt-get install -y python3-pip && \
    apt-get install -y python3.7-dev

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 1 && \
    update-alternatives --set python3 /usr/bin/python3.7
RUN rm /usr/bin/python
RUN ln -s /usr/bin/python3.7 /usr/bin/python


# install madness
#RUN python -m pip --upgrade install mkl
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
RUN apt-get update
RUN apt-get install -y intel-mkl-64bit-2020.1-102
RUN export MKLROOT=/opt/intel/compilers_and_libraries/linux/mkl

RUN wget https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2
RUN tar --bzip2 -xf boost_1_73_0.tar.bz2
RUN export BOOSTROOT=$(realpath boost_1_73_0)

RUN git clone https://github.com/dpilger26/NumCpp.git numcpp
RUN export NUMCPPROOT=$(realpath numcpp)

RUN git clone https://github.com/kottmanj/madness.git madsrc
RUN mkdir madroot

RUN echo "$NUMCPPROOT, $BOOSTROOT"

RUN cmake -D ENABLE_MKL=ON -D MKL_ROOT_DIR=/opt/intel/compilers_and_libraries/linux/mkl -D CMAKE_CXX_FLAGS="-O3 -DNDEBUG -march=native -I/app/numcpp/include/ -I/app/boost_1_73_0/"  madsrc -Bmadroot
RUN make -C madroot

RUN ln -s /madroot/src/apps/pno/pno_integrals /usr/bin/pno_integrals

# export the madness root dir, so that tq can find the executable automatically
RUN export MAD_ROOT_DIR=$(realpath madroot)

# install tequila
# PyPi version is fine here
RUN python -m pip install tequila-basic
# let's enable qulacs as well
# so that we can do most tasks beyond simply running madness in reasonable time
RUN python -m pip install qulacs

ENTRYPOINT bash
