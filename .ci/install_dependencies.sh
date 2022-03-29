#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get update
apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  cmake \
  automake \
  gcc \
  gcc-10 \
  g++-10 \
  gdb \
  git \
  python3 \
  python3-pip \
  python3-setuptools \
  python3-dev \
  wget \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libhts-dev \
  tabix \
  curl \
  libvcflib-tools \
  libcurl4-gnutls-dev \
  libssl-dev \
  samtools

update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 --slave /usr/bin/g++ g++ /usr/bin/g++-10 --slave /usr/bin/gcov gcov /usr/bin/gcov-10
pip3 install tox

if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

#_________________________ bcftools _________________________#
cd $install_root
wget -q https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xf bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2/
make
cd ..
cp -s bcftools-1.10.2/bcftools .

#________________________ vt ________________________________#
cd $install_root
git clone https://github.com/atks/vt.git vt-git
cd vt-git
git checkout 2187ff6347086e38f71bd9f8ca622cd7dcfbb40c
make
cd ..
cp -s vt-git/vt .

#________________________ minimap2 __________________________#
cd $install_root
git clone https://github.com/lh3/minimap2.git minimap2_git
cd minimap2_git
git checkout 4dfd495cc2816f67556bc6318654e572636ee40a
make
cd ..
cp -s minimap2_git/minimap2 .
cp -s minimap2_git/misc/paftools.js .
wget -q https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2
tar -jxvf k8-0.2.4.tar.bz2
cp k8-0.2.4/k8-`uname -s` k8

#________________________ racon _____________________________#
cd $install_root
git clone --recursive https://github.com/lbcb-sci/racon.git racon-git
cd racon-git
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make CC=gcc-10 CPP=g++-10 CXX=g++-10 LD=g++-10
cd ../../
cp -s racon-git/build/bin/racon .

#________________________ viridian __________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/viridian.git
cd viridian
git checkout aa9a5ee8df291114637ea49a9670f28408d47311
pip3 install .

#________________________ mummer ____________________________#
cd $install_root
wget -q https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xvf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure LDFLAGS=-static
make
make install
ldconfig
cd ..

#________________________ varifier __________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/varifier.git
cd varifier
git checkout a6938f9730d21bef166472e635b48d99e765bcd0
pip3 install .
cd ..
