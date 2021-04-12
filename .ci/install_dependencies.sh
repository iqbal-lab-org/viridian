#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  cmake \
  automake \
  gcc \
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
  samtools

if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

#________________________ minimap2 __________________________#
cd $install_root
git clone https://github.com/lh3/minimap2.git minimap2_git
cd minimap2_git
git checkout 4dfd495cc2816f67556bc6318654e572636ee40a
make
cd ..
cp -s minimap2_git/minimap2 .
cp -s minimap2_git/misc/paftools.js .
wget https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2
tar -jxvf k8-0.2.4.tar.bz2
cp k8-0.2.4/k8-`uname -s` k8

#________________________ racon _____________________________#
git clone --recursive https://github.com/lbcb-sci/racon.git racon-git
cd racon-git
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ../../
cp -s racon-git/build/bin/racon .

#________________________ viridian __________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/viridian.git
cd viridian
git checkout 1bcb9f2df0130ae6d6873ad78c385ef03642b9c0
pip3 install .

#________________________ mummer ____________________________#
cd $install_root
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xvf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure --prefix=/bioinf-tools/
make
make install
cd ..

#________________________ varifier __________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/varifier.git
cd varifier
git checkout 368079a27c36784a01a10ec17b3a2777b8013629
pip3 install .
cd ..

#________________________ qcovid ____________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/QCovid.git
cd QCovid
git checkout 598ef7f9c6e267561021888abbe90727f2784188
pip3 install .
cd ..
