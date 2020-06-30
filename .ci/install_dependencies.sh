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
  openjdk-8-jre \
  liblzma-dev \
  libcurl4-gnutls-dev \
  libbz2-dev \
  libhts-dev \
  libssl-dev \
  pkg-config \
  python3 \
  python3-pip \
  python3-setuptools \
  tabix \
  wget \
  zlib1g-dev


if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

#________________________ nextflow ____________________________#
cd $install_root
wget -qO- https://get.nextflow.io | bash
chmod 755 nextflow

#________________________ vt __________________________________#
cd $install_root
git clone https://github.com/atks/vt.git vt-git
cd vt-git
git checkout 2187ff6347086e38f71bd9f8ca622cd7dcfbb40c
make
cd ..
cp -s vt-git/vt .

#______________________vcflib _________________________________#
cd $install_root
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
make

# Why six>=1.14.0?
# See https://github.com/pypa/virtualenv/issues/1551
pip3 install tox "six>=1.14.0"

cd $install_root
git clone https://github.com/iqbal-lab-org/gramtools
cd gramtools
git checkout 9313eceb606a6fc159e4a14c168b7a6f888c5ed2
pip3 install .
