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
  gdb \
  git \
  openjdk-17-jre \
  liblzma-dev \
  libcurl4-gnutls-dev \
  libbz2-dev \
  libhts-dev \
  libssl-dev \
  libjpeg-dev \
  pkg-config \
  python-dev \
  python3 \
  python3-pip \
  python3-setuptools \
  tabix \
  libvcflib-tools \
  wget \
  zlib1g-dev


if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

pip3 install 'cluster_vcf_records==0.13.3'

#_________________________ bcftools _______________________#
cd $install_root
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xf bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2/
make
cd ..
cp -s bcftools-1.10.2/bcftools .

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


#______________________ gramtools _____________________________#
# Why six>=1.14.0?
# See https://github.com/pypa/virtualenv/issues/1551
pip3 install tox "six>=1.14.0"
cd $install_root
git clone https://github.com/iqbal-lab-org/gramtools
cd gramtools
git checkout 04c4ba717399507b643fd4b77a61c048ef2ed83f
# Note: a simple "pip3 install ." works for singularity but
# not for docker - the `gram` exectuable does not get
# put where gramtools expects to find it. The method
# below, which explicitly builds the binary, then installs
# does work ok for both docker and singularity.
mkdir cmake-build
cd cmake-build
cmake .. -DCMAKE_BUILD_TYPE=REL_WITH_ASSERTS
make gram
cd ..
pip3 install -e .
rm -rf cmake-build .git

#______________________ ivcmerge ______________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/ivcfmerge.git
cd ivcfmerge
git checkout 5819787614a263a9f35fd0c247442f092ab174ff
pip3 install .

