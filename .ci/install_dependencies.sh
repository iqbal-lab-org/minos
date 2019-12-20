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
  libbz2-dev \
  libhts-dev \
  pkg-config \
  python3 \
  python3-pip \
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


#_________________________ bcftools _______________________#
cd $install_root
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar xf bcftools-1.9.tar.bz2
rm bcftools-1.9.tar.bz2
cd bcftools-1.9/
make
cd ..
cp -s bcftools-1.9/bcftools .


#__________________________ BWA____________________________#
cd $install_root
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar xf bwa-0.7.17.tar.bz2
rm bwa-0.7.17.tar.bz2
cd bwa-0.7.17/
make
cd ..
cp -s bwa-0.7.17/bwa .


#________________________ mummer ____________________________#
cd $install_root
wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar xf mummer-4.0.0beta2.tar.gz
rm mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2
./configure
make
make install




pip3 install pysam matplotlib pandas seaborn pymummer cluster_vcf_records tox

cd $install_root
git clone https://github.com/iqbal-lab-org/gramtools
cd gramtools
git checkout b95321b7e15c0a574863b298475c880099ea6e5d
pip3 install .
