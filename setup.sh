#!/bin/bash

cwd=$(pwd)

# Resolve dependencies
echo 'Racon'
git clone --recursive https://github.com/isovic/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd $cwd

echo 'poaV2'
git clone https://github.com/tanghaibao/bio-pipeline.git
cd bio-pipeline/poaV2/
make poa
cd $cwd

echo 'minimap2'
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd $cwd

echo 'water'
wget https://users.soe.ucsc.edu/~rvolden/C3POa/EMBOSS-6.6.0_v8.tar.gz
tar -xzf EMBOSS-6.6.0_v8.tar.gz
cd EMBOSS-6.6.0
./configure --prefix=/usr/local/emboss
make
cd $cwd
