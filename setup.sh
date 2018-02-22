#!/bin/bash

cwd=$(pwd)

# Resolve dependencies
echo 'Racon'
git clone https://github.com/isovic/racon.git  && cd racon && make modules && make tools && make -j
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
tar -xzf EMBOSS-6.6.0_v8.tar.gz
cd EMBOSS-6.6.0
./configure --prefix=/usr/local/emboss
make
cd $cwd
