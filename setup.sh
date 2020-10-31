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

echo 'pyabpoa'
pip3 install --user pyabpoa

echo 'mappy'
pip3 install --user mappy

echo 'conk'
python3 -m pip install --user --upgrade wheel setuptools Cython
git clone https://github.com/rvolden/conk
cd conk && make
cd $cwd

echo 'Done'
