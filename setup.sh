#!/bin/bash

cwd=$(pwd)

# Resolve dependencies

echo 'Pip installables (scipy, numpy, pyabpoa, mappy, Cython, tqdm)'
python3 -m pip install --user --upgrade scipy numpy pyabpoa mappy Cython tqdm

echo 'conk'
python3 -m pip install --user --upgrade wheel setuptools Cython
git clone https://github.com/rvolden/conk
cd conk && make
cd $cwd

echo 'Racon'
git clone --recursive https://github.com/isovic/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd $cwd

echo 'Done'
