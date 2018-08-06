#!/bin/bash

adeptdir="$HOME/scratch/adept-2.0.5/install" 
export CMAKE_PREFIX_PATH="$adeptdir":$CMAKE_PREFIX_PATH

mkdir -v ./install
cmake -DCMAKE_INSTALL_PREFIX="./install" \
  -DCMAKE_CXX_FLAGS="-Ofast -march=native -DNDEBUG" \
  -DENABLE_CODIPACK=True \
  -DCODIPACK_DIR="$HOME/scratch/CoDiPack" \
  -DENABLE_ADEPT=True \
  -DADEPT_DIR="$adeptdir" \
  ../src
