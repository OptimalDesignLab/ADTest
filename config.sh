#!/bin/bash

mkdir -v ./install
cmake -DCMAKE_INSTALL_PREFIX="./install" \
  -DCMAKE_CXX_FLAGS="-O0 -g" \
  ../src
