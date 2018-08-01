#!/bin/bash

mkdir -v ./install
cmake -DCMAKE_INSTALL_PREFIX="./install" ../src
