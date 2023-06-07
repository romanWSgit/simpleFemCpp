#!/bin/zsh

sh clean.sh

mkdir doc_doxygen
echo "cmake:"
mkdir build && cd build
cmake ..
make
