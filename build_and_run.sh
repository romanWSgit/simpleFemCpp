#!/bin/zsh

sh clean.sh

mkdir doc_doxygen
echo "cmake:"
mkdir build && cd build
cmake ..
make

EXECUTABLE=simpleFemCpp
if [ -f "$EXECUTABLE" ]; then
    echo "$EXECUTABLE exists."
    echo "run ..." 
    ./$EXECUTABLE
else 
    echo "$EXECUTABLE does not exist."
fi
cd ..
echo "created data files, plots and pvd file:"

cd plot
open plot.svg
open plot2.svg
cd ../ParaView
open tension_rod.pvd


