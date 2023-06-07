#!/bin/zsh



# clean subfolders and build
DIR_BUILD=build
DIR_CMAKE_DEBUG=cmake-build-debug
DIR_CMAKE_RELEASE=cmake-build-release
DIR_DATA=data
DIR_PLOT=plot
DIR_PARAVIEW=ParaView

if [ -d "$DIR_BUILD" ]; then
    echo "directory $DIR_BUILD  got deleted."
    rm -rf $DIR_BUILD
fi

if [ -d "$DIR_CMAKE_DEBUG" ]; then
    echo "directory $DIR_CMAKE_DEBUG  got deleted."
    rm -rf $DIR_CMAKE_DEBUG
fi

if [ -d "$DIR_CMAKE_RELEASE" ]; then
    echo "directory $DIR_CMAKE_RELEASE  got deleted."
    rm -rf $DIR_CMAKE_RELEASE
fi

if [ -d "$DIR_DATA" ]; then
    echo "directory $DIR_DATA  got deleted."
    rm -rf $DIR_DATA
fi

if [ -d "$DIR_PLOT" ]; then
    echo "directory $DIR_PLOT  got deleted."
    rm -rf $DIR_PLOT
fi

if [ -d "$DIR_PARAVIEW" ]; then
    echo "directory $DIR_PARAVIEW  got deleted."
    rm -rf $DIR_PARAVIEW
fi

echo "cleaned"
