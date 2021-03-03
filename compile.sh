#!/bin/bash

export G4WORKDIR=$PWD

if [ ! -d build ]
then
  mkdir build
fi

cd build

rm -rf *

#source ~/opt/install/geant4-10.02.p03/share/Geant4-10.2.3/geant4make/geant4make.sh
cmake -DGeant4_DIR=$G4LIB -DCMAKE_INSTALL_PREFIX=$G4WORKDIR ../

make -j2

make install

cd ../

#./bin/Darwin-clang/Proton_MC 100000 $RANDOM Slab 300 10
