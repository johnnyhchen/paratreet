cd ../../unionfind/prefixLib
make clean
make
cd ../
make clean
make
cd ../paratreet/src
make clean
make -j8
cd ../examples
make clean
make FoFApp

