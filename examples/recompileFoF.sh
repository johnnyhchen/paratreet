cd ../../unionfind
make clean
make
cd ../paratreet/src
make clean
make -j8
cd ../examples
make clean
make FoFApp

