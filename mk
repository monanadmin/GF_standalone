#!/bin/sh
cd bin
#make clean
#/bin/cp Makefile_r4 Makefile
#/bin/cp Makefile_r8_g Makefile
rm -f gf.x
make -f Makefile_r8
sleep 2
cp gf.x ../ 
cd ..
time ./gf.x

echo "done"