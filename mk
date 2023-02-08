#!/bin/sh
cd bin
#make clean
#/bin/cp Makefile_r4 Makefile
#/bin/cp Makefile_r8_g Makefile
make -f Makefile_r8_g
cp gf.x ../ 
cd ..
time ./gf.x
