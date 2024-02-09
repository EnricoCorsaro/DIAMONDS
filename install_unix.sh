#!/bin/bash
if ! [ -x "$(command -v git)" ]; then
	echo "Error: git is not installed. Aborting..." >&2
	exit 1
fi

if ! [ -x "$(command -v cmake)" ]; then
	echo "Cmake is not installed. Trying to install it using apt-get..." >&2

	if ! [ -x "$(command -v apt-get)" ]; then
		echo "Error: apt-get is not installed. Aborting..." >&2
		exit 1
	else
		sudo apt-get install cmake
	fi
fi

if [ $flag1 -eq 1 ]; then
	echo "-----------------------------------"
	echo " Cloning and installing DIAMONDS..."
	echo "-----------------------------------"
	git clone https://github.com/EnricoCorsaro/DIAMONDS.git
	mkdir DIAMONDS/build
	cd DIAMONDS/build/
	cmake ..
	make -j 4
	echo "-----------------------------------"
	echo " Compiling and running test demo..."
	echo "-----------------------------------"
	cd ../../
	export LD_LIBRARY_PATH=${PWD}/DIAMONDS/build
	cd DIAMONDS/demos/
	g++ -o demoSingle2DGaussian demoSingle2DGaussian.cpp -L../build/ -I../include/ -ldiamonds -std=c++11
	./demoSingle2DGaussian
	cd ../../
	echo " "
	echo " "	
fi
