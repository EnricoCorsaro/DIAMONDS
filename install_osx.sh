#!/bin/bash
if ! [ -x "$(command -v git)" ]; then
	echo "Error: git is not installed. Aborting..." >&2
	exit 1
fi

if ! [ -x "$(command -v cmake)" ]; then
	echo "Cmake is not installed. Trying to install it using Homebrew..." >&2
	
	if ! [ -x "$(command -v brew)" ]; then
		echo "Error: Homebrew is not installed. Aborting..." >&2
		exit 1
	else
		sudo brew install cmake	
	fi
fi

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
cd ../demos/
clang++ -o demoSingle2DGaussian demoSingle2DGaussian.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
./demoSingle2DGaussian
cd ../../
echo " "
echo " "	