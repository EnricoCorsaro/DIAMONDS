Installation
============
The following installation guide is suited for Mac OS X users and may be subject to some variations for Linux-based OS, mostly concerning the compilation of the source files. If you encounter technical problems in installing the code, please send us an `e-mail <mailto:enrico.corsaro@inaf.it>`_ specifying the error messages (possibly sending your compilation log file) and the commands you attempted to execute.

Based on the feedback already reported by some of the users of Diamonds and the experience gained during the installation of the code in different operative systems, we provide a list of possible compilation issues and how to solve them at the bottom of this page.


Prerequisites
^^^^^^^^^^^^^

1. We advise you to first read the original code paper describing the code, `E. Corsaro & J. De Ridder 2014 A&A, 571, 71 <https://www.aanda.org/articles/aa/abs/2014/11/aa24181-14/aa24181-14.html>`_. This is a preliminary task that is necessary to become more familiar with the working scheme of the code and the meaning of its configuring parameters.

2. Before installing DIAMONDS you need to install the `CMake <http://www.cmake.org/>`_ compiler, a compiler suited for C, C++ source files that is able to recognize the most suited compiler installed in your machine, depending on the platform you have. For Mac OS X it is clang, while for Linux-based OS it is gcc. For our purposes, we recommend you to install CMake 2.8 or later. You can find the dmg file of the version 2.8.12.2 (suggested for Mavericks OS) `here <http://www.cmake.org/files/v2.8/cmake-2.8.12.2-Darwin64-universal.dmg>`_, while more recent versions are required for El Captain OS or more recent OS X (we recommend CMake version 3.5.1 or later in this case). 

    .. warning:: 
        Make sure you install the CMake command line tool as well, since you need that to compile DIAMONDS via terminal. To do so, either open the CMake app and go to Tools/Install for Command Line Use if you have installed it already, or select the option during the installation phase. To avoid further compilation issues, we also recommend to update Xcode to its latest version.

You also have the possibility to install cmake directly from the terminal. If you are running on a Mac OS X system, then execute the command

.. code:: shell
    
    $ sudo brew install cmake

If you are running on a Unix system such as Ubuntu, then use the command

.. code:: shell

    $ sudo apt-get install cmake

Alternatively, CMake can be installed automatically with the pipeline by using the installing shell script, ``install_osx.sh`` for Mac OS X, or ``install_unix.sh`` for Unix OS, that is provided in the GitHub repository of FAMED.

3. Retrieve the code package from the public GitHub repository. How to retrieve the package and a description of the content of the package are presented in the :ref:`package_content` section of this website. We also recommend to read this information before proceeding.


Shell script installation
^^^^^^^^^^^^^^^^^^^^^^^^^
If you decide to perform a shell script installation then you need to execute the shell script ``install_osx.sh`` for Mac OS X, or ``install_unix.sh`` for Unix OS. The shell scripts are made different depending on which OS you are running because different compilers are used and because Unix systems may require some additional fixes during the installation process (see also the section below). The scripts are available in the GitHub repository of DIAMONDS, for `Mac OS X <https://github.com/EnricoCorsaro/DIAMONDS/blob/master/install_osx.sh>`_ and for `Unix OS <https://github.com/EnricoCorsaro/DIAMONDS/blob/master/install_unix.sh>`_. Once you downloaded the script, place it under the main folder where you want all the software installed. Then we recommend to make it an executable by typing the terminal command (e.g. for the Mac OS version)

.. code:: shell
    
    $ chmod +x install_osx.sh

In order to start the installation from scratch, go to the directory where you want to place the DIAMONDS software and run the following command via terminal

.. code:: shell
    
    $ ./install_osx.sh

The shell script will also compile one demo of DIAMONDS and make a run test for it. If you are able to see the demo test running, then your code has been installed successfully.

Mac OS X 10.6 or later
^^^^^^^^^^^^^^^^^^^^^^
If one wants to follow a standard manual installation, the procedure is rather simple. However, we recommend you to strictly follow the steps listed below in the same order as they are indicated.
Once you have downloaded the DIAMONDS package and installed CMake following step #2 of the prerequisites, we can proceed by compiling the code. To do so, go to your Diamonds directory, then open the ASCII file CMakeLists.txt. Inside you will find the first lines commented with what you need to execute via terminal.
Simply follow the guidline below. Starting from your Diamonds directory:

.. code:: shell
    
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make -j 4

The last make option, specifies the number of jobs to run simultaneously (most suited for a 4-CPU hardware) and allows to speed up the compilation process.
When the compilation and linking is executed, after a while you should display the ending message

.. code:: shell

    $ Linking CXX shared library libdiamonds.dylib
    $ [100%] Built target diamonds

This means that you have successfully installed DIAMONDS as a library in your machine.
At this stage, you can run the demo files provided in the folder demo of the code package. To do so, simply open the ``.cpp`` files, copy the command line provided at the top, following the word "**clang:**", then copy the line in your terminal under the same demos folder and execute it.

The command line provided inside the demo source file is given by

.. code:: bash

    $ clang++ -o demoFileName demoFileName.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register

This means that the compilation is done by linking the library of the DIAMONDS code that you just created. The option ``-Wno-deprecated-register`` allows to silent the warning for the unspecified keyword register, no longer available in more recent versions of Xcode. You may want to remove this keyword if the warning is not present.

Once the demo is compiled, you can simply run the executable via terminal by typing

.. code:: bash
    
    $ ./demoFileName

where ``demoFileName`` is the name of the corresponding demo tha was compiled, as displayed after the compilation.


Linux OS
^^^^^^^^
The installation procedure for Linux OS is exactly the same as that provided for Mac OS X. However, when installing the code on a Linux-based environment we have experienced some issues with the local compiler used by CMake. The g++ should work fine in general, though listing a series of warnings related to missing typedefs used in the Eigen library provided within the code package.

After the compilation of DIAMONDS you can compile the demo files using the following command line:

.. code:: bash
    
    $ g++ -o demoFileName demoFileName.cpp -L../build/ -I../include/ -ldiamonds -std=c++11

If you are using GNU 4.6 or older, it will be approriate to replace the option ``-std=c++11`` with ``-std=c++0x``.
You can find a list of possible compilation problems and relative solutions listed in the following.

 
Missing library path
""""""""""""""""""""

When attempting to compile the demo files, the local path of the DIAMONDS library may not be recognized. If this happens, follow the guidelines below.
If the library cannot be found, the following error will occur:

.. code:: bash
    
    $ ./demoFileName
    ./demoFileName: error while loading shared libraries:
    libdiamonds.so.0: cannot open shared object file: No such file or directory

In a simple way, and since the DIAMONDS compilation is only required once, to avoid this error you can define the shell variable ``LD_LIBRARY_PATH`` to include the directory where the library is installed. For example, in the Bourne shell (``/bin/sh`` or ``/bin/bash``), the library search path can be set with the following commands:

.. code:: bash
    
    $ LD_LIBRARY_PATH=/localPath/Diamonds/build
    $ export LD_LIBRARY_PATH
    $ ./demoFileName

Alternatively you can set an environment variable and store the information in your bashrc file, so that you automatically load it and you don't have to set the library path everytime you reboot your system.

.. code:: bash
    
    export LD_LIBRARY_PATH=/localPath/Diamonds/build


Compilation failure due to hidden files starting with ._
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
As reported by some Linux users that have installed the code, another problem that may arise in the compiling phase is the presence of Mac OS X hidden files starting with ``._`` and ending with ``.cpp``, present in the source folder. These files are not meaningful in Linux OS as they are used by Mac OS to store information about tagging and comments, and must be removed in order to compile the code. If one of such a file, e.g. ``._HiddenFileName.cpp``, is present in your code directory, this will give rise to a bunch of error messages of the form:

.. code:: bash

    $ /localPath/Diamonds/source/._HiddenFileName.cpp:1:1: warning: null character(s) ignored [enabled by default]
    /localPath/Diamonds/source/._HiddenFileName.cpp:1:2: error: stray ‘\5’ in program
    /localPath/Diamonds/source/._HiddenFileName.cpp:1:2: error: stray ‘\26’ in program
    /localPath/Diamonds/source/._HiddenFileName.cpp:1:2: error: stray ‘\7’ in program

Make sure you have deleted them from your code folder. Then redo the compilation process from the beginning.


Compilation failure due to conflicts with existing MESA libraries
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

As reported by users that have installed the code in Ubuntu OS, the version of GNU, 4.9.X (or later) used during the compilation of Diamonds, does not allow to compile the demos provided in the package. This may occur for users that have installed MESA libraries in their system, thus generating conflicts in the call for the compiler used.
In particular, when attempting to compile a demo, one could display an error message similar to the following one:

.. code:: bash

    $ g++ -o demoFileName demoFileName.cpp -L../build/ -I../include/ -ldiamonds -std=c++11
    In file included from /localPathGNU/include/c++/4.9.3/bits/localefwd.h:40:0,
                     from /localPathGNU/include/c++/4.9.3/ios:41,
                     from /localPathGNU/include/c++/4.9.3/ostream:38,
                     from /localPathGNU/include/c++/4.9.3/iostream:39,
                     from demoFileName.cpp:6:
    /localPathGNU/include/c++/4.9.3/x86_64-pc-linux-gnu/bit
    /c++locale.h:52:23: error: 'uselocale' was not declared in this scope
        extern "C" __typeof(uselocale) __uselocale;
                       ^

In order to get rid of the problem it is necessary to force GNU to use the version 4.8 (or later). To do so we recommend to follow the steps below.
 
1. Restart the standard compilation procedure of Diamonds by using the following line commands (make sure you have first deleted or emptied the build folder)

.. code:: bash
    
    $ mkdir build
    $ cd build
    $ cmake -D CMAKE_CXX_COMPILER=g++-4.8 ..
    $ make -j 4
 
2. Go to the demos folder and compile the demos using the command line

.. code:: bash
    
    $ g++-4.8 -o demoFileName demoFileName.cpp -L../build/ -I../include/ -ldiamonds -std=c++11


Compilation error due to old assembler version for AMD chips
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A less common error in the compilation phase may arise for users running a Unix system in AMD chips. If the version of the assembler is too old, this may generate an error of the following type

.. code:: shell
    
    $ Assembler messages:
    $ 1316: Error: expecting string instruction after `rep'
 
The problem is that the GNU compiler is generating ``rep; ret`` instructions to avoid a performance penalty for AMD chips. Older assemblers detect this as an error.
A version of Binutils that causes the problem is the GNU assembler (Linux/GNU Binutils) 2.22.52. To solve the problem, it is necessary to install a more recent version of Binutils, namely the 2.23.52 (or later).


Windows OS 10
^^^^^^^^^^^^^
For Windows OS 10 we recommend using the free application for creating an Ubuntu virtual machine. For details on how to set up this environment, visit `Install Ubuntu on Windows 10 <https://ubuntu.com/tutorials/tutorial-ubuntu-on-windows#1-overview>`_. 

Once the Ubuntu VM is installed and running in Windows OS, simply follow the guidlines presented in the Linux OS section of this page. You can even decide to use the shell script installation with the ``install_unix.sh`` script inside the Ubuntu VM, making sure to have the basic ubuntu packages installed, which include the GCC compiler suite.
