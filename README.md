MolSim
======

The Molecular Dynamics teaching code.


Linux
-----

to build on Linux run first sh find.sh 
and then make


Windows
-------

open project in VC and compile


Mac OS X
--------

an XCode project(Mountain Lion) is provided

ICE
--------

to compile on LRZ HPC ICE use the makefile MakefileICE. Therefore you will need to have all necessary libraries(CppUnit, Log4Cxx, (GLFW), libxsd) copied into two folders /include and /lib in your home directory. To compile without GLFW give gcc the flag -DICE. Code is compatible with GCC 4.7

-------
Note:

To run the simulator, you will need to have xerces library installed

To generate docfiles, run doxygen Doxyfile.doxy

To compile sources, ensure to have the following libraries installed: 

Log4Cxx, CppUnit, Xerces 3.1.3, OpenGL, GLFW, GLEW, boost thread


Usage
-------

molsim <file> <endtime> <delta_t>	run simulation, <file> should be in the default *.txt format, endtime specifies when simulation shall stop, <delta_t> is stepsize(so simulation will run endtime / delta_t steps
molsim -help				displays help instructions on avaliable commands
molsim -test <name>			run test case with <name>, if <name> is left blank, all test cases will be run
molsim -ptest <file>			run performance tests for given xml <file> (OpenMP tests, if compiled with OpenMP)
molsim -showtests			show avaliable test cases by name
molsim --viewer				optional: run with builtin OpenGL Viewer
