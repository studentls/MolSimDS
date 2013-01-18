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
molsim -showtests			show avaliable test cases by name
molsim --viewer				optional: run with builtin OpenGL Viewer