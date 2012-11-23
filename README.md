MolSim
===

The Molecular Dynamics teaching code.


Linux
---
to build on Linux run first sh find.sh 
and then make


Windows
---
open project in VC and compile




---
Note:
You have to have xerces library installed

To generate docfiles, run doxygen Doxyfile.doxy


Usage
---
molsim <file> <endtime> <delta_t>	run simulation, <file> should be in the default *.txt format, endtime specifies when simulation shall stop, <delta_t> is stepsize(so simulation will run endtime / delta_t steps
molsim -help				displays help instructions on avaliable commands
molsim -test <name>			run test case with <name>, if <name> is left blank, all test cases will be run
molsim -showtests			show avaliable test cases by name