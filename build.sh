make clean
sh find.sh
sh xsd.sh
echo 'loading gcc 4.7'
module load gcc/4.7
echo 'loading xerces'
module load xerces/3.1
export export CPLUS_INCLUDE_PATH=/lrz/sys/libraries/xerces/3.1/include/:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lrz/sys/libraries/xerces/3.1/lib
echo 'configuring include/lib dirs'
export CPLUS_INCLUDE_PATH=$HOME/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib
make
