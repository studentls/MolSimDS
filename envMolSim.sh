echo 'loading log4cxx.'
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/software/log4cxx/apache-log4cxx-0.10.0/lib
export CPLUS_INCLUDE_PATH=$HOME/software/log4cxx/apache-log4cxx-0.10.0/include/:$CPLUS_INCLUDE_PATH

echo 'loading cppunit.'
export LD_LIBRARY_PATH=$HOME/software/cppunit/cppunit/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$HOME/software/cppunit/cppunit/include/:$CPLUS_INCLUDE_PATH

echo 'loading libxsd'
export CPLUS_INCLUDE_PATH=$HOME/software/xsd/xsd-3.3.0-x86_64-linux-gnu/libxsd:$CPLUS_INCLUDE_PATH

echo 'loading xerces.'
module load xerces/3.1
export CPLUS_INCLUDE_PATH=/lrz/sys/libraries/xerces/3.1/include/:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lrz/sys/libraries/xerces/3.1/lib

