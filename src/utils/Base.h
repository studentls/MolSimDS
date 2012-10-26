//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File MolSim.cpp
// contains main class MolSim
//------------------------------------------------------------------------------------------------

#ifndef BASE_HEADER_
#define BASE_HEADER_


//C++ defined?
#ifndef __cplusplus
#error use c++ compiler to succeed!
#endif

//defines for debug state should be later included in a separate Base.h
#ifndef DEBUG


//MSVC
#ifdef _MSC_VER 

#ifdef _DEBUG
#define DEBUG
#endif

#endif

//GCC
#ifdef __GNUC__

//doesn't define an internal DEBUG state, so it should be invoked from withhin the command line
//add -DDEBUG as an extra flag

#endif

//Intel CC
#ifdef __INTEL_COMPILER 

//same as withhin GCC

#endif

#endif

#ifdef DEBUG

//Basic LOG macro
//not done yet...
#define LOG(x) 

#endif


//typedef's

typedef int				err_type;
typedef unsigned char	byte;


//error macros&constants
#define FAILED(x)		( (x) < 0 ? true : false )
#define SUCCEEDED(x)	( (x) > 0 ? true : false )
#define E_UNKNOWN			-1
#define E_OUTOFMEMORY		-2
#define E_INVALIDPARAM		-3
#define E_INDEXOUTOFRANGE	-4
#define E_FILENOTFOUND		-5
#define E_NOTINITIALIZED	-6
#define S_OK				1

//delete macros
#define DELETE(x)		 {if((x))delete (x); (x) = NULL;}
#define DELETE_A(x)		 {if((x))delete [] (x); (x) = NULL;}


#endif