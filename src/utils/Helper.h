//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File Helper.cpp
// contains some basic helper functions
//------------------------------------------------------------------------------------------------

#ifndef HELPER_HEADER_
#define HELPER_HEADER_

#include "Base.h"

#include <iostream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>

/// specifies line width used by line(...)
static const int g_linewidth = 40;

namespace utils
{
	///
	/// prints a line
	/// @param character character to use for line
	///
	inline void	line(unsigned char character = '-')
	{
		for(int i = 0; i < g_linewidth; i++)std::cout<<character;
		std::cout<<std::endl;
	}

	///
	/// test if file exists
	/// @param filename valid path of file
	/// @return true returns if file exists
	///
	inline bool fileExists(const char *filename)
	{
		//testwise opening
		FILE *pFile = NULL;

		pFile = fopen(filename, "rb");
		
		if(pFile)fclose(pFile);
	
		//return pFile != NULL 
		return pFile != NULL;
	}

	///
	/// check if str has a valid number format
	/// allowed formats are 4 ; 4. ; 0.4 ; .5  ; .
	/// only one point . is allowed, delimiters like 400,000.0 are not allowed
	///
	/// @param str
	/// @return if str has format 4 ; 4. ; 0.4; .5; .
	///
	inline bool	strIsNumber(const char *str)
	{
		bool bPointFound = false;

		//str has invalid format if something except ' ', '.' or 0-9 is found
		const char *ptr = str;
		while(*ptr != '\0')
		{
			// 0-9 are ASCII cods 0x30 to 0x39
			if((*str < 0x30 || *str > 0x39) && *str != ' ' && *str != '.')return false;
			
			if(bPointFound)return false;
			
			if(*str == '.')bPointFound = true;

			ptr++;
		}

		return true;
	}

	/// gets file extension, e.g. for example.xml the string "xml" is returned
	/// @param str input filename
	/// @return file extension of str without dot
	inline char*		getFileExtension(char *str)
	{
		int iLastDot = -1;
		int	iChar = strlen(str);

		//letzten Backslash suchen
		while(iChar != 0)
		{
			if(str[iChar] == '.')iLastDot = iChar;
			iChar--;
		}

		//String zurückgeben
		return str + iLastDot + 1;
	}

}

#endif
