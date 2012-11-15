/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef FILE_READER_H_
#define FILE_READER_H_

#include "Logging.h"
#include "Particle.h"
#include <vector>

class FileReader {

public:
	FileReader();
	virtual ~FileReader();

	void readFile(std::vector<Particle>& particles, const char* filename);

};

#endif /* FILE_READER_H_ */
