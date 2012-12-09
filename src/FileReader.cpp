/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "utils/Vector.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

using namespace std;

FileReader::FileReader() {
}

FileReader::~FileReader() {
}


void FileReader::readFile(std::vector<Particle>& particles, const char* filename) {
	double x[] = {0,0,0};
	double v[] = {1,1,1};
	double m = 1;
	double epsilon = 0;
	double sigma = 0;
    int num_particles = 0;

    std::ifstream input_file(filename);
    string tmp_string;

    if (input_file.is_open()) {

    	getline(input_file, tmp_string);
    	LOG4CXX_INFO(fileReaderLogger, "Read line: " << tmp_string);

    	while (tmp_string.size() == 0 || tmp_string[0] == '#') {
    		getline(input_file, tmp_string);
    		LOG4CXX_INFO(fileReaderLogger, "Read line: " << tmp_string);
    	}

    	istringstream numstream(tmp_string);
    	numstream >> num_particles;
    	LOG4CXX_INFO(fileReaderLogger, "Reading " << num_particles << ".");
    	getline(input_file, tmp_string);
    	LOG4CXX_INFO(fileReaderLogger, "Read line: " << tmp_string);

    	for (int i = 0; i < num_particles; i++) {
    		istringstream datastream(tmp_string);

    		for (int j = 0; j < 3; j++) {
    			datastream >> x[j];

    		}
    		for (int j = 0; j < 3; j++) {
    			datastream >> v[j];
    		}
    		if (datastream.eof()) {
    			LOG4CXX_ERROR(fileReaderLogger, "Error reading file: eof reached unexpectedly reading from line " << i);
    			exit(-1);
    		}
    		datastream >> m;
    		datastream >> epsilon;
    		datastream >> sigma;
    		Particle p(x, v, m, epsilon, sigma);
    		particles.push_back(p);

    		getline(input_file, tmp_string);
    		cout << "Read line: " << tmp_string << endl;
    	}
    } else {
    	LOG4CXX_ERROR(fileReaderLogger, "Error: could not open file " << filename);
    	exit(-1);
    }

}

