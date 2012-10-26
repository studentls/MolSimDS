/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#ifndef XYZWRITER_H_
#define XYZWRITER_H_

#include "../Particle.h"
#include <fstream>
#include <vector>

namespace outputWriter {

class XYZWriter {

public:
	XYZWriter();

	virtual ~XYZWriter();

	void plotParticles(std::vector<Particle> particles, const std::string& filename, int iteration);

};

}

#endif /* XYZWRITER_H_ */
