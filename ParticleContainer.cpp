/*
 * ParticleContainer.cpp
 */

#include "ParticleContainer.h"

ParticleContainer::ParticleContainer(utils::Vector<Particle> particles) {
	this.particles = particles;
}


// Used to iterate over the particles
// note: This can only be used for reading operations
// not for writing operations
// due to the use of const_iterator
#define IterateOnce(particleContainer, iteratorName, code) {
	for(std::vector<Particle>::const_iterator iteratorWithDeliberatelyLongAndComplexName1 = v.begin();
		iteratorWithDeliberatelyLongAndComplexName1 != (particleContainer).end(); ++iteratorWithDeliberatelyLongAndComplexName1)
	{
		(iteratorName) = *iteratorWithDeliberatelyLongAndComplexName1;
		(code)
	}
}

// Used to iterate over the particles in pairs
// note: This can only be used for reading operations
// not for writing operations
// due to the use of const_iterator
#define IteratePairwise(particleContainer, outerIteratorName, innerIteratorName, outerCode1, innerCode, outerCpde2) {
	for(std::vector<Particle>::const_iterator iteratorWithDeliberatelyLongAndComplexName1 = v.begin();
		iteratorWithDeliberatelyLongAndComplexName1 != (particleContainer).end(); ++iteratorWithDeliberatelyLongAndComplexName1)
	{
		(outerIteratorName) = *iteratorWithDeliberatelyLongAndComplexName1;
		(outerCode1)
		for(std::vector<Particle>::const_iterator iteratorWithDeliberatelyLongAndComplexName2 = v.begin();
		iteratorWithDeliberatelyLongAndComplexName2 != (particleContainer).end(); ++iteratorWithDeliberatelyLongAndComplexName2)
			if (iteratorWithDeliberatelyLongAndComplexName1 != iteratorWithDeliberatelyLongAndComplexName2)
			{
				(innerIteratorName) = *iteratorWithDeliberatelyLongAndComplexName2;
				(innerCode)
			}
		(outerCode2)
	}
}





