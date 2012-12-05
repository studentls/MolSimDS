//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File XMLFileReaderTest.h
// contains tests for the class XMLFileReader
//------------------------------------------------------------------------------------------------
/// @author F.Dietz
/// @author L.Spiegelberg
//------------------------------------------------------------------------------------------------

#ifndef XMLFILEREADER_TEST_HEADER_
#define XMLFILEREADER_TEST_HEADER_

#include "XMLFileReader.h"
#include "utils/utils.h"

//CPP Unit
#include <cppunit/TestCase.h>
#include <cppunit/TestFixture.h>
#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>




// see http://cppunit.sourceforge.net/doc/lastest/cppunit_cookbook.html for details

/// Particle Container test class
class XMLFileReaderTest : public CppUnit::TestFixture
{	
	//fast setup macros
	CPPUNIT_TEST_SUITE(XMLFileReaderTest);
	CPPUNIT_TEST(testXMLFileReaderStability);
	CPPUNIT_TEST(testXMLFileReaderParsing);
	CPPUNIT_TEST_SUITE_END();
private:
	void	makeXMLFile(const char *filename)
	{
		FILE *pFile = fopen(filename, "w");
		CPPUNIT_ASSERT(pFile);

		fprintf(pFile,	"<?xml version=\"1.0\"?>\n" \
						"<simulationfile xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n" \
						"	   xsi:noNamespaceSchemaLocation=\"simulationfile.xsd\">\n" \
						"	  <params>\n" \
						"	<output>out</output>\n" \
    					"		<iterationsperoutput>10</iterationsperoutput>\n" \
						"	<brownianMotionFactor> 0.1 </brownianMotionFactor> \n"   
						"	<delta_t> 0.02 </delta_t>\n" \
						"	<epsilon>5.0</epsilon>\n" \
						"	<outputfmt>VTK</outputfmt>\n" \
						"	<sigma>1.0</sigma>\n"     
						"	<t_end>5.0</t_end>\n" \
						"	<t_start>0.0</t_start>\n"     
						"	<algorithm>\n" \
						"	  <List>\n" \
						"	  </List>\n" \
						"	</algorithm>\n"     
						" </params>\n" \
						" <data>\n" \
						"	<cuboid><X>0.0 1.0 0.0</X><V>1.0 0.0 0.0</V><N>20 20 1</N><h>0.1</h><m>1.0</m></cuboid>\n" \
						" </data>\n" \
						"</simulationfile>\n");


		fclose(pFile);
	}

public:
	/// function testing the behaviour if functions are not called appropriate
	void testXMLFileReaderStability()
	{
		XMLFileReader fr;

		// what happens if makeParticleContainer is called before a file is read?
		CPPUNIT_ASSERT(FAILED(fr.makeParticleContainer(NULL)));

		// what happens if a valid Pointer is given?
		ParticleContainer *pc = new ListParticleContainer();

		CPPUNIT_ASSERT(FAILED(fr.makeParticleContainer(&pc)));

		SAFE_DELETE(pc);
	}

	/// function, which will lay down a test file in the directory
	void testXMLFileReaderParsing()
	{
		char filename[] = "test.xml";

		// file present?
		LOG4CXX_INFO(generalOutputLogger, "file test.xml will be overwritten!");

		//make file
		makeXMLFile(filename);

		//now read it
		XMLFileReader fr;

		CPPUNIT_ASSERT(SUCCEEDED(fr.readFile(filename)));

		// get Container
		ParticleContainer *pc = NULL;

		CPPUNIT_ASSERT(SUCCEEDED(fr.makeParticleContainer(&pc)));

		// correct count of particles?
		// grid shall be 20x20x1 => 400 Particles
		CPPUNIT_ASSERT(pc->getParticleCount() == 400);
		
		//test value of t_end, shall be 5.0
		CPPUNIT_ASSERT(fr.getDescription().end_time == 5.0);

		CPPUNIT_ASSERT(fr.getDescription().iterationsperoutput == 10);

		SAFE_DELETE(pc);
	}
};



#endif