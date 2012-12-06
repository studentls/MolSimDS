xsd cxx-tree --output-dir ./src/XMLFile --hxx-suffix .h --cxx-suffix .cpp simulationfile.xsd
xsd cxx-tree --generate-serialization --generate-ostream --output-dir ./src/outputWriter/ --hxx-suffix .h --cxx-suffix .cpp ./src/outputWriter/vtk-unstructured.xsd
