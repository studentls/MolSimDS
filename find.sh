find src -name *.cpp -printf '%p\\\n' > files.mk.tmp
find src -name *.cxx -printf '%p\\\n' > files.mk.tmp2
echo "SOURCES=\\" > files.mk
cat files.mk.tmp >> files.mk
cat files.mk.tmp2 >> files.mk
rm files.mk.tmp
rm files.mk.tmp2

