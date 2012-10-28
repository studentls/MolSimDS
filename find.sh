find src -name *.cpp -printf '%p\\\n' > files.mk.tmp
echo "SOURCES=\\" > files.mk
cat files.mk.tmp >> files.mk
rm files.mk.tmp


