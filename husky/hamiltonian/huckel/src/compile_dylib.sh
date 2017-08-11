# compile the huckel program to a dynamic library
# that can be called within python

gcc -dynamiclib huckel.c -o hkl.dylib
mv hkl.dylib ../
