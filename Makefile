all:
	g++ -fPIC -g -c -Wall psort.cpp -fopenmp
	g++ -shared -Wl,-soname,libsort.so -o libsort.so psort.o -lc
	g++ -fopenmp -Wall psort.cpp driver_binary.cpp -o sort
clean:
	rm *.so *.o sort