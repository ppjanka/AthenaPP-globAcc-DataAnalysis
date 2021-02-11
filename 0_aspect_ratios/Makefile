CXX = mpicxx

COMPILE_FLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -O3 -std=c++11
LINK_FLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib -lhdf5

hdf5_reader.o: athena_hdf5_reader/athena.hpp athena_hdf5_reader/defs.hpp athena_hdf5_reader/hdf5_reader.hpp athena_hdf5_reader/athena_arrays.hpp athena_hdf5_reader/hdf5_reader.cpp
	${CXX} ${COMPILE_FLAGS} -c athena_hdf5_reader/hdf5_reader.cpp

aspect_ratios.o: aspect_ratios.cpp aspect_ratios.hpp
	${CXX} ${COMPILE_FLAGS} -c aspect_ratios.cpp

aspect_ratios: aspect_ratios.o hdf5_reader.o
	${CXX} -o aspect_ratios aspect_ratios.o hdf5_reader.o ${LINK_FLAGS}

PHONY: clean
clean:
	rm *.o *.pyc

PHONY: cleanall
cleanall:
	rm *.o *.pyc aspect_ratios