CXXFLAGS := -std=c++11 -Wall -Wextra -fopenmp
CXXCOMPILER := /usr/local/bin/g++-7
INC := include

main: interface.o lattice.o simulation.o rng.o
	${CXXCOMPILER}  simulation.o lattice.o rng.o interface.o -o main -lgsl -fopenmp -lc++

simulation.o: simulation.cpp ${INC}/rng.h ${INC}/lattice.h
	${CXXCOMPILER} -c ${CXXFLAGS} -o simulation.o simulation.cpp

interface.o: interface.cpp ${INC}/simulation.h
	${CXXCOMPILER} -c ${CXXFLAGS} -o interface.o interface.cpp

lattice.o: lattice.cpp 
	${CXXCOMPILER} -c ${CXXFLAGS} -o lattice.o lattice.cpp

rng.o: rng.cpp
	${CXXCOMPILER} -c ${CXXFLAGS} -o rng.o rng.cpp


clean:
	rm -rf *.o
