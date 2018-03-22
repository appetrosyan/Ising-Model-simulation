CXXFLAGS := -std=c++11 -Wall -Wextra
CXXCOMPILER := g++

main: interface.o lattice.o simulation.o rng.o
	${CXXCOMPILER}  simulation.o lattice.o rng.o interface.o -o main -lgsl

simulation.o: simulation.cpp rng.o
	${CXXCOMPILER} -c ${CXXFLAGS} -o simulation.o simulation.cpp

interface.o: interface.cpp simulation.o
	${CXXCOMPILER} -c ${CXXFLAGS} -o interface.o interface.cpp

$BINARIES}/lattice.o: lattice.cpp
	${CXXCOMPILER} -c ${CXXFLAGS} -o lattice.o lattice.cpp

rng.o: rng.cpp
	${CXXCOMPILER} -c ${CXXFLAGS} -o rng.o rng.cpp


clean:
	rm -rf *.o 
