CXXFLAGS := -std=c++11 -Wall -Wextra
CXXCOMPILER := g++
BINARIES := bin

main: simulation.o lattice.o
	${CXXCOMPILER}  ${BINARIES}/simulation.o ${BINARIES}/lattice.o ${BINARIES}/rng.o -o ${BINARIES}/main -lgsl


simulation.o: simulation.cpp rng.o
	${CXXCOMPILER} -c ${CXXFLAGS} -o ${BINARIES}/simulation.o simulation.cpp

lattice.o: lattice.cpp include/lattice.h
	${CXXCOMPILER} -c ${CXXFLAGS} -o ${BINARIES}/lattice.o lattice.cpp

rng.o: rng.cpp include/rng.h
	${CXXCOMPILER} -c ${CXXFLAGS} -o ${BINARIES}/rng.o rng.cpp


clean:
	rm -rf ${BINARIES}/* && rm -rf main
