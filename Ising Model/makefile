CXXFLAGS := -std=c++11 -Wall -Wextra

main: simulation.o thermodynamics.o lattice.o
		g++  simulation.o lattice.o  thermodynamics.o  rng.o -o main -lgsl

simulation.o: simulation.cpp lattice.h  thermodynamics.h
		g++ -c ${CXXFLAGS} -o simulation.o simulation.cpp

thermodynamics.o: thermodynamics.cpp thermodynamics.h rng.o
		g++ -c ${CXXFLAGS} -o thermodynamics.o thermodynamics.cpp

lattice.o: lattice.cpp lattice.h
		g++ -c ${CXXFLAGS} -o lattice.o lattice.cpp

rng.o: rng.cpp rng.h
	g++ -c ${CXXFLAGS} -o rng.o rng.cpp

clean:
	rm -rf *.o
