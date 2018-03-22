CXXFLAGS := -std=c++11 -Wall -Wextra

main: simulation.o lattice.o
		g++  simulation.o lattice.o rng.o -o main -lgsl

simulation.o: simulation.cpp rng.o 
		g++ -c ${CXXFLAGS} -o simulation.o simulation.cpp

lattice.o: lattice.cpp include/lattice.h
		g++ -c ${CXXFLAGS} -o lattice.o lattice.cpp

rng.o: rng.cpp include/rng.h
	g++ -c ${CXXFLAGS} -o rng.o rng.cpp

clean:
	rm -rf *.o && rm -rf *.out rm -rf *.gch && rm -rf main
