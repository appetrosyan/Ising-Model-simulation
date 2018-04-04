CXXFLAGS := -std=c++11 -Wall -Wextra -fopenmp -O3
CXX := /usr/local/bin/g++-7
INC := include


main: interface.o lattice.o simulation.o data.d
	${CXX}  simulation.o lattice.o interface.o -o main -lgsl -L/usr/local/opt/llvm/lib -fopenmp

simulation.o: simulation.cpp ${INC}/rng.h ${INC}/lattice.h
	${CXX} -c ${CXXFLAGS} -o simulation.o simulation.cpp

interface.o: interface.cpp ${INC}/simulation.h
	${CXX} -c ${CXXFLAGS} -o interface.o interface.cpp

lattice.o: lattice.cpp
	${CXX} -c ${CXXFLAGS} -o lattice.o lattice.cpp

.PHONY: clean data.d
clean:
	rm -f *.o
	rm -f ${INC}/*.gch
	rm -f main
	rm -f data/*

data.d:
	mkdir -p data/
