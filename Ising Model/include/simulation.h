#ifndef simulation_h
#define simulation_h

#include "lattice.h"
#include "rng.h"
#include <gsl/gsl_rng.h>

class simulation
{
private:
		int time=0;
		rng r = rng();
		lattice spins;
		double temperature;
		double mean_magnetisation=1;
		double mean_energy;
		double total_magnetisation;
		double total_energy;
public:

		int print_interval=1;
		simulation(int size, double temp, double J, double H);
		simulation(const simulation& other);
		simulation& operator= (const simulation& other);
		~simulation(){}
		void advance();
		// void inline visit(int, int);
		void inline print_status(FILE* fp);
		void advance(int time_steps, FILE* output);
		double compute_energy(lattice & other);
		double compute_dE(int row, int col);
		// double point_energy_with_neightbours(int row, int col);
};

#endif /* simulation_h */
