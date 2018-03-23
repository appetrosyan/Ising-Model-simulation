//
//  simulation.h
//  Ising Model
//
//  Created by Aleksandr Petrosyan on 21/03/2018.
//  Copyright © 2018 Aleksandr Petrosyan. All rights reserved.
//

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
		double temperature;
		double mean_magnetisation;
		double mean_energy;
		double total_magnetisation = 1;
public:
		lattice neu;
		lattice old;
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
};

#endif /* simulation_h */
