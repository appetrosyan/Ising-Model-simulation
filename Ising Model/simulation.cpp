#include <iostream>
#include <math.h>
#include "include/simulation.h"


using namespace std;

simulation::simulation(const simulation& other)
{
	time = other.time;
	r = rng(other.r);
	spins = lattice(other.spins);
	temperature = other.temperature;
	print_interval = other.print_interval;
}

simulation& simulation::operator= (const simulation& other)
{
	time= other.time;
	temperature = other.temperature;
	print_interval = other.print_interval;
	delete& r;
	r = other.r;
	delete & spins;
	spins = other.spins;
	return * this;
}

simulation::simulation(int new_size, double new_temp, double new_J, double new_H):
	time(0),
	spins(lattice(new_size, new_J, new_H)),
	temperature(new_temp),
	mean_energy(new_J*(-4)),
	total_magnetisation(new_size*new_size),
	total_energy(compute_energy(spins))
{}

void simulation::advance(int time_steps, FILE* output)
{
	int area = spins.get_size() * spins.get_size();
	for (int i = 0; i< time_steps; i++)
	{
		if (time % print_interval ==0)
		{
			print_status(output);
			total_magnetisation = spins.total_magnetisation();
			mean_magnetisation = total_magnetisation/area;
			total_energy = compute_energy(spins);
			mean_energy = total_energy/area;
		}
		advance();
	}
}

void simulation::advance()
{
	#pragma omp target teams distribute parallel for collapse(2)
	for (int row=0; row<spins.get_size(); row++)
	{
		for (int col=spins.get_size(); col>=0; col--)
		{
			double dE = compute_dE(row, col);
			double p = r.random_uniform();
			if (exp (-dE/temperature) > p)
			{
				// total_energy+=dE; //Doesn't work on multithreaded workloads
				// total_magnetisation-=2*spins.get(row, col);
				spins.flip(row, col);
			}
		}
	}
	time++;
}


double simulation::compute_dE(int row, int col){
	double e_0 = 2*spins.compute_point_energy(row, col);
	spins.flip(row, col);
	double e_1 = 2*spins.compute_point_energy(row, col);
	spins.flip(row, col);

	// Naive implementation
	// double spins_energy = compute_energy(spins);
	// spins.flip(row, col);
	// double new_energy = compute_energy(spins);
	// spins.flip(row, col);


	// //test
	// if((new_energy - spins_energy) !=(e_1 - e_0)){
	// 	cout<<(e_1 - e_0)<<" "<<(new_energy - spins_energy)<<endl;
	// }
	return e_1 - e_0;
}

void inline simulation::print_status(FILE* f)
{
	fprintf(f, "%4d\t%.9f\t%.9f\t%.9f\n", \
			time, mean_magnetisation, mean_energy, temperature);
}

double simulation::compute_energy(lattice& other)
{
	double energy=0;
	int max = other.get_size();
	#pragma omp target teams distribute parallel for collapse(2) reduction(+:energy)
	for (int i=0; i < max; i++)
	{
		for (int j=0; j < max; j++)
		{
			energy+=other.compute_point_energy(i, j);
		}
	}
	return energy;
}
