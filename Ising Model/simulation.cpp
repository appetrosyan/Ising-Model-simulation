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
	for (int i = 0; i< time_steps; i++)
	{
		if (time % print_interval ==0)
		{
			print_status(output);
		}
		advance();
	}
}

void simulation::advance()
{
	// cout<<"advancing--------"<<endl;
	// #pragma omp parallel for ordered reduction(+:total_energy) reduction(+:total_magnetisation)
	for (int row=0; row<spins.get_size(); row++)
	{
		for (int col=spins.get_size(); col>=0; col--)
		{
			double dE = compute_dE(row, col);
			double p = r.random_uniform();
			if (exp (-dE/temperature) > p)
			{
				// #pragma omp ordered
				total_energy+=dE;
				// cout<<total_energy<<endl;
				// cout<<"flip "<<row<<" "<<col<<":\t"<<dE<<"\t"<<exp (-dE/temperature)<<"\t"<<p<<endl;
				total_magnetisation-=2*spins.get(row, col);
				spins.flip(row, col);
			}
		}
	}
	int area = spins.get_size()*spins.get_size();
	mean_magnetisation = total_magnetisation/area;
	// //Testing code
	// double energy_estimate = compute_energy(spins);
	// double magnetisation_estimate = spins.total_magnetisation();
	// if(total_energy != energy_estimate){
	// 	cout<<total_energy<<"! = "<<energy_estimate<<":: "<<time<<endl;
	// 	total_energy = energy_estimate;
	// 	// spins.print();
	// }
	// if(total_magnetisation != magnetisation_estimate){
	// 	cout<<total_magnetisation<<"! = "<<magnetisation_estimate<<": "<<time<<endl;
	// 	total_magnetisation = magnetisation_estimate;
	// }
	// // end testing code
	mean_energy= total_energy/area;
	time++;
}


double simulation::compute_dE(int row, int col){
	//Real code
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
	return e_1-e_0;
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

	#pragma omp parallel for collapse(2) reduction(+:energy)
	for (int i=0; i<max; i++)
	{
		for (int j=0; j<max; j++)
		{
			energy+=other.compute_point_energy(i, j);
		}
	}
	return energy;
}
