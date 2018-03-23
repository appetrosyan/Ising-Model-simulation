#include <iostream>
#include <math.h>
#include "include/simulation.h"


using namespace std;

simulation::simulation(const simulation& other)
{
	time = other.time;
	r = rng(other.r);
	old = lattice(other.old);
	neu = lattice(other.neu);
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
	delete & old;
	old = other.old;
	delete & neu;
	neu = other.neu;
	return * this;
}

simulation::simulation(int size, double temp, double J, double H)
{
	this->old = lattice(size, J, H);
	this->neu = lattice(size, J, H);
	this->temperature = temp;
	this->total_magnetisation = size*size;
	this->time=0;
}

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

	#pragma omp parallel for collapse(2) reduction(+:total_magnetisation)
	for (int row=0; row<old.get_size(); row++)
	{
		for (int col=0; col<old.get_size(); col++)
		{
			double dE = compute_dE(row, col);
			double p = r.random_uniform();
			if (dE < 0 || exp (-dE/temperature) > p)
			{
				total_magnetisation-=2*neu.get(row, col);
				neu.flip(row, col);
			}
		}
	}
	int area = old.get_size()*old.get_size();
	mean_magnetisation = total_magnetisation/area;
	// mean_energy = compute_energy(neu)/area;
	old = lattice(neu);
	time++;
}

double simulation::compute_dE(int row, int col){
	//TODO implement
	return row*col;
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
