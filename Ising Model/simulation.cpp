#include <iostream>
#include <math.h>
#include "include/simulation.h"


using namespace std;

//TODO rule of three

simulation::simulation(const simulation& other){
	time = other.time;
	r = rng(other.r);
	old = lattice(other.old);
	neu = lattice(other.neu);
	temperature = other.temperature;
	print_interval = other.print_interval;
}

simulation& simulation::operator= (const simulation& other){
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
	old = lattice(size, J, H);
	neu = lattice(size, J, H);
	temperature = temp;
	time=0;
}

void simulation::advance(int time_steps, FILE* output)
{
		for (int i = 0; i< time_steps; i++)
		{
				advance();
				print_status(output);
		}
}

void simulation::advance()
{
	for (int i=0; i<old.get_size(); i++)
	{
			for (int j=0; j<old.get_size(); j++)
			{
				visit(i, j);
			}
	}
	old = lattice(neu);
	time++;
}

void simulation::visit(int row, int col)
{

	double dE = -2* old.compute_point_energy(row, col);
	double p = r.random_uniform();
	if (dE < 0 || exp (-dE/temperature) > p)
	{
		neu.flip(row, col);
	}
}

void simulation::print_status(FILE* fp){
	if (time % print_interval ==0){
		compute_means();
		fprintf(fp, "%4d\t%.9f\t%.9f\n", time, mean_magnetisation, mean_energy);
	}
}


void simulation::compute_means(){
	double total_magnetisation = 0;
	double total_energy = 0;
	int max = old.get_size();
	for (int i=0; i<max; i++){
		for (int j=0; j<max; j++){
			total_magnetisation+=old.get(i, j);
			total_energy+=old.compute_point_energy(i,j);
		}
	}
	mean_magnetisation = total_magnetisation/(max*max);
	mean_energy = total_energy/(max*max);
}
