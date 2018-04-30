#include "include/simulation.h"
#include <iostream>
#include <math.h>


using namespace std;

/*
Advances the simulation a given number of steps, and updates/prints the statistics
into the given file pointer.

Defaults to stdout.

The number of time_steps is explcitly unsigned, so that linters/IDEs remind
the end user of the file that extra care needs to be taken, as well as to allow
advancing the simulation a larger number of times.
*/
void simulation::advance(unsigned int time_steps, FILE *output) {
  unsigned int area = spin_lattice_.get_size() * spin_lattice_.get_size();
    for (unsigned int i = 0; i < time_steps; i++) {
      if (time_ % print_interval_ == 0) {
        mean_magnetisation_ = total_magnetisation_ / area;
        mean_energy_ = total_energy_ / area;
        mean_energy_squared_ =      compute_energy_squared(spin_lattice_)/area;
        print_status(output);
      }
      advance();
    }
}

/*
Advances the simulation a single step.

DOES NOT KEEP TRACK OF STATISTICS. Hence private.
*/
void simulation::advance() {
  unsigned int size = spin_lattice_.get_size();
  double dE, p;
  // # pragma omp parallel for collapse(2)
  // This slows down the python multiprocessing
  for (unsigned int row=0; row<size; row++){
    for (unsigned int col = 0; col < size; col++) {
      dE = compute_dE(row, col);
      p = r_.random_uniform();
      if (dE <0 || exp(-dE * inverse_temperature_) > p) {
        spin_lattice_.flip(row, col);
        total_magnetisation_ += 2* spin_lattice_.get(row, col);
        total_energy_ += dE;
      }
    }
  }
  time_++;
}



/*
Computes the total energy associated with spins in the spin_lattice_.
*/
double simulation::compute_energy_squared(lattice &other) {
  double energy_squared_sum = 0, point_energy;
  unsigned int size = other.get_size();
  // #pragma omp parallel for collapse(2)
  for (unsigned int i = 0; i < size; i++) {
    for (unsigned int j = 0; j < size; j++) {
      point_energy = other.compute_point_energy(i, j);
      energy_squared_sum += point_energy*point_energy;
    }
  }
  return energy_squared_sum/4;
}


void simulation::set_to_chequerboard(int step){
  if (time_ !=0){
    return;
  }else{
    unsigned int max = spin_lattice_.get_size();
    # pragma omp parallel for collapse(2)
    for (unsigned int i=0; i< max; ++i){
      for (unsigned int j=0; j<max; ++j){
        if ((i/step)%2-(j/step)%2==0){
          spin_lattice_.flip(i, j);
        }
      }
    }
    total_magnetisation_ = spin_lattice_.total_magnetisation();
  }
}
