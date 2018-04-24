#include <iostream>
#include <math.h>
#include "include/simulation.h"

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
      mean_energy_squared_ = compute_energy_squared(spin_lattice_)/area;
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
  for (unsigned int row=0; row<size; row++){
    for (unsigned int col = 0; col < size; col++) {
      double dE = compute_dE(row, col);
      double p = r_.random_uniform();
      if (dE <0 || exp(-dE / temperature_) > p) {
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

I originally used this function to test the code that tracked energy as the lattice
itself was modified, but that code turned out to be only marginally faster, and
not thread-safe. This is due to a race condition: when one thread uses a neighborhood
of a point, while another thread was computing the energy of one such point in
the neighborhood of (row, col).
*/
double simulation::compute_energy_squared(lattice &other) {
  double energy_squared_sum = 0, point_energy;
  unsigned int size = other.get_size();
  // #pragma omp target teams distribute parallel for reduction(+:energy_squared_sum)
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
