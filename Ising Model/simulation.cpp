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
      print_status(output);
      total_magnetisation_ = spin_lattice_.total_magnetisation();
      mean_magnetisation_ = total_magnetisation_ / area;
    }
    // If we don't update mean_energy_ every time, we might get incorrect
    // thermodynamic behaviour.
    total_energy_ = compute_energy(spin_lattice_);
    long double temperature_delta = total_energy_/area - mean_energy_;
    if (abs(temperature_delta) < 1/area){
      cerr<<temperature_delta<<"Reached equilibrium"<<endl;
    }
    temperature_ += temperature_delta;
    mean_energy_ = total_energy_ / area;
    advance();
  }
}

/*
Advances the simulation a single step.

DOES NOT KEEP TRACK OF STATISTICS. Hence private.
*/
void simulation::advance() {
  #pragma omp parallel for collapse(2)
  for (unsigned int row = 0; row < spin_lattice_.get_size(); row++) {
    for (unsigned int col = 0; col < spin_lattice_.get_size(); col++) {
      double dE = compute_dE(row, col);
      double p = r_.random_uniform();
      if (exp(-dE / temperature_) > p) {
        // Not thread safe. see comment in compute_energy(lattice& )
        // total_energy_+=dE;
        // total_magnetisation_-=2*spin_lattice_.get(row, col);
        spin_lattice_.flip(row, col);
      }
    }
  }
  time_++;
}

/*
Computes change in energy due to flipping one single spin.

The function returns a single-precision floating-point number, as data cannot under
most circumstances make use of greater precision than that (save J is set to a
non-machine-representable value).

The code modifies the spin lattice, as an alternative (copying the neighborhood
of a given point), would make the code run slower by a factor of 2.25
*/
float simulation::compute_dE(int row, int col) {
  float e_0 = 2 * spin_lattice_.compute_point_energy(row, col);
  spin_lattice_.flip(row, col);
  float e_1 = 2 * spin_lattice_.compute_point_energy(row, col);
  spin_lattice_.flip(row, col);
  /*
  // Naive implementation
  double spin_lattice__energy = compute_energy(spin_lattice_);
  spin_lattice_.flip(row, col);
  double new_energy = compute_energy(spin_lattice_);
  spin_lattice_.flip(row, col);

  // test
  if((new_energy - spin_lattice__energy) !=(e_1 - e_0)){
  	cout<<(e_1 - e_0)<<" "<<(new_energy - spin_lattice__energy)<<endl;
  }
  */
  return e_1 - e_0;
}
/*
Computes the total energy associated with spins in the spin_lattice_.

I originally used this function to test the code that tracked energy as the lattice
itself was modified, but that code turned out to be only marginally faster, and
not thread-safe. This is due to a race condition: when one thread uses a neighborhood
of a point, while another thread was computing the energy of one such point in
the neighborhood of (row, col).
*/
long double simulation::compute_energy(lattice &other) {
  double energy_sum = 0;
  unsigned int max = other.get_size();
  #pragma omp parallel for reduction(+ : energy_sum)
  for (unsigned int i = 0; i < max; i++) {
    for (unsigned int j = 0; j < max; j++) {
      energy_sum += other.compute_point_energy(i, j);
    }
  }
  return energy_sum;
}
