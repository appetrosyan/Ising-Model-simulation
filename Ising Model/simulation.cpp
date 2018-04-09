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
    /*
    To simulate an adiabatic Ising lattice, uncomment this
    */
    // double temperature_delta = total_energy_/area - mean_energy_;
    // if (abs(temperature_delta) < 1/area){
    //   cerr<<temperature_delta<<"! Reached equilibrium "<<endl;
    // }
    // temperature_ += temperature_delta;

    if (time_ % print_interval_ == 0) {
      total_magnetisation_ = spin_lattice_.total_magnetisation();
      mean_magnetisation_ = total_magnetisation_ / area;
      total_energy_ = compute_energy(spin_lattice_);
      mean_energy_ = total_energy_ / area;
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
  // #pragma omp parallel for collapse(2)
  for (unsigned int row = 0; row < spin_lattice_.get_size(); row++) {
    for (unsigned int col = 0; col < spin_lattice_.get_size(); col++) {
      double dE = compute_dE(row, col);
      double p = r_.random_uniform();
      // float rnd = rand() / (RAND_MAX + 1.);
      if (dE <0 || exp(-dE / temperature_) > p) {
        spin_lattice_.flip(row, col);
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
double simulation::compute_energy(lattice &other) {
  double energy_sum = 0;
  unsigned int max = other.get_size();
  #pragma omp parallel for reduction(+ : energy_sum)
  for (unsigned int i = 0; i < max; i++) {
    for (unsigned int j = 0; j < max; j++) {
      energy_sum += other.compute_point_energy(i, j);
    }
  }
  return energy_sum/2;
}


void simulation::set_to_chequerboard(int step){
  if (time_ !=0){
    return;
  }else{
    for (unsigned int i=0; i< spin_lattice_.get_size(); ++i){
      for (unsigned int j=0; j<spin_lattice_.get_size(); ++j){
        if ((i/step)%2-(j/step)%2==0){
          spin_lattice_.flip(i, j);
        }
      }
    }
  }
}


/*
Loads the lattice l into the simulation s.

*/
// void simulation::load_custom(const lattice& custom){
//   spin_lattice_ = custom;
// }
