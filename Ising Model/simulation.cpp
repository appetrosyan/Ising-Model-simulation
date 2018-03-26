#include <iostream>
#include <math.h>
#include "include/simulation.h"

using namespace std;

void simulation::advance(int time_steps, FILE *output) {
  int area = spin_lattice_.get_size() * spin_lattice_.get_size();
  for (int i = 0; i < time_steps; i++) {
    if (time_ % print_interval_ == 0) {
      print_status(output);
      total_magnetisation_ = spin_lattice_.total_magnetisation();
      mean_magnetisation_ = total_magnetisation_ / area;
      total_energy_ = compute_energy(spin_lattice_);
      mean_energy_ = total_energy_ / area;
    }
    advance();
  }
}

void simulation::advance() {
  #pragma omp parallel for collapse(2)
  for (int row = 0; row < spin_lattice_.get_size(); row++) {
    for (int col = spin_lattice_.get_size(); col >= 0; col--) {
      double dE = compute_dE(row, col);
      double p = r_.random_uniform();
      if (exp(-dE / temperature_) > p) {
        // total_energy+=dE;	//Don't work when multithreaded
        // total_magnetisation-=2*spin_lattice_.get(row, col);
        spin_lattice_.flip(row, col);
      }
    }
  }
  time_++;
}

double simulation::compute_dE(int row, int col) {
  double e_0 = 2 * spin_lattice_.compute_point_energy(row, col);
  spin_lattice_.flip(row, col);
  double e_1 = 2 * spin_lattice_.compute_point_energy(row, col);
  spin_lattice_.flip(row, col);

  // Naive implementation
  // double spin_lattice__energy = compute_energy(spin_lattice_);
  // spin_lattice_.flip(row, col);
  // double new_energy = compute_energy(spin_lattice_);
  // spin_lattice_.flip(row, col);

  // //test
  // if((new_energy - spin_lattice__energy) !=(e_1 - e_0)){
  // 	cout<<(e_1 - e_0)<<" "<<(new_energy - spin_lattice__energy)<<endl;
  // }
  return e_1 - e_0;
}

double simulation::compute_energy(lattice &other) {
  double energy_sum = 0;
  int max = other.get_size();
  #pragma omp parallel for reduction(+ : energy_sum)
  for (int i = 0; i < max; i++) {
    for (int j = 0; j < max; j++) {
      energy_sum += other.compute_point_energy(i, j);
    }
  }
  return energy_sum;
}
