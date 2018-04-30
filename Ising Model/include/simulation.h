#ifndef simulation_h
#define simulation_h

#include "lattice.h"
#include "rng.h"
#include <gsl/gsl_rng.h>

/*
The logic of the entire simulation of the Ising model of magnetism.

This simulation will run and print statistics at a given time interval.
A simulation can be advanced a single time step, or many at a time,
*/
class simulation {
private:
  unsigned int time_ = 0;  // Current time of the simulation.
  rng r_ = rng();
  lattice spin_lattice_;
  double inverse_temperature_;
  double mean_magnetisation_ = 1;
  double mean_energy_;
  double mean_energy_squared_;
  double total_magnetisation_;
  double total_energy_;
  unsigned int print_interval_ = 1;
  void advance();

public:
  void set_print_interval(unsigned int npi) { print_interval_ = npi; }
  simulation(int new_size, double new_temp, double new_J, double new_H)
      : time_(0),
        spin_lattice_(lattice(new_size, new_J, new_H)),
        inverse_temperature_(1/new_temp),
        mean_energy_(new_J * (-2)),
        total_magnetisation_(new_size * new_size),
        total_energy_(-2*new_size*new_size) {}

  void print_status(FILE *f) {
    f = f==NULL? stdout : f;
    double variance = -mean_energy_*mean_energy_ + mean_energy_squared_;
    fprintf(f, "%4d\t%e \t%e\t%e\n", time_,
            mean_magnetisation_, mean_energy_, variance);
  }
  void advance(unsigned int time_steps, FILE *output);
  double compute_energy_squared(lattice &other);

  /*
  Computes change in energy due to flipping one single spin.
  */
  double compute_dE(int row, int col) {
    return -2*spin_lattice_.compute_point_energy(row, col);
  }
  void set_to_chequerboard(int step);
};

#endif
