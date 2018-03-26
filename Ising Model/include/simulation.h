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
  long double temperature_;
  long double mean_magnetisation_ = 1;
  long double mean_energy_;
  long double total_magnetisation_;
  long double total_energy_;
  unsigned int print_interval_ = 1;
  void advance();

public:
  void set_print_interval(unsigned int new_print_interval) { print_interval_ = new_print_interval; }

  simulation(int new_size, double new_temp, double new_J, double new_H)
      : time_(0), spin_lattice_(lattice(new_size, new_J, new_H)), temperature_(new_temp),
        mean_energy_(new_J * (-4)), total_magnetisation_(new_size * new_size),
        total_energy_(compute_energy(spin_lattice_)) {}

  void print_status(FILE *f) {
    f = f==NULL? stdout : f;
    fprintf(f, "%4d\t%Le \t%Le\t%Le\n", time_, mean_magnetisation_,
            mean_energy_, temperature_);
  }
  void advance(unsigned int time_steps, FILE *output);
  long double compute_energy(lattice &other);
  float compute_dE(int row, int col);
};

#endif
