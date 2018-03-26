#ifndef simulation_h
#define simulation_h

#include "lattice.h"
#include "rng.h"
#include <gsl/gsl_rng.h>

class simulation {
private:
  int time_ = 0;
  rng r_ = rng();
  lattice spin_lattice_;
  double temperature_;
  double mean_magnetisation_ = 1;
  double mean_energy_;
  double total_magnetisation_;
  double total_energy_;

public:
  int print_interval_ = 1;
  simulation(int new_size, double new_temp, double new_J, double new_H)
      : time_(0), spin_lattice_(lattice(new_size, new_J, new_H)), temperature_(new_temp),
        mean_energy_(new_J * (-4)), total_magnetisation_(new_size * new_size),
        total_energy_(compute_energy(spin_lattice_)) {}
  void advance();
  void inline print_status(FILE *f) {
    fprintf(f, "%4d\t%.9f \t%.9f\t%.9f\n", time_, mean_magnetisation_,
            mean_energy_, temperature_);
  }
  void advance(int time_steps, FILE *output);
  double compute_energy(lattice &other);
  double compute_dE(int row, int col);
};

#endif
