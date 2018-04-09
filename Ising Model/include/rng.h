#ifndef RNG_H
#define RNG_H

#include <gsl/gsl_rng.h>
/*
This is a helper class to encapsulate the state of GNU Scientific Library's
random number generator.

While GSL_rng could be used directly, the end-user (i.e. you), would have to
keep track of memory, free it as needed, and potentially cause a memory leak.

It's a minor annoyance, but in order to make sure that the RNG's state is being
adequately preserved across multiple simulation stages.
*/
class rng {
  const gsl_rng_type *T;
  gsl_rng *r;
  unsigned long int Seed = 23410981;
public:
  rng() noexcept {
    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, Seed);
  }
  // "Rule of three": given we need any of the custom destructor, copy
  // constructor or copy-assignment operator, we need all three.
  rng(const rng &other) noexcept { r = other.r; }
  rng &operator=(const rng &other) {
    delete &r;
    r = other.r;
    return *this;
  }
  ~rng() noexcept { gsl_rng_free(r); }
  // Generates random number between 0 and 1.
  double random_uniform() { return gsl_rng_uniform(r); }
  int random_int(int max){ return gsl_rng_uniform_int(r, max);}
};

#endif
