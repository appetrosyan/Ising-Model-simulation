#include <gsl/gsl_rng.h>
#include <iostream>
#include "include/rng.h"

using namespace std;

rng::rng() noexcept
{
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
}

rng::rng(const rng& other) noexcept{
	r = other.r;
}

rng& rng::operator=(const rng& other){
	// delete& r;
	r = other.r;
	return * this;
}

rng::~rng() noexcept
{
	gsl_rng_free (r);
}

double rng::random_uniform()
{
		return gsl_rng_uniform(r);
}
