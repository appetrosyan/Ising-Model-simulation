#ifndef RNG_H
#define RNG_H

#include <gsl/gsl_rng.h>

class rng
{
	const gsl_rng_type * T;
	gsl_rng * r;
public:
	rng() noexcept {
			gsl_rng_env_setup();
			T = gsl_rng_default;
			r = gsl_rng_alloc (T);
	}
	rng(const rng& other) noexcept { r = other.r; }
	rng& operator=(const rng& other)
	{
		delete& r;
		r = other.r;
		return * this;
	}
	~rng() noexcept { gsl_rng_free (r); }
	double random_uniform() { return gsl_rng_uniform(r); }
};

#endif
