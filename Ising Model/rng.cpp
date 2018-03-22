#include <gsl/gsl_rng.h>
#include "rng.h"

rng::rng() noexcept
{
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
}

rng::~rng() noexcept
{
		gsl_rng_free (r);
}

double rng::random_uniform()
{
		return gsl_rng_uniform(r);
}
