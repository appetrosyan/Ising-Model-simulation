#ifndef RNG_H
#define RNG_H

class rng
{
		const gsl_rng_type * T;
		gsl_rng * r;
public:
		rng() noexcept;
		rng(const rng& other) noexcept;
		rng& operator=(const rng& other);
		~rng() noexcept;
		double random_uniform();
};

#endif