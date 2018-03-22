#include <iostream>
#include <gsl/gsl_rng.h>
#include <cmath>
#include "rng.h"
#include "thermodynamics.h"

using namespace std;

thermodynamics::thermodynamics(double temp) noexcept
{
	temperature = temp;
}

thermodynamics::thermodynamics() noexcept: thermodynamics::thermodynamics(273.73) {}

void thermodynamics::set_temp(double temp)
{
		temperature = temp;
}

thermodynamics& thermodynamics::operator= (const thermodynamics& old)
{
	cerr<<"assignment of thermo"<<endl;
	r = old.r;
	temperature = old.temperature;
	return* this;
}

thermodynamics::~thermodynamics() noexcept
{
		cerr<<"De allocating thermodynamics"<<endl;
		return;
}

bool thermodynamics::flip_q(double dE)
{
	double p = r.random_uniform();
	return exp (-dE/temperature) > p;
}

double thermodynamics::get_temp()
{
	return temperature;
}

// int main ()
// {
// 		rng r = rng();
// 		double acc=0;
// 		for (int i=0; i<500; i++){
// 				acc +=r.random_uniform();
// 		}
// 		printf("%f", acc/500);
// 		return 0;
// }
