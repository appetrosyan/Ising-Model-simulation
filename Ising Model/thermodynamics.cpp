//
//  thermodynamics.cpp
//  Ising Model
//
//  Created by Aleksandr Petrosyan on 21/03/2018.
//  Copyright © 2018 Aleksandr Petrosyan. All rights reserved.
//

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <cmath>

class rng{
		const gsl_rng_type * T;
		gsl_rng * r;
public:
		rng(){
				gsl_rng_env_setup();
				T = gsl_rng_default;
				r = gsl_rng_alloc (T);
		}

		~rng(){
				gsl_rng_free (r);
		}

		double random_uniform(){
				return gsl_rng_uniform(r);
		}
};


class thermod{
private:
		rng* r;
		double temperature;
public:
		thermod();
		thermod(double temperature){

		}

		void set_temperature(double temp){
				this-> temperature = temp;
		}

		bool flip_q(double dE){
				double p = r->random_uniform();
				return exp (dE/temperature) > p;
		}

		~thermod(){
				delete r;
		}
};


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
