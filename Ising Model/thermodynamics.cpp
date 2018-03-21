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
#include "rng.h"
#include "thermodynamics.h"

thermodynamics::thermodynamics(double temperature) noexcept{
	this->temperature = temperature;
}

thermodynamics::thermodynamics() noexcept: thermodynamics::thermodynamics(273.73) {}

void thermodynamics::set_temp(double temp){
		this-> temperature = temp;
}

thermodynamics::~thermodynamics() noexcept{
		delete r;
}

bool thermodynamics::flip_q(double dE){
		double p = r->random_uniform();
		return exp (dE/temperature) > p;
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
