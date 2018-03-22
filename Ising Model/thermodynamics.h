//
//  thermodynamics.hpp
//  Ising Model
//
//  Created by Aleksandr Petrosyan on 21/03/2018.
//  Copyright © 2018 Aleksandr Petrosyan. All rights reserved.
//

#ifndef thermodynamics_h
#define thermodynamics_h

#include <gsl/gsl_rng.h>
#include "rng.h"

class thermodynamics{
private:
		rng r;
    double temperature;
public:
		thermodynamics(double temperature) noexcept;
		thermodynamics(const thermodynamics& old) noexcept;
		thermodynamics() noexcept;
		~thermodynamics() noexcept;
		thermodynamics& operator=(const thermodynamics& old);
    bool flip_q(double dE);
    double get_temp();
    void set_temp(double );
};


#endif /* thermodynamics_h */
