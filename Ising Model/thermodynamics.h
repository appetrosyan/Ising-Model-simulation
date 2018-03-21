//
//  thermodynamics.hpp
//  Ising Model
//
//  Created by Aleksandr Petrosyan on 21/03/2018.
//  Copyright © 2018 Aleksandr Petrosyan. All rights reserved.
//

#ifndef thermodynamics_h
#define thermodynamics_h

class thermod{
private:
    double temperature;
public:
		thermod();
		thermod(double temperature);
    bool flip_q(double dE);
    double get_temp();
    void set_temp(double );
};


#endif /* thermodynamics_h */
