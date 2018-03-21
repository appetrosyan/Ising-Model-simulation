//
//  simulation.h
//  Ising Model
//
//  Created by Aleksandr Petrosyan on 21/03/2018.
//  Copyright © 2018 Aleksandr Petrosyan. All rights reserved.
//

#ifndef simulation_h
#define simulation_h

class simulation {
private:
    int time;
    double temperature;
    lattice l;
public:
    void advance(int time_steps);
};

#endif /* simulation_h */
