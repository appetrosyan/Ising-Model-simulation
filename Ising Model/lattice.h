//
//  lattice.h
//  Ising Model
//
//  Created by Aleksandr Petrosyan on 21/03/2018.
//  Copyright © 2018 Aleksandr Petrosyan. All rights reserved.
//

#ifndef lattice_h
#define lattice_h

class lattice{
    int size;
    short* spin;
    double J;
    double H;
public:
    lattice(int _size, double _J, double _H);
    void print();
    short get(int row, int col);
    void flip(int row, int col);
    int get_size();
    double compute_point_energy(int row, int col);
};
#endif /* lattice_h */
