//
//  simluation.cpp
//  Ising Model
//
//  Created by Aleksandr Petrosyan on 21/03/2018.
//  Copyright © 2018 Aleksandr Petrosyan. All rights reserved.
//

#include <iostream>
#include "lattice.h"
#include "thermodynamics.h"

using namespace std;

class simulation {
private:
		int time;
		lattice* l;
		thermodynamics* t_d;
public:

		simulation(int size, double temperature, double J, double H){
			this->t_d = new thermodynamics(temperature);
			this->l = new lattice(size, J, H);
			this->time=0;
		}

		void advance(){
				lattice out = lattice(*l);
				for (int i=0; i<l->get_size(); i++){
						for (int j=0; j<l->get_size(); j++){
								visit(i, j, out);
						}
				}
				time++;
		}

		void advance(int time_steps){
				for (int i = 0; i< time_steps; i++){
						advance();
				}
		}

		void visit(int row, int col, lattice out){

				double dE = -2* l->compute_point_energy(row, col);
				if (dE <0 || t_d->flip_q(dE)) out.flip(row, col);
		}


};

int main(){
	simulation* s = new simulation(4, 273.73, 1.0, 0.0);
	s->advance();
	return 0;
}
