#include <iostream>
#include <math.h>
#include "include/lattice.h"
#include "include/rng.h"
#include <gsl/gsl_rng.h>

using namespace std;

class simulation
{
private:
		int time;
		rng r = rng();
		lattice old;
		double temperature;
public:
		lattice neu;
		// thermodynamics t_d;
		simulation(int size, double temp, double J, double H)
		{
			// cerr<<"Entered constructor"<<endl;
			// t_d = thermodynamics(temperature);
			// cerr<<"Passed thermo assignment";
			old = lattice(size, J, H);
			// old. print();
			neu = lattice(size, J, H);
			temperature = temp;
			// neu. print();
			time=0;
		}

		~simulation(){
			cerr<<"De Allocating simulation"<<endl;
		}

		void advance()
		{
			// old.print();
			// neu.print();
			for (int i=0; i<old.get_size(); i++)
			{
					for (int j=0; j<old.get_size(); j++)
					{
						visit(i, j);
					}
			}
			old = lattice(neu);
			time++;
			// cerr<<" ---"<<time<<endl;
		}

		void advance(int time_steps)
		{
				for (int i = 0; i< time_steps; i++)
				{
						advance();
						cout<<mean_magnetisation()<<endl;
				}
		}

		void visit(int row, int col)
		{

			double dE = -2* old.compute_point_energy(row, col);
			double p = r.random_uniform();
			if (dE < 0 || exp (-dE/temperature) > p)
			{
				neu.flip(row, col);
				// cerr<<"Flip"<<endl;
			}
			// cerr<<" "<<t_d.get_temp()<<t_d.flip_q(dE);
		}

		double mean_magnetisation(){
			double accumulator=0;
			int max = old.get_size();
			for (int i=0; i<max; i++){
				for (int j=0; j<max; j++){
					accumulator+=old.get(i, j);
				}
			}
			return accumulator/(max*max);
		}


};

int main()
{
	simulation s = simulation(100, 1, 0.3, 0.0);
	// k.print();
	// s.neu.print();
	s.advance(15);
	// delete& s;
	// cerr<<endl<<endl<<"Done"<<endl;
	return 0;
}
