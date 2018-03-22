#include <iostream>
#include "lattice.h"
#include "thermodynamics.h"
#include "rng.h"

using namespace std;

class simulation
{
private:
		int time;
public:
		thermodynamics t_d;
		lattice old;
		lattice neu;
		simulation(int size, double temperature, double J, double H)
		{
			t_d = thermodynamics(temperature);
			old = lattice(size, J, H);
			old. print();
			neu = lattice(size, J, H);
			neu. print();
			time=0;
		}

		~simulation(){
			cerr<<"De Allocating simulation"<<endl;
		}

		void advance()
		{
			// cerr<<"calling advance";
			for (int i=0; i<old.get_size(); i++)
			{
					for (int j=0; j<old.get_size(); j++)
					{
						// cerr<<i<<" "<<j<<endl;
						visit(i, j);
					}
			}
			time++;
		}

		void advance(int time_steps)
		{
				for (int i = 0; i< time_steps; i++)
				{
						advance();
				}
		}

		void visit(int row, int col)
		{
			// cerr<<"called visit "<<row<<" "<<col<<endl;
			double dE = -2* old.compute_point_energy(row, col);
			// cerr<<"de = "<<dE<<endl;
			if (dE < 0 || t_d.flip_q(dE)){
				neu.flip(row, col);
				// cerr<<"Flip"<<endl;
			}
			// cerr<<" "<<t_d.get_temp()<<t_d.flip_q(dE);

		}


};

int main()
{
	simulation s = simulation(20, 1, 1.0, 0.0);
	// lattice l = lattice (40, 1.0, 0.0);
	// l.print();
	// lattice k = lattice (30, 1.0, 0.0);
	// k.print();
	s.advance(500);
	s.neu.print();
	// cerr<<endl<<endl<<"Done"<<endl;
	return 0;
}
