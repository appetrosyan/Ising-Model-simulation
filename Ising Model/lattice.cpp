#include <iostream>
#include "lattice.h"
#include <vector>
using namespace std;


lattice::lattice(int _size, double _J, double _H) noexcept
{
		size = _size;
		spin = new short[size * size];
		for(int i=0; i<size*size; i++)
		{
				spin[i] = 1;
		}
		J = _J;
		H = _H;
}

lattice::lattice(const lattice& old) noexcept
{
		size = old.size;
		spin = new short[size*size];
		for(int i=0; i<size*size; i++)
		{
				this->spin[i] = old.spin[i];
		}
		this->J = old.J;
		this->H = old.H;
}

lattice::lattice() noexcept{
	cerr<<"Spin is NULL"<<endl;
	size = 0;
	spin = NULL;
	J = 1.0;
	H = 0.0;
}

lattice::~lattice(){
	cerr<<"De allocating lattice"<<endl;
	// delete[] spin;
}

void lattice::print()
{
		cout<<size<<endl<<endl;
		int area = size*size;
		for(int i=0; i< area; i++)
		{
				cout<<spin[i]<<' ';
				if (i%size ==size-1 ) cout<<endl;
		}
		cout<<endl;
}

short lattice::get(int row, int col)
{
		if(row<0 || row >=size)
		{
				cerr<<"Out of bounds"<<row;
				return 0;
		}
		if(col<0 || col>=size)
		{
				cerr<<"Out of bounds"<<col;
				return 0;
		}
		return spin[row*size+col];
}

void lattice::flip(int row, int col)
{
		if(row<0 || row >=size)
		{
				cerr<<"Out of bounds:"<<row;
				return;
		}
		if(col<0 || col>=size)
		{
				cerr<<"Out of bounds:"<<col;
				return;
		}
		spin[row*size + col] *=-1;
}


int lattice::get_size()
{
	return size;
}

double lattice::compute_point_energy(int row, int col)
{
		int accumulator=0;
		if (row >= size || col >= size)
		{
				cerr<<"Out of bounds"<<endl;
				return 0.0;
		}
		if(row > 0)
		{
				accumulator+=get(row-1, col);
		}
		else
		{
				accumulator+=get(size-1,col);
		}
		if (row < size - 1)
		{
				accumulator+=get(row+1, col);
		}
		else
		{
				accumulator+=get(0, col);
		}
		if (col>0)
		{
				accumulator+=get(row, col-1);
		}
		else
		{
				accumulator+=get(row, size-1);
		}
		if (col< size-1)
		{
				accumulator+=get(row, col+1);
		}
		else
		{
				accumulator+=get(row, 0);
		}
		return -get(row, col)*accumulator*J - H * get(row, col);
}


// int main() {
// 		lattice k = lattice (40, 1.0, 0.0);
// 		// lattice m = lattice (30, 1.0, 0.0);
// 		// m.compute_point_energy(0,0);
// 		// k.print();
// 		// m.print();
// 		// delete &k;
// 		return 0;
// }
