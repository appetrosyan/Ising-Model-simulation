#include <iostream>
#include "include/lattice.h"
#include <vector>
#include <new>
using namespace std;

int inline to_periodic(int, int, int);

lattice::lattice(int new_size, double new_J, double new_H) noexcept:
size(new_size),
spin(new vector <short> (size*size, 1)),
J(new_J),
H(new_H)
{}

lattice::lattice(const lattice& old) noexcept
{
	size = old.size;
	spin = new vector<short> (size*size, 1);
	#pragma omp for
	for(int i=0; i<size*size; i++)
	{
			this->spin->at(i) = old.spin->at(i);
	}
	J = old.J;
	H = old.H;
}

lattice& lattice::operator= (const lattice& other)
{
	size = other.size;
	delete spin;
	spin = other.spin;
	J = other.J;
	H = other.H;
	return * this;
}

lattice::lattice() noexcept:
size(0),
spin(NULL),
J(0),
H(0)
{}

lattice::~lattice(){
}

char inline symbol(int in)
{
	return in==-1?'-':'*';
}

void lattice::print()
{
	int area = size*size;
	for(int i=0; i< area; i++)
	{
		cout<<symbol(spin->at(i))<<' ';
		if (i%size ==size-1 ) cout<<endl;
	}
	cout<<endl;
}


short lattice::get(int row, int col)
{
	return spin->at(to_periodic(row, col, size));
}

void lattice::flip(int row, int col)
{
	spin->at(to_periodic(row, col, size)) *=-1;
}


int lattice::get_size()
{
	return size;
}

double lattice::compute_point_energy(int row, int col)
{
	int accumulator=get(row+1, col)+get(row-1, col) + get(row, col-1) + get(row, col+1);
	return -get(row, col)*(accumulator*J + H );
}

int inline abs(int in){
	return in>0? in:-in;
}

int inline to_periodic(int row, int col, int sz){
	if(row<0 || row>=sz){
		row = abs(sz -abs(row));
	}
	if(col<0 || col>=sz){
		col = abs(sz -abs(col));
	}
	return row* sz + col;
}




// int main() {
// 	int len = 5;
// 	lattice k = lattice (len, 1.0, 0.0);
// 	// for(int i =0; i<len; i++){
// 	// 	k.flip(i,i);
// 	// }
// 	k.flip(len/2, len/2);
// 	k.flip(1,6);
// 	cout<<k.compute_point_energy(0,0)<<endl;
// 	cout<<k.compute_point_energy(1,1)<<endl;
// 	k.print();
//
// 	return 0;
// }
