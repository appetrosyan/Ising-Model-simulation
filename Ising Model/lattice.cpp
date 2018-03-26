#include <iostream>
#include "include/lattice.h"

using namespace std;

/*
Copy assignment operator, too long to include in the header.
*/
lattice &lattice::operator=(const lattice &other) {
  size_ = other.size_;
  spins_ = other.spins_;
  J_ = other.J_;
  H_ = other.H_;
  delete spins_;
  return *this;
}

void lattice::print() {
  unsigned int area = size_ * size_;
  for (unsigned int i = 0; i < area; i++) {
    cout << to_symbol(spins_->at(i));
    if (i % size_ == size_ - 1)
      cout << endl;
  }
  cout << endl;
}

/*
Computes the energy associated with a spin at the given point.

It is explicitly float as that would allow the compiler to make use of multiple
registers instead of keeping track of unneeded precision.  (typically J, H ~ 1).
*/
float lattice::compute_point_energy(int row, int col) {
  int accumulator = get(row + 1, col) + get(row - 1, col) + get(row, col - 1) +
                    get(row, col + 1);
  return -get(row, col) * (accumulator * J_ + H_);
}

/*
Computes total magnetisation in O(n^2). Thread safe
*/
int lattice::total_magnetisation() {
  int sum = 0;
  #pragma omp parallel for reduction(+ : sum)
  for (unsigned int i = 0; i < size_ * size_; i++) {
    sum += spins_->at(i);
  }
  return sum;
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
