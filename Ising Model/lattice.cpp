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
  short top = get(row? row-1: size_-1, col);
  short bottom = get((row+1)%size_, col);
  short left = get(row, col? col-1: size_-1);
  short right = get(row, (col+1)%size_) ;
  int accumulator = top + bottom +left +right;
  return -get(row, col) * (accumulator * J_ + H_);
}

/*
Computes total magnetisation in O(n^2). Thread safe
*/
int lattice::total_magnetisation() {
  int sum = 0;
  #pragma omp target teams distribute parallel for reduction(+ : sum)
  for (unsigned int i = 0; i < size_ * size_; i++) {
    sum += spins_->at(i);
  }
  return sum;
}
