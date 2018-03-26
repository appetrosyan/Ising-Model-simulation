#include <iostream>
#include "include/lattice.h"

using namespace std;

lattice &lattice::operator=(const lattice &other) {
  size_ = other.size_;
  spins_ = other.spins_;
  J_ = other.J_;
  H_ = other.H_;
  delete spins_;
  return *this;
}

void lattice::print() {
  int area = size_ * size_;
  for (int i = 0; i < area; i++) {
    cout << symbol(spins_->at(i));
    if (i % size_ == size_ - 1)
      cout << endl;
  }
  cout << endl;
}



double lattice::compute_point_energy(int row, int col) {
  int accumulator = get(row + 1, col) + get(row - 1, col) + get(row, col - 1) +
                    get(row, col + 1);
  return -get(row, col) * (accumulator * J_ + H_);
}

int lattice::total_magnetisation() {
  int sum = 0;
  #pragma omp parallel for reduction(+ : sum)
  for (int i = 0; i < size_ * size_; i++) {
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
