#ifndef lattice_h
#define lattice_h

#include <vector>

char inline symbol(int in) { return in == 0 ? '-' : '*'; }

int inline abs(int in) { return in > 0 ? in : -in; }

int inline to_periodic(int row, int col, int size) {
  if (row < 0 || row >= size)
    row = abs(size - abs(row));
  if (col < 0 || col >= size)
    col = abs(size - abs(col));
  return row * size + col;
}

class lattice {
private:
  int size_;
  std::vector<short> *spins_; //vector<bool> would be more space efficient, but it would not allow multithreading
  double J_;
  double H_;

public:
  lattice() noexcept : size_(0), spins_(NULL), J_(1.0), H_(0.0) {}
  lattice(int new_size, double new_J, double new_H) noexcept
      : size_(new_size), spins_(new std::vector<short>(size_ * size_, 1)),
        J_(new_J), H_(new_H) {}

  lattice(const lattice &other) noexcept
      : lattice(other.size_, other.J_, other.H_) {
#pragma omp parallel for
    for (int i = 0; i < size_ * size_; i++)
      spins_->at(i) = other.spins_->at(i);
  }
  lattice &operator=(const lattice &);

  ~lattice() { delete spins_; }
  void print();
  short get(int row, int col) {
    return spins_->at(to_periodic(row, col, size_));
  }
  int get_size() { return size_; }
  void flip(int row, int col) {
    spins_->at(to_periodic(row, col, size_)) *= -1;
  }
  int total_magnetisation();
  double compute_point_energy(int row, int col);
};

#endif
