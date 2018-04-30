#ifndef lattice_h
#define lattice_h

#include <cmath>
#include <vector>

/* Converts spin up/down to easily printable symbols. */
char inline to_symbol(int in) { return in == -1 ? '-' : '+'; }

/* Converts given pair of indices to those with periodic boundary conditions. */
int inline to_periodic(int row, int col, int size) {
  if (row < 0 || row >= size)
    row = abs(size - abs(row));
  if (col < 0 || col >= size)
    col = abs(size - abs(col));
  return row * size + col;
}

class lattice {
private:
  unsigned int size_;
  // vector<bool> would be more space efficient, but it would not allow
  // multithreading
  std::vector<short> *spins_;
  float J_;
  float H_;

public:
  lattice() noexcept : size_(0), spins_(NULL), J_(1.0), H_(0.0) {}
  lattice(int new_size, double new_J, double new_H) noexcept
      : size_(new_size), spins_(new std::vector<short>(size_ * size_, 1)),
        J_(new_J), H_(new_H) {}
  lattice(const lattice &other) noexcept
      : lattice(other.size_, other.J_, other.H_) {
    for (unsigned int i = 0; i < size_ * size_; i++)
      spins_->at(i) = other.spins_->at(i);
  }
  lattice &operator=(const lattice &);

  ~lattice() { delete spins_; }
  void print() ;
  short get(int row, int col) const{
    return spins_->at(to_periodic(row, col, size_));
  }
  unsigned int get_size()  { return size_; }
  void flip(int row, int col) { spins_->at(to_periodic(row, col, size_)) *= -1; }
  int total_magnetisation() ;
  float compute_point_energy(int row, int col);
};

#endif
