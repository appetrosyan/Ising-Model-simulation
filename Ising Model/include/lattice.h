#ifndef lattice_h
#define lattice_h
#include <vector>

class lattice
{
private:
	int size;
	std::vector<short>*  spin;
	double J;
	double H;
public:
	lattice(int _size, double _J, double _H) noexcept;
	lattice(const lattice & old) noexcept;
	lattice() noexcept;
	lattice& operator=(const lattice& );
	~lattice();
	void print();
	short get(int row, int col);
	void flip(int row, int col);
	int get_size();
	int total_magnetisation();
	double compute_point_energy(int row, int col);
};

#endif
