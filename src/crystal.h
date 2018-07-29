#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <vector>
#include "atom.h"
#include "../../GSL-lib/src/matrix.h"
#include "../../GSL-lib/src/vector.h"

class Crystal {

	void calc_Kn(size_t num);
	void calc_Rn(size_t num);

	double calc_volume();
public:
	std::vector<GSL::Vector> Rn_vecs;
	std::vector<GSL::Vector> Kn_vecs;
	double scale;
	GSL::Matrix lattice;
	std::vector<Atom> atoms;
	double volume;

	Crystal();
	Crystal(double& a);
	Crystal(const double& a);
	Crystal(double& a, double& b, double& c);
	Crystal(const double& a, const double& b, const double& c);
	Crystal(GSL::Vector& a, GSL::Vector& b, GSL::Vector& c);
	Crystal(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c);

	void add_atoms(const std::vector<Atom>& v);

	GSL::Vector& get_Rn(const size_t& i);
};
#endif //CRYSTAL_H
