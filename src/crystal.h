#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <vector>
#include "lattice.h"
#include "atom.h"
#include "utils.h"
#include "../../GSL-lib/src/matrix.h"
#include "../../GSL-lib/src/vector.h"

class Crystal {

	double calc_volume();
	double calc_bz_volume();
public:
	std::vector<GSL::Vector> Rn_vecs;
	std::vector<GSL::Vector> Kn_vecs;
	Lattice lat;
	std::vector<Atom> atoms;

	size_t calc_nk(double tol, double kappa, lm l);
	size_t calc_nr(double tol, double kappa, lm l);

	double calc_Kmax(double tol, double kappa, lm l);
	double calc_Rmax(double tol, double kappa, lm l);

	void set_Kn(double Kmax);
	void set_Rn(double Rmax);

	void calc_Kn(size_t num);
	void calc_Rn(size_t num);

	double volume, bz_volume;

	Crystal();
	Crystal(const double& a);
	Crystal(const double& a, const double& b, const double& c);
	Crystal(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c);

	void add_atoms(const std::vector<Atom>& v);

	GSL::Vector& get_Rn(const size_t& i);
	std::vector<std::vector<Atom>> calc_nearest_neighbours();
};
#endif //CRYSTAL_H
