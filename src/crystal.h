#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <vector>
#include "lattice.h"
#include "atom.h"
#include "utils.h"
#include "../../GSL-lib/src/matrix.h"
#include "../../GSL-lib/src/vector.h"

/***************************************************************************//**
* Class used to describe the structure of a crystal.\n
* Contains:\n
* __Rn_vecs__ - Lattice vectors of the crystal (up to a maximum lentgh).\n
* __Kn_vecs__ - Reciprocal lattice vectors of the crystal (up to a maximum length).\n
* __lat__ - The lattice of the crystal.\n
* __atoms__ - The atoms in the unit cell of the crystal (the base of the crystal).\n
*******************************************************************************/
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

	void calc_Kn(size_t num);
	void calc_Rn(size_t num);

	double volume, bz_volume;

	Crystal();
	Crystal(double& a);
	Crystal(const double& a);
	Crystal(double& a, double& b, double& c);
	Crystal(const double& a, const double& b, const double& c);
	Crystal(GSL::Vector& a, GSL::Vector& b, GSL::Vector& c);
	Crystal(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c);

	void add_atoms(const std::vector<Atom>& v);

	GSL::Vector& get_Rn(const size_t& i);
	std::vector<std::vector<Atom>> calc_nearest_neighbours();
};
#endif //CRYSTAL_H
