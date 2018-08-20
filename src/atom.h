#ifndef ATOM_H
#define ATOM_H

#include <gsl/gsl_math.h>
#include "../../GSL-lib/src/vector.h"
#include "log_mesh.h"

/***************************************************************************//**
* A class for representing atoms in the cell.\n
* Contains:\n
* __Z__ - Nuclear charge\n
* __AS__ - Atomic sphere radius\n
* __MT__ - Muffin tin radius\n
* __pos__ - position, in cartesian coordinates\n
* __mesh__ - logarithmic mesh to use in intraatomic calculations\n
*******************************************************************************/

class Atom {
	public:
	// Nuclear charge
	int Z;
	// Muffin tin radius (MT) &
	// Atomic sphere radius (AS)
	double MT, AS;

	// Atomic position, cartesian coordinates
	GSL::Vector pos;
	Logarithmic_mesh mesh;
	void set_Z(const int Z);
	//! Get nuclear charge
	int get_Z();
	//! Set atom position (cartersian)
	void set_pos(const GSL::Vector &r);
	//! Set muffin tin radius
	void set_MT(double mt);
	//! Set atomic sphere radius
	void set_AS(double as);
	//! Get atomic position (cartesian)
	GSL::Vector get_pos();
	//! Get muffin tin radius
	double get_MT();
	//! Get atomic sphere radius
	double get_AS();

	Atom();
	Atom(Logarithmic_mesh &mesh, GSL::Vector &r);
	Atom(const Logarithmic_mesh &mesh, const GSL::Vector &r);
	Atom(double mt, double as, double z, Logarithmic_mesh &mesh,
		GSL::Vector &r);
	Atom(const double mt, const double as, const double z,
		const Logarithmic_mesh &mesh, const GSL::Vector &r);

	friend bool operator==(const Atom &a, const Atom &b);
	friend bool operator!=(const Atom &a, const Atom &b);

};

bool operator==(const Atom &a, const Atom &b);
bool operator!=(const Atom &a, const Atom &b);
#endif //ATOM_H
