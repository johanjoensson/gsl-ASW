#ifndef ATOM_H
#define ATOM_H

#include <gsl/gsl_math.h>
#include <GSLpp/vector.h>
#include <log_mesh.h>

/*******************************************************************************
* A class for representing atoms in the cell.\n
* @param Z Nuclear charge\n
* @param AS Atomic sphere radius\n
* @param MT Muffin tin radius\n
* @param pos Position, in cartesian coordinates\n
* @param mesh Logarithmic mesh to use in intraatomic calculations\n
*******************************************************************************/
class Atom {
private:
	// Nuclear charge
	size_t Z_m;
	// Muffin tin radius (MT) &
	// Atomic sphere radius (AS)
	double MT_m, AS_m;

public:
	//! Get nuclear charge
	size_t Z() const {return Z_m;}
	size_t& Z() {return Z_m;}
	//! Get muffin tin radius
	double MT() const {return MT_m;}
	double& MT() {return MT_m;}
	//! Get atomic sphere radius
	double AS() const {return AS_m;}
	double& AS() {return AS_m;}

	Atom() = default;
	Atom(const Atom&) = default;
	Atom(Atom&&) = default;

	Atom(const double mt, const double as, const size_t z)
	 : Z_m(z), MT_m(mt), AS_m(as)
	{}


	~Atom() = default;

	Atom& operator=(const Atom&) = default;
	Atom& operator=(Atom&&) = default;

	friend bool operator==(const Atom &a, const Atom &b);
	friend bool operator!=(const Atom &a, const Atom &b);

};

inline bool operator==(const Atom &a, const Atom &b)
{
	return a.MT() == b.MT() && a.AS() == b.AS() && a.Z() == b.Z();
}

inline bool operator!=(const Atom &a, const Atom &b)
{
	return !(a == b);
}

#endif //ATOM_H
