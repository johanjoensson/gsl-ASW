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
	public:
	// Nuclear charge
	size_t Z;
	// Muffin tin radius (MT) &
	// Atomic sphere radius (AS)
	double MT, AS;

	// Atomic position, cartesian coordinates
	GSL::Vector pos;
//	Logarithmic_mesh mesh;
	void set_Z(const size_t Z_n);
	//! Get nuclear charge
	size_t get_Z() const;
	//! Set atom position (cartersian)
	void set_pos(const GSL::Vector &r);
	//! Set muffin tin radius
	void set_MT(double mt);
	//! Set atomic sphere radius
	void set_AS(double as);
	void set_mesh(const Logarithmic_mesh& mesh);
	//! Get atomic position (cartesian)
	GSL::Vector get_pos() const;
	//! Get muffin tin radius
	double get_MT() const;
	//! Get atomic sphere radius
	double get_AS() const;


/*
	Atom() = default;
	Atom(const Atom&) = default;
	Atom(Atom&&) = default;
	~Atom() = default;
*/
	Atom(/*const Logarithmic_mesh &mesh,*/ const GSL::Vector &r);
	Atom(const double mt, const double as, const size_t z,
		/*const Logarithmic_mesh &mesh,*/ const GSL::Vector &r);


	// Atom& operator=(const Atom&) = default;
	// Atom& operator=(Atom&&) = default;
	friend bool operator==(const Atom &a, const Atom &b);
	friend bool operator!=(const Atom &a, const Atom &b);

};

bool operator==(const Atom &a, const Atom &b);
bool operator!=(const Atom &a, const Atom &b);
#endif //ATOM_H
