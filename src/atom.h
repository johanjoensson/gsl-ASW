#ifndef ATOM_H
#define ATOM_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
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
	// Nuclear charge
	int Z;
	// Muffin tin radius (MT) &
	// Atomic sphere radius (AS)
	double MT, AS;

	// Atomic position, cartesian coordinates
	gsl_vector *pos;
	// Logarithmic mesh associated with the atom
	Logarithmic_mesh mesh;
	public:
	//! Get nuclear charge
	int get_Z(); 
	//! Set atom position (cartersian)
	void set_pos(gsl_vector &r); 
	//! Set muffin tin radius
	void set_MT(double mt); 
	//! Set atomic sphere radius
	void set_AS(double as); 
	//! Get atomic position (cartesian)
	gsl_vector get_pos();   
	//! Get muffin tin radius
	double get_MT(); 
	//! Get atomic sphere radius
	double get_AS(); 

	Atom(gsl_vector &r, Logarithmic_mesh &mesh);
	Atom(double mt, double as, double z, gsl_vector &r, Logarithmic_mesh &mesh);

	~Atom();
};
#endif //ATOM_H