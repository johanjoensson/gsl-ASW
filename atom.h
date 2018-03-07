#ifndef ATOM_H
#define ATOM_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include "log_mesh.h"

class Atom {
	// Nuclear charge
	int Z;
	// Muffin tin radius (MT) &
	// Atomic sphere radius (AS)
	double MT, AS;

	gsl_vector *pos;
	Logarithmic_mesh mesh;
	public:
	int get_Z();
	gsl_vector get_pos(gsl_vector pos);
	void set_MT(double mt);
	void set_AS(double as);

	gsl_vector get_pos();
	double get_MT();
	double get_AS();

	Atom(gsl_vector &r, Logarithmic_mesh &mesh);
	Atom(double mt, double as, double z, gsl_vector &r, Logarithmic_mesh &mesh);

	~Atom();
};
#endif //ATOM_H
