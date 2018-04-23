#include "atom.h"

int Atom::get_Z()
{
	return Z;
}

void Atom::set_MT(double mt)
{
	this->MT = mt;
}

void Atom::set_AS(double as)
{
	this->AS = as;
}

void Atom::set_pos(gsl_vector &r)
{
	gsl_vector_memcpy(this->pos, &r);
}
gsl_vector Atom::get_pos()
{
	return *(this->pos);
}
	
double Atom::get_MT()
{
	return this->MT;
}

double Atom::get_AS()
{
	return this->AS;
}

Atom::Atom(double mt, double as, double z, gsl_vector &r, Logarithmic_mesh &mesh)
	: mesh(mesh) 
{
	this->MT = mt;
	this->AS = as;
	this->Z = z;
	this->pos = gsl_vector_alloc(3*sizeof(double)); 
	gsl_vector_memcpy(this->pos, &r);
}

Atom::Atom(gsl_vector &r, Logarithmic_mesh &mesh)
	: Atom(1, 1, 0, r, mesh)
{
}

Atom::~Atom()
{
	gsl_vector_free(this->pos);
}
