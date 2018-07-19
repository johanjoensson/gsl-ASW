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

void Atom::set_pos(GSL::Vector &r)
{
	this->pos.copy(r);
}
GSL::Vector Atom::get_pos()
{
	return (this->pos);
}

double Atom::get_MT()
{
	return this->MT;
}

double Atom::get_AS()
{
	return this->AS;
}

Atom::Atom(double mt, double as, double z, GSL::Vector &r, Logarithmic_mesh &mesh)
	: pos(), mesh(mesh)
{
	this->MT = mt;
	this->AS = as;
	this->Z = z;
	this->pos.copy(r);
}

Atom::Atom(GSL::Vector &r, Logarithmic_mesh &mesh)
	: Atom(1, 1, 0, r, mesh)
{
}
