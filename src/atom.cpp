#include "atom.h"

int Atom::get_Z()
{
	return Z;
}

void Atom::set_MT(const double mt)
{
	this->MT = mt;
}

void Atom::set_AS(const double as)
{
	this->AS = as;
}

void Atom::set_pos(const GSL::Vector &r)
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

Atom::Atom(double mt, double as, double z, GSL::Vector &r)
	: pos(r)
{
	this->MT = mt;
	this->AS = as;
	this->Z = z;
}

Atom::Atom(const double mt, const double as, const double z, const GSL::Vector &r)
	: pos(r)
{
	this->MT = mt;
	this->AS = as;
	this->Z = z;
}

Atom::Atom(GSL::Vector &r)
	: Atom(1, 1, 0, r)
{
}

Atom::Atom(const GSL::Vector &r)
	: Atom(1, 1, 0, r)
{
}
