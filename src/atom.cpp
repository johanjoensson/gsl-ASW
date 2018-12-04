#include "atom.h"

int Atom::get_Z()
{
	return Z;
}

void Atom::set_Z(const int Z)
{
	this->Z = Z;
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

void Atom::set_mesh(const Logarithmic_mesh& mesh)
{
	this->mesh = Logarithmic_mesh(mesh);
}

GSL::Vector Atom::get_pos()
{
	return this->pos;
}

double Atom::get_MT()
{
	return this->MT;
}

double Atom::get_AS()
{
	return this->AS;
}

Atom::Atom()
 : Z(), MT(), AS(), pos(3), mesh()
{}

Atom::Atom(double mt, double as, int z, Logarithmic_mesh &mesh,
	GSL::Vector &r)
	: pos(r), mesh(mesh)
{
	this->MT = mt;
	this->AS = as;
	this->Z = z;
}

Atom::Atom(const double mt, const double as, const int z,
	const Logarithmic_mesh &mesh, const GSL::Vector &r)
	: pos(r), mesh(mesh)
{
	this->MT = mt;
	this->AS = as;
	this->Z = z;
}

Atom::Atom(Logarithmic_mesh &mesh, GSL::Vector &r)
	: Atom(1, 1, 0, mesh, r)
{}

Atom::Atom(const Logarithmic_mesh &mesh, const GSL::Vector &r)
	: Atom(1, 1, 0, mesh, r)
{}

bool operator==(const Atom &a, const Atom &b)
{
	if(a.pos == b.pos && a.Z == b.Z){
		return true;
	}else{
		return false;
	}
}

bool operator!=(const Atom &a, const Atom &b)
{
	return !(a == b);
}
