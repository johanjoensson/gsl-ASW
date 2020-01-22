#include "atom.h"

size_t Atom::get_Z() const
{
	return Z;
}

void Atom::set_Z(const size_t Z_n)
{
	this->Z = Z_n;
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
	this->pos = r;
}

void Atom::set_mesh(const Logarithmic_mesh& mesh_n)
{
	this->mesh = mesh_n;
}

GSL::Vector Atom::get_pos() const
{
	return this->pos;
}

double Atom::get_MT() const
{
	return this->MT;
}

double Atom::get_AS() const
{
	return this->AS;
}

/*
Atom::Atom()
 : Z(), MT(), AS(), pos(3), mesh()
{}
*/

Atom::Atom(const double mt, const double as, const size_t z,
	const Logarithmic_mesh &mesh_n, const GSL::Vector &r)
	: Z(z), MT(mt), AS(as), pos(r), mesh(mesh_n)
{}

Atom::Atom(const Logarithmic_mesh &mesh_n, const GSL::Vector &r)
	: Atom(1, 1, 0, mesh_n, r)
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
