#ifndef STRUCTURE_CONST_H
#define STRUCTURE_CONST_H
#include <iostream>
//#include <gsl/gsl_vector.h>
#include "utils.h"
#include "bloch_sum.h"
#include "GSLpp/vector.h"
#include "GSLpp/complex.h"

/***************************************************************************//**
* A class for representing structure constants\n
* Contains:\n
* __l_int__ - Maximum orbital angular momentum to be included for Bessel
* expansions\n __l_low__ - Maximumm orbital angular momentum to be included for
* Hankel expansions\n __l1__, __l2__ - Orbital angular momenta to couple via the
* structure constant\n
* __kappa__ - Energy parameter used (usually kappa^2 = -0.015)\n
* __r__ - Position of atom\n
* __val__ - Value of the structure constant\n
* __dk_val__ - Value of energy derivative of the structure constant\n
*******************************************************************************/
class Structure_constant{
		int l_int;
	    lm l1, l2;
		double kappa;
	public:
		Structure_constant(int l_int_n, double kappa_n, lm l1_n, lm l2_n);
		Structure_constant(int l_int_n, lm l1_n, lm l2_n);
		Structure_constant(lm l1_n, lm l2_n);
		Structure_constant();

		double operator()(const GSL::Vector& r) const;
		double dot(const GSL::Vector& r) const;
		friend std::ostream& operator << ( std::ostream&,
			const Structure_constant& );
};

class Bloch_summed_structure_constant{
		int l_int;
		double kappa;
		Crystal c;
	    lm l1, l2;
	public:
		Bloch_summed_structure_constant():l_int(), kappa(), c(), l1(), l2()
			{};
		Bloch_summed_structure_constant(int l_int, double kappa,
			const Crystal& c, lm l1, lm l2);
		Bloch_summed_structure_constant(int l_int, const Crystal& c, lm l1,
			lm l2);
		Bloch_summed_structure_constant(const Crystal& c, lm l1, lm l2);

		GSL::Complex operator()(const GSL::Vector& tau, const GSL::Vector& kp)const;
		GSL::Complex dot(const GSL::Vector& tau, const GSL::Vector& kp)const;

		friend std::ostream& operator << ( std::ostream&,
			const Bloch_summed_structure_constant& );
};

std::ostream& operator << ( std::ostream&, const Structure_constant& );
std::ostream& operator << ( std::ostream&,
	const Bloch_summed_structure_constant& );

#endif //STRUCTURE_CONST_H
