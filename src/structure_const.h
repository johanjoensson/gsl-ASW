#ifndef STRUCTURE_CONST_H
#define STRUCTURE_CONST_H
#include <iostream>
#include <gsl/gsl_vector.h>
#include "utils.h"
#include "bloch_sum.h"
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/complex.h"

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
		int l_int, l_low;
	        lm l1, l2;
		double kappa;
		GSL::Vector r;
	public:
		Structure_constant(int l_low, int l_int, double kappa, lm l1, lm l2,
			GSL::Vector& r);
		Structure_constant(int l_low, int l_int, lm l1, lm l2, GSL::Vector& r);
		Structure_constant(int l_low, lm l1, lm l2, GSL::Vector& r);
		Structure_constant(lm l1, lm l2, GSL::Vector& r);
		Structure_constant();

		double val;
		double dk_val;

		friend std::ostream& operator << ( std::ostream&,
			const Structure_constant& );
};

class Bloch_summed_structure_constant{
		int l_low, l_int;
		double kappa;
		Crystal c;
	    lm l1, l2;
	public:
		Bloch_summed_structure_constant(int l_low, int l_int, double kappa,
			Crystal& c, lm l1, lm l2);
		Bloch_summed_structure_constant(int l_low, int l_int, Crystal& c, lm l1,
			lm l2);
		Bloch_summed_structure_constant(int l_low, Crystal& c, lm l1, lm l2);
		Bloch_summed_structure_constant(Crystal& c, lm l1, lm l2);
		Bloch_summed_structure_constant();

		GSL::Complex evaluate(const GSL::Vector& tau, const GSL::Vector& kp);
		GSL::Complex dot_evaluate(const GSL::Vector& tau,
			const GSL::Vector& kp);

		friend std::ostream& operator << ( std::ostream&,
			const Bloch_summed_structure_constant& );
};

std::ostream& operator << ( std::ostream&, const Structure_constant& );
std::ostream& operator << ( std::ostream&,
	const Bloch_summed_structure_constant& );

#endif //STRUCTURE_CONST_H
