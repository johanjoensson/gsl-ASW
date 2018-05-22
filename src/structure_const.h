#ifndef STRUCTURE_CONST_H
#define STRUCTURE_CONST_H
#include <iostream>
#include <gsl/gsl_vector.h>
#include "spherical_fun.h"

/***************************************************************************//**
* A class for representing structure constants\n
* Contains:\n
* __l_int__ - Maximum orbital angular momentum to be included for Bessel expansions\n
* __l_low__ - Maximumm orbital angular momentum to be included for Hankel expansions\n
* __l1__, __l2__ - Orbital angular momenta to couple via the structure constant\n
* __kappa__ - Energy parameter used (usually kappa^2 = -0.015)\n
* __r__ - Position of atom\n
* __val__ - Value of the structure constant\n
* __dk_val__ - Value of energy derivative of the structure constant\n
*******************************************************************************/
class Structure_constant{
		int l_int, l_low;
	        lm l1, l2;
		double kappa;
		gsl_vector r;
	public:
		Structure_constant(int l_low, int l_int, double kappa, lm l1, lm l2, gsl_vector r);
		Structure_constant(int l_low, int l_int, lm l1, lm l2, gsl_vector r);
		Structure_constant(int l_low, lm l1, lm l2, gsl_vector r);
		Structure_constant(lm l1, lm l2, gsl_vector r);

		double val;
		double dk_val;

		//! Functionality for printing structure constants
		friend std::ostream& operator << ( std::ostream&, const Structure_constant& );
};

std::ostream& operator << ( std::ostream&, const Structure_constant& );

#endif //STRUCTURE_CONST_H
