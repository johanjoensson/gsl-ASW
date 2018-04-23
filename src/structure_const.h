#ifndef STRUCTURE_CONST_H
#define STRUCTURE_CONST_H
#include <iostream>
#include <gsl/gsl_vector.h>
#include "spherical_fun.h"

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

		friend std::ostream& operator << ( std::ostream&, const Structure_constant& );  
};

std::ostream& operator << ( std::ostream&, const Structure_constant& );  

#endif //STRUCTURE_CONST_H
