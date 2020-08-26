#include <iostream>
#include <algorithm>
#include <cmath>
#include "gaunt.h"
#include "spherical_fun.h"
#include "GSLpp/special_functions.h"


#include <chrono>

bool triangle(lm l1, lm l2, lm l3)
{
	if(l3.l >= std::abs(l1.l - l2.l) && l3.l <= l1.l + l2.l){
		return true;
	}
	return false;
}
GSL::Result gaunt(lm l1, lm l2, lm l3)
{

	if((l1.l + l2.l + l3.l) % 2 != 0){
		return GSL::Result(0, 0);
	}else if(!triangle(l1, l2, l3)){
		return GSL::Result(0, 0);
	}else if(l3.m + l2.m  != l1.m){
		return GSL::Result(0, 0);
	}

	double C = std::sqrt((2.*l1.l + 1)*(2.*l2.l + 1)*(2.*l3.l + 1)/(4*M_PI));
	C = std::abs(l1.m) % 2 == 0 ? C : -C;

	return  C*GSL::wigner_3j(l1.l, l2.l, l3.l,  0,    0,    0  )*
	  			GSL::wigner_3j(l1.l, l2.l, l3.l, -l1.m, l2.m, l3.m);
}

/* Not sure if this is needed anymore...
 * Or indeed if it ever was needed... */
GSL::Result real_gaunt(lm l1, lm l2, lm l3)
{
	int m_max = std::max(l1.m, std::max(l2.m, l3.m));
	int c = 1;
	GSL::Result res;

	if (m_max == l1.m){
		if (l1.m % 2 != 0){
			c = -1;
		}
		res = c*0.5*std::sqrt(2)*gaunt(lm {l1.l, -l1.m}, l2, l3);
	}else if(m_max == l2.m){
		if (l2.m % 2 != 0){
			c = -1;
		}
		res = c*0.5*std::sqrt(2)*gaunt(l1, lm {l2.l, -l2.m}, l3);
	}else if(m_max == l3.m){
		if (l3.m % 2 != 0){
			c = -1;
		}
		res = c*0.5*std::sqrt(2)*gaunt(l1, l2, lm {l3.l, -l3.m});
	}

	return res;
}
