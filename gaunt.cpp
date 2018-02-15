#include <iostream>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "gaunt.h"
#include "spherical_fun.h"

double triangle_fun(int l1, int l2, int l3)
{
	if (l1 > l2 + l3){
		return 0.;
	}else if (l2 > l3 + l1){
		return 0.;
	}else if (l3 > l1 + l2){
		return 0.;
	}

	unsigned long int a = factorial(l1 + l2 - l3);
	unsigned long int b = factorial(l1 - l2 + l3);
	unsigned long int c = factorial(-l1 + l2 + l3);
	unsigned long int d = factorial(l1 + l2 + l3 + 1);
	return sqrt(a*b*c/double(d));
}

double wigner_3j(int l1, int l2, int l3, int m1, int m2, int m3)
{
	if (m1 + m2 + m3 != 0){
		return 0.;
	}
	int c = 1;
	if ((l1 - l2 - m3) % 2 != 0){
		c = -1;
	}

	unsigned long int f1p, f1m, f2p, f2m, f3p, f3m = 0.;

	f1p = factorial(l1 + m1);
	f1m = factorial(l1 - m1);
	f2p = factorial(l2 + m2);
	f2m = factorial(l2 - m2);
	f3p = factorial(l3 + m3);
	f3m = factorial(l3 - m3);
	double a = sqrt(f1p*f1m*f2p*f2m*f3p*f3m);

	double sum, recipr = 0;
	unsigned long int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 = 0;
	int kmin = std::max(l2 - l3 - m1, std::max(l1 - l3 + m2,0));
	int kmax = std::min(l1 + l2 - l3, std::min(l1 - m1, l2 + m2));
	for (int k = kmin; k <= kmax; k++){
		tmp1 = factorial(k);
		tmp2 = factorial(l1+l2-l3-k);
		tmp3 = factorial(l1-m1-k);
		tmp4 = factorial(l2+m2-k);
		tmp5 = factorial(l3-l2+m1+k);
		tmp6 = factorial(l3-l1-m2+k);
		recipr = 1./(tmp1*tmp2*tmp3*tmp4*tmp5*tmp6);
		if(k % 2 == 0){
			sum += recipr;
		}else{
			sum -= recipr;
		}
	}
	double res = triangle_fun(l1,l2,l3)*c*a*sum;
//	std::cout << "wigner_3j(" << l1 <<", " << l2 << ", " << l3 << ", " << m1 << ", " << m2 << ", " << m3 << ") = " << res << std::endl;
	return res;
}

double wigner_3j_0(int l1, int l2, int l3){
	if ((l1 + l2 + l3) % 2 != 0){
		return 0.;
	}
	int g = (l1 + l2 + l3)/2;

	int c = 1;
	if (g % 2 != 0){
		c = -1;
	}

	unsigned long int f1, f2, f3, fg, f2g = 0;
	double if2g, if1, if2, if3 = 0;
	f1 = factorial(-l1 + l2 + l3);
	f2 = factorial(l1 - l2 + l3);
	f3 = factorial(l1 + l2 - l3);
	f2g = factorial(l1 + l2 + l3 + 1);
	if2g = 1./f2g;

	fg = factorial(g);
	if1 = 1./factorial(g - l1);
	if2 = 1./factorial(g - l2);
	if3 = 1./factorial(g - l3);

	double res = sqrt(f1*f2*f3*if2g)*fg*if1*if2*if3*c;

	return res;

}

double gaunt(int l1, int l2, int l3, int m1, int m2, int m3)
{
	if ((l1 + l2 + l3) % 2 != 0){
		return 0.;
	}else if(m1 + m2 + m3 != 0){
		return 0.;
	}else if (l1 > l2 + l3){
		return 0.;
	}else if (l2 > l3 + l1){
		return 0.;
	}else if (l3 > l1 + l2){
		return 0.;
	}else if (abs(m1) > l1){
		return 0.;
	}else if (abs(m2) > l2){
		return 0.;
	}else if (abs(m3) > l3){
		return 0.;
	}

	int c = 1;
	if (m1 % 2 != 0) {
		c = -1;
	}
	double C = sqrt((2*l1 + 1)*(2*l2 + 1)*(2*l3 + 1)/(4*M_PI));

	double res = c*C*wigner_3j_0(l1, l2, l3)*wigner_3j(l1, l2, l3, -m1, m2, m3);


	return res;
}
double real_gaunt(int l1, int l2, int l3, int m1, int m2, int m3)
{
	int m_max = std::max(m1, std::max(m2, m3));
	int c = 1;
	double res = 0.;

	if (m_max == m1){
		if (m1 % 2 != 0){
			c = -1;
		}
		res = c*0.5*sqrt(2)*gaunt(l1, l2, l3, -m1, m2, m3);
	}else if(m_max == m2){
		if (m2 % 2 != 0){
			c = -1;
		}
		res = c*0.5*sqrt(2)*gaunt(l1, l2, l3, m1, -m2, m3);
	}else if(m_max == m3){
		if (m3 % 2 != 0){
			c = -1;
		}
		res = c*0.5*sqrt(2)*gaunt(l1, l2, l3, m1, m2, -m3);
	}

	return res;
}

/*
int main()
{
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			for(int k = 0; k < 4; k++){
				for(int m1 = -i; m1 <= i; m1++){
				for(int m2 = -j; m2 <= j; m2++){
				for(int m3 = -k; m3 <= k; m3++){
					std::cout << "gaunt((" << i << ", " << m1 << "), (" << j << ", " << m2 << "), (" << k << ", " << m3 << ")) = " << gaunt(i,j,k,m1,m2,m3) << std::endl;
					std::cout << "real_gaunt((" << i << ", " << m1 << "), (" << j << ", " << m2 << "), (" << k << ", " << m3 << ")) = " << real_gaunt(i,j,k,m1,m2,m3) << std::endl;
				}
				}
				}
			}
		}
	}

	return 0;
}
*/
