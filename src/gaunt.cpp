#include <iostream>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_sf_coupling.h>
#include "gaunt.h"
#include "spherical_fun.h"

/* This function should not be used */
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
	return std::sqrt(a*b*c/double(d));
}

/* Use the gsl-provided version instead */
double wigner_3j(lm l1, lm l2, lm l3)
{
	if (l1.m + l2.m + l3.m != 0){
		return 0.;
	}
	int c = 1;
	if ((l1.l - l2.l - l3.m) % 2 != 0){
		c = -1;
	}

	unsigned long int f1p, f1m, f2p, f2m, f3p, f3m = 0.;

	f1p = factorial(l1.l + l1.m);
	f1m = factorial(l1.l - l1.m);
	f2p = factorial(l2.l + l2.m);
	f2m = factorial(l2.l - l2.m);
	f3p = factorial(l3.l + l3.m);
	f3m = factorial(l3.l - l3.m);
	double a = std::sqrt(f1p*f1m*f2p*f2m*f3p*f3m);

	double sum = 0., recipr = 0.;
	unsigned long int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 = 0;
	int kmin = std::max(l2.l - l3.l - l1.m, std::max(l1.l - l3.l + l2.m,0));
	int kmax = std::min(l1.l + l2.l - l3.l, std::min(l1.l - l1.m, l2.l + l2.m));
	for (int k = kmin; k <= kmax; k++){
		tmp1 = factorial(k);
		tmp2 = factorial(l1.l + l2.l - l3.l - k);
		tmp3 = factorial(l1.l - l1.m - k);
		tmp4 = factorial(l2.l + l2.m - k);
		tmp5 = factorial(l3.l - l2.l + l1.m + k);
		tmp6 = factorial(l3.l - l1.l - l2.m + k);
		recipr = 1./(tmp1*tmp2*tmp3*tmp4*tmp5*tmp6);
		if(k % 2 == 0){
			sum += recipr;
		}else{
			sum -= recipr;
		}
	}
	double res = triangle_fun(l1.l, l2.l, l3.l)*c*a*sum;
	return res;
}

/* Use the gsl-provided version instead */
double wigner_3j_0(lm l1, lm l2, lm l3){
	if ((l1.l + l2.l + l3.l) % 2 != 0){
		return 0.;
	}
	int g = (l1.l + l2.l + l3.l)/2;

	int c = 1;
	if (g % 2 != 0){
		c = -1;
	}

	unsigned long int f1, f2, f3, fg, f2g = 0;
	double if2g, if1, if2, if3 = 0;
	f1 = factorial(-l1.l + l2.l + l3.l);
	f2 = factorial(l1.l - l2.l + l3.l);
	f3 = factorial(l1.l + l2.l - l3.l);
	f2g = factorial(l1.l + l2.l + l3.l + 1);
	if2g = 1./f2g;

	fg = factorial(g);
	if1 = 1./factorial(g - l1.l);
	if2 = 1./factorial(g - l2.l);
	if3 = 1./factorial(g - l3.l);

	double res = c*std::sqrt(f1*f2*f3*if2g)*fg*if1*if2*if3;

	return res;

}

double gaunt(lm l1, lm l2, lm l3)
{
	if ((l1.l + l2.l + l3.l) % 2 != 0){
		return 0.;
	}else if(l1.m + l2.m + l3.m != 0){
		return 0.;
	}else if (l1.l > l2.l + l3.l){
		return 0.;
	}else if (l2.l > l3.l + l1.l){
		return 0.;
	}else if (l3.l > l1.l + l2.l){
		return 0.;
	}else if (std::abs(l1.m) > l1.l){
		return 0.;
	}else if (std::abs(l2.m) > l2.l){
		return 0.;
	}else if (std::abs(l3.m) > l3.l){
		return 0.;
	}

	double C = std::sqrt((2*l1.l + 1)*(2*l2.l + 1)*(2*l3.l + 1)/(4*M_PI));

	
//	double res = C*wigner_3j_0(l1, l2, l3)*wigner_3j(l1, l2, l3);
	double res = C*gsl_sf_coupling_3j(2*l1.l, 2*l2.l, 2*l3.l, 0, 0, 0)*gsl_sf_coupling_3j(2*l1.l, 2*l2.l, 2*l3.l, 2*l1.m, 2*l2.m, 2*l3.m);

	return res;
}

/* Not sure if this is needed anymore...
 * Or indeed if it ever was needed... */
double real_gaunt(lm l1, lm l2, lm l3)
{
	int m_max = std::max(l1.m, std::max(l2.m, l3.m));
	int c = 1;
	double res = 0.;

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


/* Proof that the gsl-provided functions shoud be used */
void test_wigner()
{
	double my, gsl = 0;
	for(int i=0;i < 20; i++){
		for(int j = 0; j < 20; j++){
			for(int k = 0; k < 20; k++){
				for(int a = -i; a <= i; a++){
					for(int b = -j; b <= j; b++){
						for(int c = -k; c <= k; c++){
							my = wigner_3j(lm {i, a}, lm {j, b}, lm {k, c});
							gsl = gsl_sf_coupling_3j(2*i, 2*j, 2*k, 2*a, 2*b, 2*c);
							if(my - gsl > 1E-6){
							    std::cout << lm {i, a} << lm {j, b} << lm {k, c}<< "\t" << my;
							    std::cout << "\t" << gsl;
							    std::cout << "\t" << (gsl - my)/gsl << std::endl; 
							}
						}
					}
				}
			}
		}

	}
}

