#include <iostream>
#include <algorithm>
#include <cmath>
#include "gaunt.h"
#include "spherical_fun.h"
#include "../../GSL-lib/src/special_functions.h"

#include <gsl/gsl_sf_log.h>

#include <chrono>

// Efficient impolementation of factorial here!
// Using logarithms and continued fractions
double stieltjes_ln_factorial(int n)
{
	// Numerical constants
	double a0 = 1./12, a1 = 1./30, a2 = 53./210, a3 = 195./371;
	double a4 = 22999./22737, a5 = 29944523./19733142;
	double a6 = 109535241009./48264275462;

	double N =  n+1;
	return (1./2)*gsl_sf_log(2*M_PI)+(N-1./2)*gsl_sf_log(N)-N +
	 a0/(N+a1/(N+a2/(N+a3/(N+a4/(N+a5/(N+a6/N)))))) ;

}
// If you want to use a factorial, then use this one!
// Works for n <= 170
double stieltjes_factorial(int n)
{
	unsigned long int shift = 1;
	int m = n;
	while( m < 8){
		shift *= m;
		++m;
	}
	double r = GSL::exp(stieltjes_ln_factorial(n)).val;
	if(n < 8){
		r = (r * n)/(shift * m);
	}
	return r;
}

// In general try to avoid using this one
double triangle_fun(int l1, int l2, int l3)
{
	if (l1 > l2 + l3){
		return 0.;
	}else if (l2 > l3 + l1){
		return 0.;
	}else if (l3 > l1 + l2){
		return 0.;
	}

	double a = stieltjes_factorial(l1 + l2 - l3);
	double b = stieltjes_factorial(l1 - l2 + l3);
	double c = stieltjes_factorial(-l1 + l2 + l3);
	double d = stieltjes_factorial(l1 + l2 + l3 + 1);
	return std::sqrt(a*b*c/d);
}

// Use the gsl-provided version instead
// This one is occasionally just as fast as the GSL one, but mostly it is waay
// slower!
double wigner_3j(lm l1, lm l2, lm l3)
{
	if (l1.m + l2.m + l3.m != 0){
		return 0.;
	}
	int c = 1;
	if ((l1.l - l2.l - l3.m) % 2 != 0){
		c = -1;
	}

	double f1p = 0., f1m = 0., f2p = 0., f2m = 0., f3p = 0., f3m = 0.;

	f1p = stieltjes_factorial(l1.l + l1.m);
	f1m = stieltjes_factorial(l1.l - l1.m);
	f2p = stieltjes_factorial(l2.l + l2.m);
	f2m = stieltjes_factorial(l2.l - l2.m);
	f3p = stieltjes_factorial(l3.l + l3.m);
	f3m = stieltjes_factorial(l3.l - l3.m);
	double a = std::sqrt(f1p*f1m*f2p*f2m*f3p*f3m);

	double sum = 0., recipr = 0.;
	double tmp1 = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0., tmp5 = 0., tmp6 = 0;
	int kmin = std::max(l2.l - l3.l - l1.m, std::max(l1.l - l3.l + l2.m,0));
	int kmax = std::min(l1.l + l2.l - l3.l, std::min(l1.l - l1.m, l2.l + l2.m));
	for (int k = kmin; k <= kmax; k++){
		tmp1 = 1./stieltjes_factorial(k);
		tmp2 = 1./stieltjes_factorial(l1.l + l2.l - l3.l - k);
		tmp3 = 1./stieltjes_factorial(l1.l - l1.m - k);
		tmp4 = 1./stieltjes_factorial(l2.l + l2.m - k);
		tmp5 = 1./stieltjes_factorial(l3.l - l2.l + l1.m + k);
		tmp6 = 1./stieltjes_factorial(l3.l - l1.l - l2.m + k);
		recipr = tmp1*tmp2*tmp3*tmp4*tmp5*tmp6;
		if(k % 2 == 0){
			sum += recipr;
		}else{
			sum -= recipr;
		}
	}
	return triangle_fun(l1.l, l2.l, l3.l)*c*a*sum;
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

	fg = stieltjes_factorial(g);
	if1 = 1./factorial(g - l1.l);
	if2 = 1./factorial(g - l2.l);
	if3 = 1./factorial(g - l3.l);

	double res = c*std::sqrt(f1*f2*f3*if2g)*fg*if1*if2*if3;

	return res;

}

double gaunt(lm l1, lm l2, lm l3)
{
	// Avoid evaluation of the result is going to be zero
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

	double res = C*GSL::wigner_3j(l1.l, l2.l, l3.l, 0, 0, 0).val*
	  GSL::wigner_3j(l1.l, l2.l, l3.l, l1.m, l2.m, l3.m).val;

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
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();


	double dt_gsl = 0, dt_me = 0.;
	double my, gsl = 0;
	for(int i=0;i < 10; i++){
		for(int j = 0; j < 10; j++){
			for(int k = 0; k < 10; k++){
				for(int a = -i; a <= i; a++){
					for(int b = -j; b <= j; b++){
						for(int c = -k; c <= k; c++){

							t1 = std::chrono::high_resolution_clock::now();
							my = wigner_3j(lm {i, a}, lm {j, b}, lm {k, c});
							t2 = std::chrono::high_resolution_clock::now();
							dt_me = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

							t1 = std::chrono::high_resolution_clock::now();
							gsl = GSL::wigner_3j(2*i, 2*j, 2*k, 2*a, 2*b, 2*c).val;
							t2 = std::chrono::high_resolution_clock::now();
							dt_gsl = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

							std::cout << "My function: dt =\t" << dt_me << std::endl;
							std::cout << "GSL function: dt =\t" << dt_gsl << std::endl;
							if(std::abs(my - gsl) > 0.){
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
