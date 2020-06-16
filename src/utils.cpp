#include "utils.h"
#include "GSLpp/special_functions.h"
#include "ewald_int.h"
#include <iostream>

double calc_eta(const double vol)
{
	return GSL::exp(GSL::log(6.5) + 2./3*GSL::log(4*M_PI/3) -
			2./3*GSL::log(vol)).val;
}

double bisect_f(const std::function<double(double)>& f, const double tol, const double x_min, const double x_max)
{
	double x = 0;
	double xl = x_min, xu = x_max;
	while(xu - xl > tol){
	    x = (xu + xl)/2;
		if(f(x) > 0){
			xl = x;
		}else{
			xu = x;
		}
	}

	return x;
}

double calc_Rmax(const double vol, const double kappa, const lm &l, const double tol)
{
    double eta = calc_eta(vol);

	Ewald_integral I;
	double r = 1;

    auto f = [=](const double x)
    {
        double res;
        res = (l.l*GSL::log(x) + GSL::log(I.ewald_int(kappa, eta, l, x)) -
        GSL::log(I.ewald_int(kappa, eta, l, 1.)) - GSL::log(tol)).val;
        return res;
    };

	while(f(r) > 0)
	{
		r += 2./sqrt(eta);
	}
	r = bisect_f(f, tol, r - 2./sqrt(eta), r);

	return r;
}

double calc_Kmax(const double vol, const double kappa, const lm &l, const double tol)
{
    double eta = calc_eta(vol);

	double k = 1;

    auto f = [=](const double x)
    {
        double res;
        res = -x*x + eta*(l.l*GSL::log(x) + GSL::log(1 + kappa*kappa) -
        GSL::log(x*x + kappa*kappa) - GSL::log(tol)).val + 1;
        return res;
    };

	while(f(k) > 0){
		k += sqrt(eta);
	}
	k = bisect_f(f, tol, k - sqrt(eta), k);

	return k;
}
