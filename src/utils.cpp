#include "utils.h"
#include "GSLpp/special_functions.h"
#include "ewald_int.h"
#include <iostream>

std::ostream& operator<<(std::ostream& os, const lm& l)
{
    return os << "(" << l.l << ", " << l.m <<")";
}

bool operator==(const lm &a, const lm &b)
{
    return (a.l == b.l) && (a.m == b.m);
}

bool operator!=(const lm &a, const lm &b)
{
    return !(a == b);
}

double lerp(double x, double x0, double x1, double v0, double v1)
{
    return 1./(x1 - x0) * ((x1 - x)*v0 + (x - x0)*v1);
}

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
	I.set_kappa(kappa);
	I.set_ewald_param(eta);
	double r = 1;

    auto f = [l, I, tol](const double x)
    {
        double res;
        res = (l.l*GSL::log(x) + GSL::log(I.ewald_int(l, x)) -
        GSL::log(I.ewald_int(l, 1.)) - GSL::log(tol)).val;
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

	Ewald_integral I;
	I.set_kappa(kappa);
	I.set_ewald_param(eta);
	double k = 1;

    auto f = [l, I, kappa, eta, tol](const double x)
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
