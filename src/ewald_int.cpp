#include "ewald_int.h"
#include "GSLpp/special_functions.h"
#include "GSLpp/basic_math.h"
#include "GSLpp/complex.h"
#include <cmath>

#include <iostream>

double Ewald_integral::bar_ew_int(const double kappa, const double eta, const lm l, const double r) const
{
    const double a = 0.5*sqrt(eta)*r;
    const double a2 = 0.25*eta*r*r;
    const double b = -0.25*kappa*kappa*r*r;
    if(l.l == 0){
        return  M_SQRTPI/4*( GSL::exp(-kappa*r)*GSL::erfc(a - kappa/sqrt(eta)) ).val +
                M_SQRTPI/4*( GSL::exp( kappa*r)*GSL::erfc(a + kappa/sqrt(eta)) ).val;
    }else if(l.l == -1){
	if(std::abs(kappa) > 1e-10){
		return  M_SQRTPI/(2*kappa*r)*(GSL::exp(-kappa*r)*GSL::erfc(a - kappa/sqrt(eta))).val +
			-M_SQRTPI/(2*kappa*r)*(GSL::exp( kappa*r)*GSL::erfc(a + kappa/sqrt(eta))).val;
	}else{
		return (1./a*GSL::exp(-a2) - M_SQRTPI*GSL::erfc(a)).val;
	}
    }
    return   (2.*l.l - 1)/2 * bar_ew_int(kappa, eta, lm {l.l - 1, l.m}, r) +
            -b*bar_ew_int(kappa, eta, lm {l.l - 2, l.m}, r) +
             GSL::pow_int(a, 2*l.l - 1)/2 * GSL::exp(-a2 + b/a2).val;
}

double Ewald_integral::ewald_int(const double kappa, const double eta, const lm l, const double r) const
{
    return bar_ew_int(kappa, eta, l, r)/GSL::pow_int(r, 2*l.l + 1);
}

double Ewald_integral::bar_comp_ew_int(const double kappa, const double eta, const lm l, const double r) const
{
    double a = sqrt(eta)*r/2;
    double a2 = GSL::pow_int(a, 2);
    double b = -kappa*kappa*r*r/4;

    if(l.l == 0){
	    return   M_SQRTPI/4.*(GSL::exp(-kappa*r)*GSL::erfc(-0.5*sqrt(eta)*r + kappa/sqrt(eta))).val +
		    -M_SQRTPI/4.*(GSL::exp(kappa*r)*GSL::erfc(0.5*sqrt(eta)*r + kappa/sqrt(eta))).val;
    }else if(l.l == -1){
		return  M_SQRTPI/(2*kappa*r)*(GSL::exp(-kappa*r)*GSL::erfc(-0.5*sqrt(eta)*r + kappa/sqrt(eta))).val +
			M_SQRTPI/(2*kappa*r)*(GSL::exp(kappa*r)*GSL::erfc(0.5*sqrt(eta)*r + kappa/sqrt(eta))).val;
    }
    return   (2*l.l - 1)/2. * bar_comp_ew_int(kappa, eta, lm {l.l - 1, l.m}, r) +
            -b*bar_comp_ew_int(kappa, eta, lm {l.l - 2, l.m}, r) +
            -0.5*GSL::pow_int(a, 2*l.l - 1) * GSL::exp(-a2 + b/a2).val;
}

double Ewald_integral::comp_ewald_int(const double kappa, const double eta, const lm l, const double r) const
{
    return bar_comp_ew_int(kappa, eta, l, r)/GSL::pow_int(r, 2*l.l+1);
}

std::vector<double> Ewald_integral::evaluate(const double kappa, const double eta, const lm l, const Logarithmic_mesh &mesh) const
{
    double r = 0.;

    std::vector<double> res(mesh.size(), 0.);

    for(unsigned int i = 0; i < mesh.size(); i++){
        r = mesh.r(i);
        res[i] = ewald_int(kappa, eta, l, r);
    }
    return res;
}
std::vector<double> Ewald_integral::evaluate_comp(const double kappa, const double eta, const lm l, const Logarithmic_mesh &mesh) const
{
    double r = 0.;

    std::vector<double> res(mesh.size(), 0.);

    for(unsigned int i = 0; i < mesh.size(); i++){
        r = mesh.r(i);
        res[i] = comp_ewald_int(kappa, eta, l, r);
    }
    return res;
}
