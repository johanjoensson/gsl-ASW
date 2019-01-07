#include "ewald_int.h"
#include "GSLpp/special_functions.h"
#include "GSLpp/basic_math.h"
#include <cmath>

#include <iostream>

Ewald_integral::Ewald_integral()
 : ewald_param(1.), kappa(sqrt(0.015))
{}

void Ewald_integral::set_ewald_param(double eta)
{
	this->ewald_param = eta;
}

void Ewald_integral::set_kappa(double kappa_n)
{
	this->kappa = kappa_n;
}

double Ewald_integral::bar_ew_int(lm l, double r)
{

    double res = 0;

    double eta = this->ewald_param;
    double t1 = 0., t2 = 0., t3 = 0.;
    double a = sqrt(eta)*r/2.;
    double b = -kappa*kappa*r*r/4.;
    if(l.l == 0){
        t1 = GSL::exp(-kappa*r).val + GSL::exp(kappa*r).val;
        t2 = GSL::exp(-kappa*r).val*GSL::erf(0.5*sqrt(eta)*r - kappa/sqrt(eta)).val;
        t3 = GSL::exp(kappa*r).val*GSL::erf(0.5*sqrt(eta)*r + kappa/sqrt(eta)).val;
        // t1 = GSL::exp(-kappa*r).val*GSL::erfc(1./2. *sqrt(eta)*r - kappa/sqrt(eta)).val;
        // t2 = GSL::exp(kappa*r).val*GSL::erfc(1./2. *sqrt(eta)*r + kappa/sqrt(eta)).val;

        res = M_SQRTPI/4*(t1 - t2 - t3);
    }else if(l.l == -1){
        t1 = (GSL::exp(-kappa*r) - GSL::exp(kappa*r)).val;
        t2 = GSL::exp(-this->kappa*r).val*GSL::erf(sqrt(eta)*r/2. - this->kappa/sqrt(eta)).val;
        t3 = GSL::exp(this->kappa*r).val*GSL::erf(sqrt(eta)*r/2. + this->kappa/sqrt(eta)).val;
        // t1 = GSL::exp(-kappa*r).val*GSL::erfc(1./2. *sqrt(eta)*r - kappa/sqrt(eta)).val;
        // t2 = GSL::exp(kappa*r).val*GSL::erfc(1./2. *sqrt(eta)*r + kappa/sqrt(eta)).val;

        res = M_SQRTPI/(2*kappa*r)*(t1 - t2 + t3);
    }else if(l.l >= 1){
        t1 = (2*l.l - 1)/2. * bar_ew_int(lm {l.l - 1, l.m}, r);
        t2 = -b*bar_ew_int(lm {l.l - 2, l.m}, r);
        t3 = GSL::pow_int(a, (2*l.l-1))/2. * GSL::exp(-a*a + b/(a*a)).val;
        res = t1 + t2 + t3;
    }
    return res;
}

double Ewald_integral::ewald_int(lm l, double r)
{
    return bar_ew_int(l, r)/GSL::pow_int(r, 2*l.l + 1);
}

double Ewald_integral::bar_comp_ew_int(lm l, double r)
{
    double res = 0.;
    double eta = this->ewald_param;
    double t1 = 0., t2 = 0., t3 = 0.;
    double a = sqrt(eta)*r/2;
    double b = -kappa*kappa*r*r/4;

    if(l.l == 0){
        t1 = GSL::exp(-kappa*r).val - GSL::exp(kappa*r).val;
        t2 = GSL::exp(-kappa*r).val*GSL::erf(-0.5*sqrt(eta)*r + kappa/sqrt(eta)).val;
        t3 = GSL::exp(kappa*r).val*GSL::erf(0.5*sqrt(eta)*r + kappa/sqrt(eta)).val;

        res = M_SQRTPI/4 * (t1 - t2 + t3);
    }else if(l.l == -1){
        t1 = GSL::exp(-kappa*r).val + GSL::exp(kappa*r).val;
        t2 = GSL::exp(-kappa*r).val*GSL::erf(-0.5*sqrt(eta)*r + kappa/sqrt(eta)).val;
        t3 = GSL::exp(kappa*r).val*GSL::erf(0.5*sqrt(eta)*r + kappa/sqrt(eta)).val;

        res = M_SQRTPI/(2*kappa*r) * (t1 - t2 - t3);
    }else if(l.l >= 1){
        t1 = (2*l.l - 1)/2. * bar_comp_ew_int(lm {l.l - 1, l.m}, r);
        t2 = -b*bar_comp_ew_int(lm {l.l - 2, l.m}, r);
        t3 = -0.5*GSL::pow_int(a, (2*l.l-1)) * GSL::exp(-a*a + b/(a*a)).val;
        res = t1 + t2 + t3;
    }
    return res;
}

double Ewald_integral::comp_ewald_int(lm l, double r)
{
    return bar_comp_ew_int(l, r)/GSL::pow_int(r, 2*l.l+1);
}

std::vector<double> Ewald_integral::evaluate(lm l, Logarithmic_mesh &mesh)
{
    double r = 0.;

    std::vector<double> res(mesh.r.size(), 0.);

    for(unsigned int i = 0; i < mesh.r.size(); i++){
        r = mesh.r[i];
        res[i] = ewald_int(l,r);
    }
    return res;
}
std::vector<double> Ewald_integral::evaluate_comp(lm l, Logarithmic_mesh &mesh)
{
    double r = 0.;

    std::vector<double> res(mesh.r.size(), 0.);

    for(unsigned int i = 0; i < mesh.r.size(); i++){
        r = mesh.r[i];
        res[i] = comp_ewald_int(l,r);
    }
    return res;
}
