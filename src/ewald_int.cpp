#include "ewald_int.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <iostream>

Ewald_integral::Ewald_integral()
{
    this->kappa = sqrt(0.015);
    this->ewald_param = 1;
}

void Ewald_integral::set_ewald_param(double eta)
{
	this->ewald_param = eta;
}

void Ewald_integral::set_kappa(double kappa)
{
	this->kappa = kappa;
}

double Ewald_integral::bar_ew_int(lm l, double r){

    double res = 0;

    double eta = this->ewald_param, kappa = this->kappa;
    double t1 = 0., t2 = 0., t3 = 0.;
    double a = sqrt(eta)*r/2.;
    double b = -kappa*kappa*r*r/4.;
    if(l.l == 0){
        t1 = exp_gsl(-kappa*r) + exp_gsl(kappa*r);
        t2 = exp_gsl(-kappa*r)*gsl_sf_erf(0.5*sqrt(eta)*r - kappa/sqrt(eta));
        t3 = exp_gsl(kappa*r)*gsl_sf_erf(0.5*sqrt(eta)*r + kappa/sqrt(eta));

        res = M_SQRTPI/4*(t1 - t2 - t3);
    }else if(l.l == -1){
        t1 = exp_gsl(-kappa*r) - exp_gsl(kappa*r);
        t2 = exp_gsl(-this->kappa*r)*gsl_sf_erf(sqrt(eta)*r/2 - this->kappa/sqrt(eta));
        t3 = exp_gsl(this->kappa*r)*gsl_sf_erf(sqrt(eta)*r/2 + this->kappa/sqrt(eta));

        res = M_SQRTPI/(2*kappa)*(t1 - t2 + t3)/r;
    }else if(l.l >= 1){
        t1 = (2*l.l - 1)/2. * bar_ew_int(lm {l.l - 1, l.m}, r);
        t2 = -b*bar_ew_int(lm {l.l - 2, l.m}, r);
        if(-a*a + b/(a*a) > -700){
            t3 = gsl_pow_int(a, (2*l.l-1))/2 * exp_gsl(-a*a + b/(a*a));
        }else{
            t3 = 0;
        }
        res = t1 + t2 + t3;
    }
    return res;
}

double Ewald_integral::ewald_int(lm l, double r)
{
    return bar_ew_int(l, r)/gsl_pow_int(r, 2*l.l + 1);
}

double Ewald_integral::bar_comp_ew_int(lm l, double r){
    double res = 0.;
    double eta = this->ewald_param;
    double kappa = this->kappa;
    double t1 = 0., t2 = 0., t3 = 0.;
    double a = sqrt(eta)*r/2;
    double b = -kappa*kappa*r*r/4;

    if(l.l == 0){
        t1 = exp_gsl(-kappa*r) - exp_gsl(kappa*r);
        t2 = exp_gsl(-kappa*r)*gsl_sf_erf(-0.5*sqrt(eta)*r + kappa/sqrt(eta));
        t3 = exp_gsl(kappa*r)*gsl_sf_erf(0.5*sqrt(eta)*r + kappa/sqrt(eta));

        res = M_SQRTPI/4 * (t1 - t2 + t3);
    }else if(l.l == -1){
        t1 = exp_gsl(-kappa*r) + exp_gsl(kappa*r);
        t2 = exp_gsl(-kappa*r)*gsl_sf_erf(-0.5*sqrt(eta)*r + kappa/sqrt(eta));
        t3 = exp_gsl(kappa*r)*gsl_sf_erf(0.5*sqrt(eta)*r + kappa/sqrt(eta));

        res = M_SQRTPI/(2*kappa) * (t1 - t2 - t3)/r;
    }else if(l.l >= 1){
        t1 = (2*l.l - 1)/2. * bar_comp_ew_int(lm {l.l - 1, l.m}, r);
        t2 = -b*bar_comp_ew_int(lm {l.l - 2, l.m}, r);
        t3 = -gsl_pow_int(a, (2*l.l-1))/2 * exp_gsl(-a*a + b/(a*a));
        res = t1 + t2 + t3;
    }
    return res;
}

double Ewald_integral::comp_ewald_int(lm l, double r)
{
    return bar_comp_ew_int(l, r)/gsl_pow_int(r, 2*l.l+1);
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
