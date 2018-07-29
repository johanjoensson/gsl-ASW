#include "bloch_sum.h"
#include "ewald_int.h"

Bloch_sum::Bloch_sum(const lm l, const double kappa, const Crystal c)
 : l(l), kappa(kappa), c(c)
{
	eta = GSL::exp(GSL::log(6.5) + 2./3*GSL::log(4*M_PI/3) -
	2./3*GSL::log(c.volume)).val;
}

GSL::Complex Bloch_sum::calc_d1(GSL::Vector& tau, GSL::Vector& kp)
{
    GSL::Complex d1(0,0), e(0,0);
    GSL::Vector kn;
    GSL::Result tmp;
    GSL::Complex c (-1.0, 0.0);
    double k = 0, dot = 0;
    for(size_t i = 0; i < this->c.Kn_vecs.size(); i++){
        kn = kp + this->c.Kn_vecs[i];
        k = kn.norm();
        dot = GSL::dot(kn, tau);
        e = GSL::exp(GSL::Complex(0., dot));
        tmp = GSL::pow_int(k, l.l)*cubic_harmonic(l, kn)*
        GSL::exp((-kappa*kappa - k*k)/eta)/(kappa*kappa - k*k);
        d1 += e*tmp.val;
    }
    if(l.l % 2 == 0 && l.l/2 % 2 == 1){
        c = GSL::Complex(1.0, 0.0);
    }else if(l.l%2 != 0 && l.l/2 == 0){
        c = GSL::Complex(0., 1.0);
    }else if(l.l%2 != 0 && l.l/2 != 0){
        c = GSL::Complex(0., -1.0);
    }
    d1 *= 4*M_PI/this->c.volume;
    return d1;
}

GSL::Complex Bloch_sum::calc_d2(GSL::Vector& tau, GSL::Vector& kp)
{
	Ewald_integral I;
	I.set_kappa(this->kappa);
	I.set_ewald_param(this->eta);

    GSL::Complex d2 = GSL::Complex(0., 0.), e = GSL::Complex(0,0);
    double t = 0, dot = 0;
    GSL::Result tmp;
    GSL::Vector tau_mu = tau;
    for(size_t i = 0; i < this->c.Rn_vecs.size(); i++){
        tau_mu = tau - this->c.Rn_vecs[i];
        t = tau_mu.norm();
        dot = GSL::dot(kp, tau_mu);
        tmp = GSL::pow_int(t, l.l)*cubic_harmonic(l, tau_mu)*I.ewald_int(l, t);
        e = GSL::exp(GSL::Complex(0., -dot));
        d2 += e*tmp.val;
    }
    d2 *= 2./M_SQRTPI*GSL::pow_int(2, l.l)*GSL::exp(GSL::Complex(0., GSL::dot(kp,tau)));
    return d2;
}

GSL::Complex Bloch_sum::calc_d3(GSL::Vector& tau)
{
    if(tau == GSL::Vector(3) && l.l == 0){
        return GSL::Complex(0., 0.);
    }
    GSL::Result d1 = -kappa/(2*M_SQRTPI)*GSL::erf(-kappa/sqrt(eta)) -
    sqrt(eta)/(2*M_PI)*GSL::exp(-kappa*kappa/eta) - kappa/(2*M_SQRTPI);
    return GSL::Complex(d1.val, 0);
}
