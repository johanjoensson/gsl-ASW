#include "bloch_sum.h"
#include "ewald_int.h"
#include "envelope_fun.h"

#include <cmath>
#include <string>

Bloch_sum::Bloch_sum(const lm l_n, const double kappa_n, const Crystal& c_n)
 : l(l_n), kappa(kappa_n), c(c_n),
   eta(GSL::exp(GSL::log(6.5) + 2./3*GSL::log(4*M_PI/3.) -
   2./3*GSL::log(c.volume)).val)
{}


GSL::Complex Bloch_sum::calc_d1(const GSL::Vector& tau, const GSL::Vector& kp)
{
    GSL::Complex d1(0,0), e(0,0);
    GSL::Vector kn;
    GSL::Result tmp;
    GSL::Complex c_fac(-1.0, 0.0);

    double k = kp.norm(), dot = GSL::dot(kp, tau);
    // Loop over all K-vectors, including (0, 0,0)
    for(auto Kn : c.Kn_vecs){
        kn = kp + Kn;
        k = kn.norm();
        dot = GSL::dot(kn, tau);
        e = GSL::exp(GSL::Complex(0., dot));
        tmp = GSL::pow_int(k, l.l)*cubic_harmonic(l, kn)*
              GSL::exp((-kappa*kappa - k*k)/eta)/(-kappa*kappa - k*k);
        d1 += e*tmp.val;
    }
    if(l.l % 4 == 1){
        c_fac = GSL::Complex(0.0, 1.0);
    }else if(l.l%4 == 2){
        c_fac = GSL::Complex(1.0, 0.0);
    }else if(l.l%4 == 3){
        c_fac = GSL::Complex(0., -1.0);
    }
    d1 *= 4*M_PI*c_fac/this->c.volume;
    return d1;
}

GSL::Complex Bloch_sum::calc_d1_dot(const GSL::Vector& tau, const GSL::Vector& kp)
{
    GSL::Complex d1_dot(0,0), e(0,0);
    GSL::Vector kn;
    GSL::Result tmp;
    GSL::Complex c_fac(1.0, 0.0);
    double k = kp.norm(), dot = GSL::dot(kp, tau);
    // Loop over all K-vectors, including (0, 0,0)
    for(auto Kn : c.Kn_vecs){
        kn = kp + Kn;
        k = kn.norm();
        dot = GSL::dot(kn, tau);
        e = GSL::exp(GSL::Complex(0., dot));
        tmp = GSL::pow_int(k, l.l)*cubic_harmonic(l, kn)*
              GSL::exp((-kappa*kappa - k*k)/eta)/
              GSL::pow_int(-kappa*kappa - k*k, 2);
        d1_dot += e*tmp.val;
    }
    if(l.l % 4 == 1){
        c_fac = GSL::Complex(0.0, -1.0);
    }else if(l.l%4 == 2){
        c_fac = GSL::Complex(-1.0, 0.0);
    }else if(l.l%4 == 3){
        c_fac = GSL::Complex(0., 1.0);
    }
    d1_dot *= 4*M_PI*c_fac/this->c.volume;
    d1_dot += 1./this->eta * calc_d1(tau, kp);
    return d1_dot;
}

GSL::Complex Bloch_sum::calc_d2(const GSL::Vector& tau, const GSL::Vector& kp)
{
	Ewald_integral I;
	I.set_kappa(this->kappa);
	I.set_ewald_param(this->eta);

    GSL::Complex d2 = GSL::Complex(0., 0.), e = GSL::Complex(0,0);
    double t = 0, dot = 0;
    GSL::Result tmp;
    GSL::Vector tau_mu;
    // Loop over all lattice vectors
    for(auto Rn : c.Rn_vecs){
        // Do not add unit cell contribution at (0, 0, 0)
        if(Rn == GSL::Vector(3) && tau == GSL::Vector(3)){
            continue;
        }
        tau_mu = tau - Rn;
        t = tau_mu.norm();
        dot = GSL::dot(kp, tau_mu);
        tmp = GSL::pow_int(t, l.l)*cubic_harmonic(l, tau_mu)*I.ewald_int(l, t);
        e = GSL::exp(GSL::Complex(0., -dot));
        d2 += e*tmp.val;
    }
    d2 *= GSL::pow_int(2, l.l)*GSL::exp(GSL::Complex(0., GSL::dot(kp,tau)));
    return d2;
}

GSL::Complex Bloch_sum::calc_d2_dot(const GSL::Vector& tau, const GSL::Vector& kp)
{
	Ewald_integral I;
	I.set_kappa(this->kappa);
	I.set_ewald_param(this->eta);

    GSL::Complex d2_dot = GSL::Complex(0., 0.), e = GSL::Complex(0,0);
    double t = 0, dot = 0;
    GSL::Result tmp;
    GSL::Vector tau_mu = tau;
    // Loop over all lattice vectors
    for(auto Rn : c.Rn_vecs){
        // Do not add unit cell contribution at (0, 0, 0)
        if(Rn == GSL::Vector(3) && tau == GSL::Vector(3)){
            continue;
        }
        tau_mu = tau - Rn;
        t = tau_mu.norm();
        dot = GSL::dot(kp, tau_mu);
        tmp = GSL::pow_int(t, l.l)*cubic_harmonic(l, tau_mu)*
              I.ewald_int(lm {l.l - 1, 0}, t);
        e = GSL::exp(GSL::Complex(0., -dot));
        d2_dot += e*tmp.val;
    }
    d2_dot *= 1./4*GSL::pow_int(2, l.l)*
              GSL::exp(GSL::Complex(0., GSL::dot(kp,tau)));
    return d2_dot;
}

GSL::Complex Bloch_sum::calc_d3(const GSL::Vector& tau)
{
    // delta(tau) * delta(l, 0)
    if(tau != GSL::Vector(3) || l.l != 0){
        return GSL::Complex(0., 0.);
    }
    GSL::Result d3 = -kappa/(4*M_SQRTPI)*(GSL::erfc(-kappa/std::sqrt(eta)) -
		     GSL::erfc(kappa/std::sqrt(eta))) - 
		     std::sqrt(eta)/(2*M_PI)*GSL::exp(-kappa*kappa/eta) +
		     kappa/(2*M_SQRTPI);

    return GSL::Complex(d3.val, 0);
}

GSL::Complex Bloch_sum::calc_d3_dot(const GSL::Vector& tau)
{
    if(tau != GSL::Vector(3) || l.l != 0){
        return GSL::Complex(0., 0.);
    }
    // delta(tau) * delta(l, 0)
    GSL::Result d3_dot = 1/(8*M_SQRTPI*kappa) * (GSL::erfc(-kappa/std::sqrt(eta)) - 
        GSL::erfc(kappa/std::sqrt(eta))) - 
	1./(4*M_SQRTPI*kappa);
    return GSL::Complex(d3_dot.val, 0);
}

GSL::Complex Bloch_sum::hankel_envelope(const GSL::Vector& tau, const GSL::Vector& kp)
{
    return calc_d1(tau, kp) + calc_d2(tau, kp) + calc_d3(tau);
}

GSL::Complex Bloch_sum::hankel_envelope_dot(const GSL::Vector& tau, const GSL::Vector& kp)
{
    return calc_d1_dot(tau, kp) + calc_d2_dot(tau, kp) + calc_d3_dot(tau);
}
