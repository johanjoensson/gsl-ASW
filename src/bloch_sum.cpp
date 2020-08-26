#include "bloch_sum.h"
#include "ewald_int.h"
#include "envelope_fun.h"

#include <cmath>
#include <string>

void Bloch_sum::Container::add(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp)
{
    Bloch_sum D;
    values_m.insert({{l, kappa, tau, kp}, D(l, kappa, c, tau, kp)});
    dot_values_m.insert({{l, kappa, tau, kp}, D.dot(l, kappa, c, tau, kp)});
}

GSL::Complex Bloch_sum::Container::get(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp)
{
    auto it = values_m.find({l, kappa, tau, kp});
    if(it != values_m.end()){
        return it->second;
    }
    it = values_m.find({l, kappa, -tau, -kp});
    if(it != values_m.end()){
        return (l.l % 2) == 0 ? it->second : -it->second;
    }
    it = values_m.find({l, kappa, tau, -kp});
    if(it != values_m.end()){
        return it->second.conjugate();
    }
    it = values_m.find({l, kappa, -tau, kp});
    if(it != values_m.end()){
        return (l.l % 2) == 0 ? it->second.conjugate() : -it->second.conjugate();
    }
    Bloch_sum D;
    this->add(l, kappa, c, tau, kp);
    return D(l, kappa, c, tau, kp);
}

GSL::Complex Bloch_sum::Container::get_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp)
{
    auto it = dot_values_m.find({l, kappa, tau, kp});
    if(it != dot_values_m.end()){
        return it->second;
    }
    it = dot_values_m.find({l, kappa, -tau, -kp});
    if(it != dot_values_m.end()){
        return (l.l % 2) == 0 ? it->second : -it->second;
    }
    it = dot_values_m.find({l, kappa, tau, -kp});
    if(it != dot_values_m.end()){
        return it->second.conjugate();
    }
    it = dot_values_m.find({l, kappa, -tau, kp});
    if(it != dot_values_m.end()){
        return (l.l % 2) == 0 ? it->second.conjugate() : -it->second.conjugate();
    }
    Bloch_sum D;
    this->add(l, kappa, c, tau, kp);
    return D.dot(l, kappa, c, tau, kp);
}

void Bloch_sum::Container::recalculate_all(const Crystal_t<3, Atom>& c)
{
    auto old_values = values_m;
    auto old_dot_values = dot_values_m;

    values_m.clear();
    dot_values_m.clear();

    for(const auto& it : old_values){
        auto l = std::get<0>(it.first);
        auto kappa = std::get<1>(it.first);
        auto tau = std::get<2>(it.first);
        auto kp = std::get<3>(it.first);
        this->add(l, kappa, c, tau, kp);
    }
    for(const auto& it : old_dot_values){
        auto l = std::get<0>(it.first);
        auto kappa = std::get<1>(it.first);
        auto tau = std::get<2>(it.first);
        auto kp = std::get<3>(it.first);
        this->add(l, kappa, c, tau, kp);
    }
}

GSL::Complex Bloch_sum::calc_d1(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const
{
    GSL::Complex d1(0,0), e;
    GSL::Result tmp;
    GSL::Complex c_fac(-1.0, 0.0);

    const double eta = calc_eta(c.volume());
    double k, dot;
    // Loop over all K-vectors
    for(const auto& Kn : c.Kn_vecs()){
        if((Kn + kp).norm2() < 1e-16){
            continue;
        }
        dot = tau.dot(Kn + kp);
        k = (Kn + kp).norm();
        e = GSL::exp(GSL::Complex(0., dot));
        tmp = GSL::pow_int(k, l.l)*cubic_harmonic(l, (Kn + kp)/k)*
              GSL::exp((-kappa*kappa - k*k)/eta)/(-kappa*kappa - k*k);
        d1 += e*tmp.val;
    }
    if(l.l % 4 == 1){
        c_fac = GSL::Complex(0.0, 1.0);
    }else if(l.l % 4 == 2){
        c_fac = GSL::Complex(1.0, 0.0);
    }else if(l.l % 4 == 3){
        c_fac = GSL::Complex(0., -1.0);
    }
    d1 *= 4*M_PI*c_fac/c.volume();
    return d1;
}

GSL::Complex Bloch_sum::calc_d1_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const
{
    GSL::Complex d1_dot(0,0), e(0,0);
    GSL::Result tmp;
    GSL::Complex c_fac(1.0, 0.0);
    double k = kp.norm<double>(), dot = kp.dot(tau);
    const double eta = calc_eta(c.volume());
    for(const auto& Kn : c.Kn_vecs()){
        // GSL::Vector kn(Kn + kp);
        if((Kn + kp).norm2() < 1e-16){
            continue;
        }
        // k = kn.norm<double>();
        k = (Kn + kp).norm();
        dot = tau.dot(Kn + kp);
        e = GSL::exp(GSL::Complex(0., dot));
        tmp = GSL::pow_int(k, l.l)*cubic_harmonic(l, (Kn + kp)/k)*
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
    d1_dot *= 4*M_PI*c_fac/c.volume();
    d1_dot += 1./eta * calc_d1(l, kappa, c, tau, kp);
    return d1_dot;
}

GSL::Complex Bloch_sum::calc_d2(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const
{
	Ewald_integral I;
    GSL::Complex d2 = GSL::Complex(0., 0.), e = GSL::Complex(0,0);
    double t = 0, dot = 0;
    GSL::Result tmp;
    // Loop over all lattice vectors
    for(const auto Rn : c.Rn_vecs()){
        if((tau - Rn).norm2() < 1e-16){
            continue;
        }
        t = (tau - Rn).norm();
        dot = kp.dot(tau - Rn);
        tmp = GSL::pow_int(t, l.l)*cubic_harmonic(l, (tau - Rn)/t)*I.ewald_int(kappa, calc_eta(c.volume()), l, t);
        e = GSL::exp(GSL::Complex(0., -dot));
        d2 += e*tmp.val;
    }
    d2 *= GSL::pow_int(2, l.l)*GSL::exp(GSL::Complex(0, kp.dot(tau)));
    return d2;
}

GSL::Complex Bloch_sum::calc_d2_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const
{
	Ewald_integral I;
    GSL::Complex d2_dot = GSL::Complex(0., 0.), e = GSL::Complex(0,0);
    double t = 0, dot = 0;
    GSL::Result tmp;
    // Loop over all lattice vectors
    for(const auto Rn : c.Rn_vecs()){
        // Do not add unit cell contribution at (0, 0, 0)
        if((tau - Rn).norm2() < 1e-16){
            continue;
        }
        t = (tau - Rn).norm<double>();
        dot = kp.dot(tau - Rn);
        tmp = GSL::pow_int(t, l.l)*cubic_harmonic(l, (tau - Rn)/t)*
              I.ewald_int(kappa, calc_eta(c.volume()), lm {l.l - 1, 0}, t);
        e = GSL::exp(GSL::Complex(0., -dot));
        d2_dot += e*tmp.val;
    }
    d2_dot *= GSL::pow_int(2, l.l-2)*GSL::exp(GSL::Complex(0., kp.dot(tau)));
    return d2_dot;
}

GSL::Complex Bloch_sum::calc_d3(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau) const
{
    // delta(tau) * delta(l, 0)
    if(tau != GSL::Vector(3) || l.l != 0){
        return GSL::Complex(0., 0.);
    }
    const double eta = calc_eta(c.volume());
    GSL::Result d3 = -kappa/(4*M_SQRTPI)*(GSL::erfc(-kappa/std::sqrt(eta)) -
		     GSL::erfc(kappa/std::sqrt(eta))) -
		     std::sqrt(eta)/(2*M_PI)*GSL::exp(-kappa*kappa/eta) +
		     kappa/(2*M_SQRTPI);

    return GSL::Complex(d3.val, 0);
}

GSL::Complex Bloch_sum::calc_d3_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau) const
{
    if(tau != GSL::Vector(3) || l.l != 0){
        return GSL::Complex(0., 0.);
    }
    const double eta = calc_eta(c.volume());
    // delta(tau) * delta(l, 0)
    GSL::Result d3_dot = 1./(8*M_SQRTPI*kappa) * (GSL::erfc(-kappa/std::sqrt(eta)) -
        GSL::erfc(kappa/std::sqrt(eta))) -
	       1./(4*M_SQRTPI*kappa);

    return GSL::Complex(d3_dot.val, 0);
}

GSL::Complex Bloch_sum::hankel_envelope(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const
{
    return  calc_d1(l, kappa, c, tau, kp) +
            calc_d2(l, kappa, c, tau, kp) +
            calc_d3(l, kappa, c, tau);
}

GSL::Complex Bloch_sum::hankel_envelope_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const
{
    return  calc_d1_dot(l, kappa, c, tau, kp) +
            calc_d2_dot(l, kappa, c, tau, kp) +
            calc_d3_dot(l, kappa, c, tau);
}
