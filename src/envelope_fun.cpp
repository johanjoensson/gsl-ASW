#include "envelope_fun.h"
#include "spherical_fun.h"

double Envelope_Hankel::barred_fun(const double x) const
{
    if(l_m.l < 0){
        Envelope_Hankel h(center_m, lm {-l_m.l - 1, l_m.m}, kappa_m);
        return GSL::pow_int(kappa_m, 2*l_m.l + 1)*h.barred_fun(x);
    }
    if(std::abs(kappa_m) < 1e-16){
        if(l_m.l == 0){
            return 1./x;
        }
        return (GSL::doublefact(static_cast<unsigned int>(2*l_m.l - 1)).val/GSL::pow_int(x, l_m.l + 1));
    }
    Hankel_function h(l_m);
    return GSL::pow_int(kappa_m, l_m.l + 1)*h(kappa_m*x);
}

double Envelope_Bessel::barred_fun(const double x) const
{
    if(l_m.l < 0){
        Envelope_Bessel jp1(center_m, lm {l_m.l  + 1, l_m.m}, kappa_m), jp2(center_m, lm {l_m.l + 2, l_m.m}, kappa_m);
        return (2*l_m.l + 3)*jp1.barred_fun(x)/x + kappa_m*kappa_m*jp2.barred_fun(x);
    }
    if(std::abs(kappa_m) < 1e-16){
        return GSL::pow_int(x, l_m.l)/GSL::doublefact(static_cast<unsigned int>(2*l_m.l + 1)).val;
    }
    Bessel_function j(l_m);
    return GSL::pow_int(kappa_m, -l_m.l)*j(kappa_m*x);
}

double Envelope_Neumann::barred_fun(const double x) const
{
    if(l_m.l < 0){
        Envelope_Neumann np1(center_m, lm {l_m.l + 1, l_m.m}, kappa_m), np2(center_m, lm {l_m.l + 2, l_m.m}, kappa_m);
        return ((2*l_m.l + 3)/x*np1.barred_fun(x) + np2.barred_fun(x))/(kappa_m*kappa_m);
    }
    if(std::abs(kappa_m) < 1e-16){
        if(l_m.l == 0){
            return 1./x;
        }
        return (GSL::doublefact(static_cast<unsigned int>(2*l_m.l - 1)).val/GSL::pow_int(x, l_m.l + 1));
    }
    Neumann_function n(l_m);
    return -GSL::pow_int(kappa_m, l_m.l+1)*n(kappa_m*x);
}

double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Hankel& H2, const double rs)
{
    if(H1.l_m != H2.l_m){
        return 0.;
    }else if(H1.center_m != H2.center_m){
        return 0.;
    }

    double res = 0.;
    Envelope_Hankel H1p1 {H1.center_m, lm {H1.l_m.l + 1, H1.l_m.m}, H1.kappa_m};
    Envelope_Hankel H2p1 {H2.center_m, lm {H2.l_m.l + 1, H2.l_m.m}, H2.kappa_m};
    Envelope_Hankel H1m1 {H1.center_m, lm {H1.l_m.l - 1, H1.l_m.m}, H1.kappa_m};
    Envelope_Hankel H2m1 {H2.center_m, lm {H2.l_m.l - 1, H2.l_m.m}, H2.kappa_m};


    if(H1.kappa_m != H2.kappa_m){
        res = H1.barred_fun(rs)*H2p1.barred_fun(rs) -
              H1p1.barred_fun(rs)*H2.barred_fun(rs);
        res *= 1./(-H1.kappa_m*H1.kappa_m + H2.kappa_m*H2.kappa_m);
    }else{
    	res = 0.5*rs*( H1m1.barred_fun(rs)*H2p1.barred_fun(rs) -
	      H1.barred_fun(rs)*H2.barred_fun(rs) );
    }

    res *= rs*rs;
    return res;
}

double atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs)
{
    double res = 0.;
    Envelope_Hankel H1p1 {H1.center_m, lm {H1.l_m.l + 1, H1.l_m.m}, H1.kappa_m};
    Envelope_Bessel J2p1 {J2.center_m, lm {J2.l_m.l + 1, J2.l_m.m}, J2.kappa_m};
    Envelope_Hankel H1m1 {H1.center_m, lm {H1.l_m.l - 1, H1.l_m.m}, H1.kappa_m};
    Envelope_Bessel J2m1 {J2.center_m, lm {J2.l_m.l - 1, J2.l_m.m}, J2.kappa_m};

    if(H1.center_m != J2.center_m){
        return 0.;
    }

    if(H1.kappa_m != J2.kappa_m){
        res = rs*rs*(H1.barred_fun(rs)*J2m1.barred_fun(rs) +
               H1.kappa_m*H1.kappa_m*H1m1.barred_fun(rs)*J2.barred_fun(rs)) - 1;
        res *= 1./(-H1.kappa_m*H1.kappa_m + J2.kappa_m*J2.kappa_m);
    }else{
        res = rs*rs*rs*(2*H1.barred_fun(rs)*J2.barred_fun(rs)  +
              H1.kappa_m*H1.kappa_m*H1m1.barred_fun(rs)*J2p1.barred_fun(rs) +
              1./(H1.kappa_m*H1.kappa_m)*H1p1.barred_fun(rs)*J2m1.barred_fun(rs)) -
              (2*H1.l_m.l + 1)/(H1.kappa_m*H1.kappa_m);
        res *= 1./4;
    }

    return res;
}

double atomic_integral(const Envelope_Bessel& J1, const Envelope_Hankel& H2, const double rs)
{
    return atomic_integral(H2, J1, rs);
}

double atomic_integral(const Envelope_Bessel& J1, const Envelope_Bessel& J2, const double rs)
{
    double res = 0;
    Envelope_Bessel J1p1 {J1.center_m, lm {J1.l_m.l + 1, J1.l_m.m}, J1.kappa_m};
    Envelope_Bessel J1m1 {J1.center_m, lm {J1.l_m.l - 1, J1.l_m.m}, J1.kappa_m};
    Envelope_Bessel J2p1 {J2.center_m, lm {J2.l_m.l + 1, J2.l_m.m}, J2.kappa_m};
    Envelope_Bessel J2m1 {J2.center_m, lm {J2.l_m.l - 1, J2.l_m.m}, J2.kappa_m};

    if(J1.center_m != J2.center_m){
        return 0.;
    }

    if(J1.kappa_m != J2.kappa_m){
        res = J1.barred_fun(rs)*J2m1.barred_fun(rs) -
              J1m1.barred_fun(rs)*J2.barred_fun(rs);
        res *= 1./(-J1.kappa_m*J1.kappa_m + J2.kappa_m*J2.kappa_m);
    }else{
        res = 0.5*rs*(J1.barred_fun(rs)*J2.barred_fun(rs) -
              J1m1.barred_fun(rs)*J2p1.barred_fun(rs) );
    }

    res *= rs*rs;
    return res;
}
