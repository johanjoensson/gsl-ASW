#include "envelope_fun.h"
#include "spherical_fun.h"

Envelope_function::Envelope_function()
 : center(), l(), kappa()
{}

Envelope_function::Envelope_function(const Atom& center, lm l, double kappa)
 : center(center), l(l), kappa(kappa)
{}

/*
Envelope_function::Envelope_function(Envelope_function &a)
 : center(a.center), l(a.l), kappa(a.kappa)
{}

Envelope_function::Envelope_function(Envelope_function &&a)
 : center(), l(), kappa()
{
    std::swap(this->center, a.center);
    std::swap(this->l, a.l);
    std::swap(this->kappa, a.kappa);
}

double Envelope_function::operator()(GSL::Vector r)
{
    double x = (r - center.pos).norm();
    return (this->barred_fun(x)*cubic_harmonic(l,r-center.pos)).val;
}

Envelope_function& Envelope_function::operator=(const Envelope_function &a)
{
    this->center = a.center;
    this->l = a.l;
    this->kappa = a.kappa;
    return *this;
}

Envelope_function& Envelope_function::operator=(Envelope_function &&a)
{
    std::swap(this->center,a.center);
    std::swap(this->l, a.l);
    std::swap(this->kappa, a.kappa);
    return *this;
}
*/
double Envelope_Hankel::barred_fun(double x)
{
    Hankel_function h(l);
    return GSL::pow_int(kappa, l.l+1)*h(kappa*x);
}

double Envelope_Bessel::barred_fun(double x)
{
    Bessel_function j(l);
    return GSL::pow_int(kappa, -l.l)*j(kappa*x);
}

double Envelope_Neumann::barred_fun(double x)
{
    Neumann_function n(l);
    return GSL::pow_int(kappa, l.l+1)*n(kappa*x);
}

double off_atomic_integral(Envelope_Hankel& H1, Envelope_Hankel& H2)
{
    double res = 0.;
    double rs = H1.center.mesh.r.back();
    Envelope_Hankel H1m1, H2m1, H1p1, H2p1;
    H1p1 = Envelope_Hankel(H1.center, lm {H1.l.l + 1, H1.l.m}, H1.kappa);
    H2p1 = Envelope_Hankel(H2.center, lm {H2.l.l + 1, H2.l.m}, H2.kappa);
    H1m1 = Envelope_Hankel(H1.center, lm {H1.l.l - 1, H1.l.m}, H1.kappa);
    H2m1 = Envelope_Hankel(H2.center, lm {H2.l.l - 1, H2.l.m}, H2.kappa);

    if(H1.l.l != H2.l.l){
        return 0.;
    }else if(H1.center != H2.center){
        return 0.;
    }

    if(H1.kappa != H2.kappa){
        res = H1.barred_fun(rs)*H2p1.barred_fun(rs) -
              H1p1.barred_fun(rs)*H2.barred_fun(rs);
        res *= 1./(H1.kappa*H1.kappa - H2.kappa*H2.kappa);
    }else{
        res = H1.barred_fun(rs)*H1.barred_fun(rs) -
              H1m1.barred_fun(rs)*H2p1.barred_fun(rs);
        res *= rs/2;
    }

    res *= rs*rs;
    return res;
}

double atomic_integral(Envelope_Hankel& H1, Envelope_Bessel& J2)
{
    double res = 0.;
    double rs = H1.center.mesh.r.back();
    Envelope_Hankel H1m1, H1p1;
    Envelope_Bessel J2m1, J2p1;
    H1p1 = Envelope_Hankel(H1.center, lm {H1.l.l + 1, H1.l.m}, H1.kappa);
    J2p1 = Envelope_Bessel(J2.center, lm {J2.l.l + 1, J2.l.m}, J2.kappa);
    H1m1 = Envelope_Hankel(H1.center, lm {H1.l.l - 1, H1.l.m}, H1.kappa);
    J2m1 = Envelope_Bessel(J2.center, lm {J2.l.l - 1, J2.l.m}, J2.kappa);

    if(H1.center != J2.center){
        return 0.;
    }

    if(H1.kappa != J2.kappa){
        res = rs*rs*(H1.barred_fun(rs)*J2m1.barred_fun(rs) -
              H1.kappa*H1.kappa*H1m1.barred_fun(rs)*J2.barred_fun(rs)) - 1;
        res *= 1./(H1.kappa*H1.kappa - J2.kappa*J2.kappa);
    }else{
        res = rs*rs*rs*(2*H1.barred_fun(rs)*J2.barred_fun(rs)-
              H1.kappa*H1.kappa*H1m1.barred_fun(rs)*J2p1.barred_fun(rs) -
              1./(H1.kappa*H1.kappa)*H1p1.barred_fun(rs)*J2m1.barred_fun(rs)) +
              (2*H1.l.l + 1)/(H1.kappa*H1.kappa);
        res *= 1./4;
    }

    return res;
}

double atomic_integral(Envelope_Bessel& J1, Envelope_Hankel& H2)
{
    return atomic_integral(H2, J1);
}

double atomic_integral(Envelope_Bessel& J1, Envelope_Bessel& J2)
{
    double res = 0;
    double rs = J1.center.mesh.r.back();
    Envelope_Bessel J1m1, J2m1, J1p1, J2p1;
    J1p1 = Envelope_Bessel(J1.center, lm {J1.l.l + 1, J1.l.m}, J1.kappa);
    J1m1 = Envelope_Bessel(J1.center, lm {J1.l.l - 1, J1.l.m}, J1.kappa);
    J2p1 = Envelope_Bessel(J2.center, lm {J2.l.l + 1, J2.l.m}, J2.kappa);
    J2m1 = Envelope_Bessel(J2.center, lm {J2.l.l - 1, J2.l.m}, J2.kappa);

    if(J1.center != J2.center){
        return 0.;
    }

    if(J1.kappa != J2.kappa){
        res = J1.barred_fun(rs)*J2m1.barred_fun(rs) -
              J1m1.barred_fun(rs)*J2.barred_fun(rs);
        res *= 1./(J1.kappa*J1.kappa - J2.kappa*J2.kappa);
    }else{
        res = J1.barred_fun(rs)*J1.barred_fun(rs) -
              J1m1.barred_fun(rs)*J1p1.barred_fun(rs);
        res *= rs/2;
    }

    res *= rs*rs;
    return res;
}
