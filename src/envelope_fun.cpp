#include <envelope_fun.h>
#include <spherical_fun.h>

double Envelope_Hankel::barred_fun(const double x) const
{
    if(l_m.l < 0){
        Envelope_Hankel h(lm {-l_m.l - 1, l_m.m}, kappa_m);
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
        Envelope_Bessel jp1(lm {l_m.l  + 1, l_m.m}, kappa_m), jp2(lm {l_m.l + 2, l_m.m}, kappa_m);
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
        Envelope_Neumann np1(lm {l_m.l + 1, l_m.m}, kappa_m), np2(lm {l_m.l + 2, l_m.m}, kappa_m);
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
    }

    double k1 = H1.kappa(), k2 = H2.kappa();
    if(k1 != k2){
        return -1./(k2*k2 - k1*k1) * wronskian(H1, H2, rs);
    }else{
        Envelope_Hankel Hp1 {lm {H1.l_m.l + 1, H1.l_m.m}, H1.kappa_m};
        Envelope_Hankel Hm1 {lm {H1.l_m.l - 1, H1.l_m.m}, H1.kappa_m};
    	return -rs*rs*rs/2  *  (H1.barred_fun(rs)*H2.barred_fun(rs) -
                                Hp1.barred_fun(rs)*Hp1.barred_fun(rs));
    }
}

double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs)
{
    if(H1.l_m != J2.l_m){
        return 0.;
    }

    double k1 = H1.kappa(), k2 = J2.kappa();
    if(k1 != k2){
        return -wronskian(H1, J2, rs)/(k2*k2 - k1*k1);
    }else{
    	return 0.;
    }
}

double atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs)
{
    double k1 = H1.kappa(), k2 = J2.kappa();
    if(k1 != k2){
        return (wronskian(H1, J2, rs) - 1)/(k2*k2 - k1*k1);
    }else{
        Envelope_Hankel Hm1 {lm {H1.l_m.l - 1, H1.l_m.m}, k1};
        Envelope_Hankel Hp1 {lm {H1.l_m.l - 1, H1.l_m.m}, k1};
        Envelope_Bessel Jm1 {lm {J2.l_m.l - 1, J2.l_m.m}, k2};
        Envelope_Bessel Jp1 {lm {J2.l_m.l - 1, J2.l_m.m}, k2};

        double res =    rs*rs*rs*(2*H1.barred_fun(rs)*J2.barred_fun(rs)  +
                        k1*k1*Hm1.barred_fun(rs)*Jp1.barred_fun(rs));
	    if(std::abs(H1.kappa_m) > 1e-10){
            res += rs*rs*rs * 1./(H1.kappa_m*H1.kappa_m)*Hp1.barred_fun(rs)*Jm1.barred_fun(rs)
            -(2*H1.l_m.l + 1)/(k1*k1);
	    }
        res *= 1./4;
        return res;
    }
}

double atomic_integral(const Envelope_Bessel& J1, const Envelope_Hankel& H2, const double rs)
{
    return atomic_integral(H2, J1, rs);
}

double atomic_integral(const Envelope_Bessel& J1, const Envelope_Bessel& J2, const double rs)
{
    double k1 = J1.kappa(), k2 = J2.kappa();
    if(k1 != k2){
        return wronskian(J1, J2, rs)/(k2*k2 - k1*k1);
    }else{
        Envelope_Bessel Jm1 {lm {J1.l().l - 1, J1.l().m}, k1};
        Envelope_Bessel Jp1 {lm {J1.l().l + 1, J1.l().m}, k1};

        return  rs*rs*rs/2*(J1.barred_fun(rs)*J1.barred_fun(rs) -
                Jm1.barred_fun(rs)*Jp1.barred_fun(rs) );
    }
}

//Wronskian
double wronskian(const Envelope_Hankel& f, const Envelope_Hankel& g, const double rs)
{
	if(std::abs(f.kappa() - g.kappa()) < 1e-10 && f.l().l == g.l().l){
		return 0;
	}
	const Envelope_Hankel fp1(lm{f.l().l + 1, f.l().m}, f.kappa());
	const Envelope_Hankel gp1(lm{g.l().l + 1, g.l().m}, g.kappa());
	return rs*rs*(fp1.barred_fun(rs)*g.barred_fun(rs) - gp1.barred_fun(rs)*f.barred_fun(rs));
}

double wronskian(const Envelope_Hankel& f, const Envelope_Bessel& g, const double rs)
{
	if(std::abs(f.kappa() - g.kappa()) < 1e-10 && f.l().l == g.l().l){
		return 1;
	}
    double k2 = g.kappa();
	const Envelope_Hankel fp1(lm{f.l().l + 1, f.l().m}, f.kappa());
	const Envelope_Hankel gp1(lm{g.l().l + 1, g.l().m}, g.kappa());
	return rs*rs*(fp1.barred_fun(rs)*g.barred_fun(rs) + k2*k2*gp1.barred_fun(rs)*f.barred_fun(rs));
}

double wronskian(const Envelope_Bessel& f, const Envelope_Bessel& g, const double rs)
{
	if(std::abs(f.kappa() - g.kappa()) < 1e-10 && f.l().l == g.l().l){
		return 0;
	}
	const Envelope_Hankel fm1(lm{f.l().l - 1, f.l().m}, f.kappa());
	const Envelope_Hankel gm1(lm{g.l().l - 1, g.l().m}, g.kappa());
	return rs*rs*(f.barred_fun(rs)*gm1.barred_fun(rs) - fm1.barred_fun(rs)*g.barred_fun(rs));
}
