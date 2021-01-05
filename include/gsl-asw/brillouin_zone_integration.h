
#ifndef BRILLOUIN_ZONE_INTEGRATION_H
#define BRILLOUIN_ZONE_INTEGRATION_H

#include "GSLpp/vector.h"
#include "GSLpp/matrix.h"
#include "GSLpp/special_functions.h"
#include <vector>
#include <algorithm>
#include <numeric>

class Brillouin_zone_integrator
{
protected:
    Brillouin_zone_integrator() = default;
public:
    Brillouin_zone_integrator(const Brillouin_zone_integrator&) = default;
    Brillouin_zone_integrator(Brillouin_zone_integrator&&) = default;

    virtual ~Brillouin_zone_integrator() = default;

    Brillouin_zone_integrator& operator=(const Brillouin_zone_integrator&) = default;
    Brillouin_zone_integrator& operator=(Brillouin_zone_integrator&&) = default;

    virtual double dos(const double E) const {return 0*E;}
};

class Brillouin_sampler : public Brillouin_zone_integrator
{
private:
    GSL::Matrix energies_m;
protected:
    inline virtual double delta(double x) const = 0;
    inline virtual double theta(double x) const = 0;
public:
    Brillouin_sampler() = delete;
    Brillouin_sampler(const Brillouin_sampler&) = default;
    Brillouin_sampler(Brillouin_sampler&&) = default;

    Brillouin_sampler(const GSL::Matrix& energies)
     : Brillouin_zone_integrator(), energies_m{energies}
    {}
    virtual ~Brillouin_sampler() = default;

    double dos(const double E) const override
    {
        size_t n_k = energies_m.dim().first;
        double res = 0;
        for (auto& row : energies_m){
            res = std::accumulate(row.begin(), row.end(), res,
            [=](const double acc, const double enk)
            {
                return acc + delta(E - enk);
            });
        }
        return res/static_cast<double>(n_k);
    }
};

class Simple_sampler : public Brillouin_sampler
{
private:
    double tol_m, h_m;
    inline double delta(const double x) const
    {
        return std::abs(x) < tol_m ? 1/h_m : 0;
    }
    inline double theta(const double x) const
    {
        return x <= 0? 1 : 0;
    }

public:
    Simple_sampler(const GSL::Matrix& energies, const double tol = 1e-12, const double h = 1)
     : Brillouin_sampler(energies), tol_m(tol), h_m(h)
    {}
    ~Simple_sampler() = default;
};

class Fermi_sampler : public Brillouin_sampler
{
private:
    double sigma_m;
    inline double delta(const double x) const
    {
        return 1./(2*sigma_m)*1./(cosh(x/sigma_m) + 1);
    }
    inline double theta(const double x) const
    {
        return 1./(exp(x/sigma_m) + 1);
    }

public:
    Fermi_sampler(const GSL::Matrix& energies, const double sigma)
     : Brillouin_sampler(energies), sigma_m(sigma)
    {}
    ~Fermi_sampler() = default;
};

class Gaussian_sampler : public Brillouin_sampler
{
private:
    double sigma_m;
    inline double delta(const double x) const
    {
        return 1./(sigma_m*M_SQRT2*M_SQRTPI)*GSL::exp(-0.5*GSL::pow_int(x/sigma_m, 2)).val;
    }
    inline double theta(const double x) const
    {
        return 1./2 * (1 - GSL::erf(x/sigma_m).val);
    }

public:
    Gaussian_sampler(const GSL::Matrix& energies, const double sigma)
     : Brillouin_sampler(energies), sigma_m(sigma)
    {}
    ~Gaussian_sampler() = default;
};

class Lorentzian_sampler : public Brillouin_sampler
{
private:
    double sigma_m;
    inline double delta(const double x) const
    {
        return 1./(sigma_m*M_PI)*1./(1 + GSL::pow_int(x/sigma_m, 2));
    }
    inline double theta(const double x) const
    {
        return 1./2 - atan(x/sigma_m)*M_1_PI;
    }

public:
    Lorentzian_sampler(const GSL::Matrix& energies, const double sigma)
     : Brillouin_sampler(energies), sigma_m(sigma)
    {}
    ~Lorentzian_sampler() = default;
};



class MethfesselPaxton_sampler : public Brillouin_sampler
{
private:
    double sigma_m;
    unsigned int order_m;
    inline double alpha(const unsigned int n) const
    {
        double sign = n % 2 == 0 ? 1 : -1;
        return sign/(GSL::fact(n).val*GSL::pow_uint(4, n)*M_SQRTPI);
    }
    inline double delta(const double x) const
    {
        std::vector<unsigned int> i(2*order_m + 1);
        std::iota(i.begin(), i.end(), 0);
        std::vector<double> alphas(2*order_m + 1);
        std::transform(i.begin(), i.end(), alphas.begin(),
            [this](const unsigned int n){ return n % 2 == 0 ? alpha(n/2) : 0;}
        );
        return 1./sigma_m*(GSL::hermite_phys_series(static_cast<int>(2*order_m + 1), x/sigma_m, alphas)*GSL::exp(-GSL::pow_int(x/sigma_m, 2))).val;
    }
    inline double theta(const double x) const
    {
        std::vector<unsigned int> i(2*order_m + 1);
        std::iota(i.begin(), i.end(), 0);
        std::vector<double> alphas(2*order_m + 1);
        std::transform(i.begin(), i.end(), alphas.begin(),
            [this](const unsigned int n){ return n % 2 == 0 ? 0  : alpha((n - 1)/2);}
        );
        return (GSL::hermite_phys_series(static_cast<int>(2*order_m + 1), x/sigma_m, alphas)*GSL::exp(-GSL::pow_int(x/sigma_m, 2))).val;
    }

public:
    MethfesselPaxton_sampler(const GSL::Matrix& energies, const double sigma, const unsigned int order = 1)
     : Brillouin_sampler(energies), sigma_m(sigma), order_m(order)
    {}
    ~MethfesselPaxton_sampler() = default;
};

class ColdSmearing_sampler : public Brillouin_sampler
{
private:
    double sigma_m;
    inline double delta(const double x) const
    {
        return 1./sigma_m*1./M_SQRTPI*GSL::exp(-GSL::pow_int(x/sigma_m - M_SQRT1_2, 2)).val*(2 - M_SQRT2*x/sigma_m);
    }
    inline double theta(const double x) const
    {
        return 0.5*(M_SQRT2/M_SQRTPI*GSL::exp(-GSL::pow_int(x/sigma_m, 2) - M_SQRT2*x/sigma_m - 0.5)
            + 1 - GSL::erf(x/sigma_m + M_SQRT1_2)).val;
    }

public:
    ColdSmearing_sampler(const GSL::Matrix& energies, const double sigma)
     : Brillouin_sampler(energies), sigma_m(sigma)
    {}
    ~ColdSmearing_sampler() = default;
};



#endif // BRILLOUIN_ZONE_INTEGRATION_H
