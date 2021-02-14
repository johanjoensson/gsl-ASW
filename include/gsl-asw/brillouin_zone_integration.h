
#ifndef BRILLOUIN_ZONE_INTEGRATION_H
#define BRILLOUIN_ZONE_INTEGRATION_H

#include "GSLpp/vector.h"
#include "GSLpp/matrix.h"
#include "GSLpp/special_functions.h"
#include "GSLpp/complex.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>

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
    virtual GSL::Complex measure_observable(const double E, const GSL::Matrix A) const
    {
        return 0*E*A[0][0];
    }
};

class Brillouin_sampler : public Brillouin_zone_integrator
{
private:
    GSL::Matrix energies_m;
protected:
    inline virtual double delta(double x) const = 0;
    inline virtual double theta(double x) const = 0;

    double real_weight(const double enk, const double E)
     const
    {
        size_t n_k = energies_m.dim().first;

        return 1./static_cast<double>(n_k) * 1./(E - enk);
    }

    double imaginary_weight(const double enk, const double E)
     const
    {
        size_t n_k = energies_m.dim().first;

        return 1./static_cast<double>(n_k) * delta(E - enk);
    }

public:
    Brillouin_sampler() = delete;
    Brillouin_sampler(const Brillouin_sampler&) = default;
    Brillouin_sampler(Brillouin_sampler&&) = default;

    Brillouin_sampler(const GSL::Matrix& energies)
     : Brillouin_zone_integrator(), energies_m{energies}
    {}
    virtual ~Brillouin_sampler() = default;

    double measure_real(const double E, const GSL::Matrix& A) const
    {
        double res = 0;
        if(A.dim().first != energies_m.dim().first){
            throw(std::runtime_error("Error in measure_real. Matrix A has been "
            "measured at a dfferents number of k-points than the Kohn-Sham energies."));
        }
        if(A.dim().second != energies_m.dim().second){
            throw(std::runtime_error("Error in measure_real. Matrix A has been "
            "measured at a dfferents number of bands than the Kohn-Sham energies."));
        }

        for(size_t row = 0; row < energies_m.dim().first; row++){
            auto zipped = std::vector<std::pair<double, double>>(energies_m[row].size());
            std::transform(energies_m[row].begin(), energies_m[row].end(), A[row].begin(), zipped.begin(),
                [](const double enk, const double Ank)
                {
                    return std::pair<double, double>{enk, Ank};
                });
            res = std::accumulate(zipped.begin(), zipped.end(), res,
                [=](const double acc, const std::pair<double, double> eA)
                {
                    double enk = eA.first, Ank = eA.second;
                    return acc + Ank*real_weight(E, enk);
                });
        }
        return res;
    }
    double measure_imaginary(const double E, const GSL::Matrix& A) const
    {
        double res = 0;
        if(A.dim().first != energies_m.dim().first){
            throw(std::runtime_error("Error in measure_imaginary. Matrix A has been "
            "measured at a dfferents number of k-points than the Kohn-Sham energies."));
        }
        if(A.dim().second != energies_m.dim().second){
            throw(std::runtime_error("Error in measure_imaginary. Matrix A has been "
            "measured at a dfferents number of bands than the Kohn-Sham energies."));
        }

        for(size_t row = 0; row < energies_m.dim().first; row++){
            auto zipped = std::vector<std::pair<double, double>>(energies_m[row].size());
            std::transform(energies_m[row].begin(), energies_m[row].end(), A[row].begin(), zipped.begin(),
                [](const double enk, const double Ank)
                {
                    return std::pair<double, double>{enk, Ank};
                });
            res = std::accumulate(zipped.begin(), zipped.end(), res,
                [=](const double acc, const std::pair<double, double> eA)
                {
                    double enk = eA.first, Ank = eA.second;
                    return acc + Ank*imaginary_weight(E, enk);
                });
        }
        return res;
    }

    double dos(const double E) const override
    {
        GSL::Matrix A{energies_m.dim().first, energies_m.dim().second};
        for(auto row : A){
            for(auto& elem : row){
                elem = 1.;
            }
        }
        return measure_imaginary(E, A);
    }

    GSL::Complex measure_observable(const double E, const GSL::Matrix A) const override
    {
        return {measure_real(E, A), -M_PI*measure_imaginary(E, A) };
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
