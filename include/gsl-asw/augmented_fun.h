#ifndef AUGMENTED_FUN_H
#define AUGMENTED_FUN_H
#include "utils.h"
#include "log_mesh.h"
#include "atomic_quantity.h"
#include <numerical-mesh/numerical-mesh.h>
#include "GSLpp/vector.h"

class Augmented_function{
protected:
    double En_m;
public:
    lm l_m;
    double kappa_m;
    Spin s_m;
    Exponential_mesh<1, double> mesh_m;

    double S_m;

    // std::vector<double> val_m;

    double operator()(const GSL::Vector& r) const;

    Augmented_function() = default;
    Augmented_function(const Augmented_function&) = default;
    Augmented_function(Augmented_function&&) = default;
    virtual ~Augmented_function(){}

    Augmented_function(const lm l, const double kappa, const Spin s,
        const Exponential_mesh<1, double>& mesh);

    Augmented_function& operator=(const Augmented_function&) = default;
    Augmented_function& operator=(Augmented_function&&) = default;

    friend bool operator==(const Augmented_function &a, const Augmented_function &b);
    friend bool operator!=(const Augmented_function &a, const Augmented_function &b);

    virtual void update(std::vector<double>& v, const double en, const bool core)
     = 0;

    lm l() const {return l_m;}
    lm& l() {return l_m;}

    double kappa() const {return kappa_m;}
    double& kappa() {return kappa_m;}

    Spin s() const {return s_m;}
    Spin& s() {return s_m;}

    Exponential_mesh<1, double> mesh() const {return mesh_m;}
    Exponential_mesh<1, double>& mesh() {return mesh_m;}

    double S() const {return S_m;}
    // double& S() {return S_m;}

    // std::vector<double> val() const {return val_m;}
    // std::vector<double>& val() {return val_m;}

    double En() const {return En_m;}
    double& En() {return En_m;}
};

// double augmented_integral(const Augmented_function &a, const Augmented_function &b);

// Compare equality of different augmented funcitons
// Two agumented functions are equal if they are centered on the same site and have identical quantum numbers
bool operator==(const Augmented_function &a, const Augmented_function &b);
bool operator!=(const Augmented_function &a, const Augmented_function &b);

class Augmented_Hankel: public Augmented_function{
public:
    double& EH(){return En_m;}
    double EH() const {return En_m;}

    double SH() const {return S_m;}

    Augmented_Hankel() = default;
    Augmented_Hankel(const Augmented_Hankel&) = default;
    Augmented_Hankel(Augmented_Hankel&&) = default;
    ~Augmented_Hankel(){}

    Augmented_Hankel(const lm l, const double kappa, const Spin s,
        const Exponential_mesh<1, double>& mesh);

    Augmented_Hankel& operator=(const Augmented_Hankel& a) = default;
    Augmented_Hankel& operator=(Augmented_Hankel&& a) = default;

    void update(std::vector<double>& v, const double en, const bool core) override;
};

class Augmented_Bessel: public Augmented_function{
public:
    double& EJ(){return En_m;}
    double EJ() const {return En_m;}
    double SJ() const {return S_m;}
    Augmented_Bessel() = default;
    Augmented_Bessel(const Augmented_Bessel&) = default;
    Augmented_Bessel(Augmented_Bessel&&) = default;
    ~Augmented_Bessel(){}

    Augmented_Bessel(const lm l, const double kappa, const Spin s,
        const Exponential_mesh<1, double>& mesh);


    Augmented_Bessel& operator=(const Augmented_Bessel& a) = default;
    Augmented_Bessel& operator=( Augmented_Bessel&& a) = default;

    void update(std::vector<double>& v, const double en, const bool core) override;
};

namespace std {
	template<>
	struct hash<Augmented_Hankel>{
		size_t operator()(const Augmented_Hankel &f) const
		{
            size_t res = 0;
            std::hash<lm> hlm;
            std::hash<double> hkappa;

			res ^= hlm(f.l()) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= hkappa(f.kappa()) + 0x9e3779b9 + (res<< 6) + (res>> 2);

            return res;
		}
	};
	template<>
	struct hash<Augmented_Bessel>{
		size_t operator()(const Augmented_Bessel &f) const
		{
            size_t res = 0;
            std::hash<lm> hlm;
            std::hash<double> hkappa;

			res ^= hlm(f.l()) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= hkappa(f.kappa()) + 0x9e3779b9 + (res<< 6) + (res>> 2);

            return res;
		}
	};
}

class Hankel_container
{
private:
    std::vector<Augmented_Hankel> functions;
public:
    Hankel_container():functions(){};
    void add_function(const Augmented_Hankel&);
    Augmented_Hankel& get_function(const lm& l, const double& kappa,
        const Spin& s);
    size_t get_index(const lm& l, const double& kappa,
        const Spin& s) const;
        Augmented_Hankel& operator()(const lm& l, const double& kappa, const Spin& s)
        {return this->get_function(l, kappa, s);}

    std::vector<Augmented_Hankel>::iterator begin()
    {
        return functions.begin();
    }
    std::vector<Augmented_Hankel>::iterator end()
    {
        return functions.end();
    }

    std::vector<Augmented_Hankel>::const_iterator cbegin()
    {
        return functions.cbegin();
    }
    std::vector<Augmented_Hankel>::const_iterator cend()
    {
        return functions.cend();
    }
    Augmented_Hankel& front()
    {
        return functions.front();
    }
    Augmented_Hankel& back()
    {
        return functions.back();
    }
    const Augmented_Hankel& front() const
    {
        return functions.front();
    }
    const Augmented_Hankel& back() const
    {
        return functions.back();
    }

    size_t size() const
    {return functions.size();}
};

class Bessel_container
{
private:
    std::vector<Augmented_Bessel> functions;
public:
    Bessel_container():functions(){};
    void add_function(const Augmented_Bessel&);
    Augmented_Bessel& get_function(const lm& l, const double& kappa,
        const Spin& s);
    Augmented_Bessel& get_function(size_t index){ return functions[index];}
    size_t get_index(const lm& l, const double& kappa,
        const Spin& s) const;

    lm max_lm() const;
    lm min_lm() const;

    std::vector<Augmented_Bessel>::iterator begin()
    {
        return functions.begin();
    }
    std::vector<Augmented_Bessel>::iterator end()
    {
        return functions.end();
    }
    std::vector<Augmented_Bessel>::const_iterator begin() const
    {
        return functions.begin();
    }
    std::vector<Augmented_Bessel>::const_iterator end() const
    {
        return functions.end();
    }
    std::vector<Augmented_Bessel>::const_iterator cbegin()
    {
        return functions.cbegin();
    }
    std::vector<Augmented_Bessel>::const_iterator cend()
    {
        return functions.cend();
    }
    Augmented_Bessel& front()
    {
        return functions.front();
    }
    Augmented_Bessel& back()
    {
        return functions.back();
    }
    const Augmented_Bessel& front() const
    {
        return functions.front();
    }
    const Augmented_Bessel& back() const
    {
        return functions.back();
    }
    size_t size() const
    {return functions.size();}
};

// Use Simpsons rule to calculate integrals over atomic spheres of augmented functions
double augmented_integral(const Augmented_Hankel &a, const Augmented_Hankel &b);
double augmented_integral(const Augmented_Bessel &a, const Augmented_Bessel &b);
double augmented_integral(const Augmented_Hankel &a, const Augmented_Bessel &b);
double augmented_integral(const Augmented_Bessel &a, const Augmented_Hankel &b);

#endif // AUGMENTED_FUN_H
