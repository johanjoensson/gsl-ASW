#ifndef AUGMENTED_FUN_H
#define AUGMENTED_FUN_H
#include "utils.h"
#include "log_mesh.h"
#include "atomic_quantity.h"
#include "GSLpp/vector.h"

class Augmented_function{
protected:
    double En;
public:
    lm l;
    double kappa, radius;
    GSL::Vector center;
    Logarithmic_mesh mesh;
    std::vector<double> val;

    double operator()(const GSL::Vector& r) const;

    Augmented_function(): En(), l(), kappa(), radius(), center(), mesh(), val(){}
    Augmented_function(const Augmented_function&) = default;
    Augmented_function(Augmented_function&&) = default;
    virtual ~Augmented_function() = default;

    Augmented_function(const lm l, const double kappa,
        const GSL::Vector& center, const Logarithmic_mesh& mesh);

    Augmented_function& operator=(const Augmented_function&) = default;
    Augmented_function& operator=(Augmented_function&&) = default;

    friend bool operator==(const Augmented_function &a, const Augmented_function &b);
    friend bool operator!=(const Augmented_function &a, const Augmented_function &b);

    virtual void update(std::vector<double>& v, const double en, const bool core)
     = 0;
};

// Use Simpsons rule to calculate integrals over atomic spheres of augmented functions
double augmented_integral(const Augmented_function &a, const Augmented_function &b);

// Compare equality of different augmented funcitons
// Two agumented functions are equal if they are centered on the same site and have identical quantum numbers
bool operator==(const Augmented_function &a, const Augmented_function &b);
bool operator!=(const Augmented_function &a, const Augmented_function &b);

class Augmented_Hankel: public Augmented_function{
public:
    double& EH(){return En;}
    double EH() const {return En;}
    Augmented_Hankel() = default;
    Augmented_Hankel(const Augmented_Hankel& a) = default;
    Augmented_Hankel(Augmented_Hankel&& a) = default;
    ~Augmented_Hankel() = default;

    Augmented_Hankel(const lm l, const double kappa,
        const GSL::Vector& center, const Logarithmic_mesh& mesh);


    Augmented_Hankel& operator=(const Augmented_Hankel& a) = default;
    Augmented_Hankel& operator=(Augmented_Hankel&& a) = default;

    void update(std::vector<double>& v, const double en, const bool core) override;
};

class Augmented_Bessel: public Augmented_function{
public:
    double& EJ(){return En;}
    double EJ() const {return En;}
    Augmented_Bessel() = default;
    Augmented_Bessel(const Augmented_Bessel& a) = default;
    Augmented_Bessel(Augmented_Bessel&& a) = default;
    ~Augmented_Bessel() = default;

    Augmented_Bessel(const lm l, const double kappa, const GSL::Vector& center, const Logarithmic_mesh& mesh);


    Augmented_Bessel& operator=(const Augmented_Bessel& a) = default;
    Augmented_Bessel& operator=(Augmented_Bessel&& a) = default;

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

			res ^= hlm(f.l) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= hkappa(f.kappa) + 0x9e3779b9 + (res<< 6) + (res>> 2);

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

			res ^= hlm(f.l) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= hkappa(f.kappa) + 0x9e3779b9 + (res<< 6) + (res>> 2);

            return res;
		}
	};
}

class Hankel_container
{
private:
    std::unordered_set<Augmented_Hankel> functions;
private:
    void add_function(const Augmented_Hankel&);
    Augmented_Hankel& get_function(const lm& l, const double& kappa,
        const spin& s);

    std::unordered_set<Augmented_Hankel>::iterator begin()
    {
        return functions.begin();
    }
    std::unordered_set<Augmented_Hankel>::iterator end()
    {
        return functions.end();
    }

    std::unordered_set<Augmented_Hankel>::const_iterator cbegin()
    {
        return functions.cbegin();
    }
    std::unordered_set<Augmented_Hankel>::const_iterator cend()
    {
        return functions.cend();
    }
};

class Bessel_container
{
private:
    std::unordered_set<Augmented_Bessel> functions;
private:
    void add_function(const Augmented_Bessel&);
    Augmented_Hankel& get_function(const lm& l, const double& kappa,
        const spin& s);

    std::unordered_set<Augmented_Bessel>::iterator begin()
    {
        return functions.begin();
    }
    std::unordered_set<Augmented_Bessel>::iterator end()
    {
        return functions.end();
    }

    std::unordered_set<Augmented_Bessel>::const_iterator cbegin()
    {
        return functions.cbegin();
    }
    std::unordered_set<Augmented_Bessel>::const_iterator cend()
    {
        return functions.cend();
    }
};

#endif // AUGMENTED_FUN_H
