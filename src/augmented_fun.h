#ifndef AUGMENTED_FUN_H
#define AUGMENTED_FUN_H
#include "utils.h"
#include "log_mesh.h"
#include "atomic_quantity.h"
#include "../../GSL-lib/src/vector.h"

class Augmented_function{
public:
    int n;
    lm l;
    double kappa, radius;
    GSL::Vector center;
    Logarithmic_mesh mesh;
    std::vector<double> val;

    double operator()(const GSL::Vector r);

    Augmented_function();
    Augmented_function(const Augmented_function&);
    Augmented_function(Augmented_function&&);
    Augmented_function(const int n, const lm l, const double kappa,
        const GSL::Vector center, const Logarithmic_mesh mesh);

    virtual ~Augmented_function();
    Augmented_function& operator=(const Augmented_function&);
    Augmented_function& operator=(Augmented_function&&);

    friend bool operator==(const Augmented_function &a, const Augmented_function &b);
    friend bool operator!=(const Augmented_function &a, const Augmented_function &b);

    virtual void update(std::vector<double> v, const double en, const bool core)
     = 0;
};

// Use Simpsons rule to calculate integrals over atomic spheres of augmented functions
double augmented_integral(Augmented_function &a, Augmented_function &b);

// Compare equality of different augmented funcitons
// Two agumented functions are equal if they are centered on the same site and have identical quantum numbers
bool operator==(const Augmented_function &a, const Augmented_function &b);
bool operator!=(const Augmented_function &a, const Augmented_function &b);

class Augmented_Hankel: public Augmented_function{
public:
    double EH;
    Augmented_Hankel();
    Augmented_Hankel(const Augmented_Hankel& a);
    Augmented_Hankel(Augmented_Hankel&& a);
    Augmented_Hankel(const int n, const lm l, const double kappa,
        const GSL::Vector center, const Logarithmic_mesh mesh);

    ~Augmented_Hankel();

    Augmented_Hankel& operator=(const Augmented_Hankel& a);
    Augmented_Hankel& operator=(Augmented_Hankel&& a);

    void update(std::vector<double> v, const double en, const bool core);
};

class Augmented_Bessel: public Augmented_function{
public:
    double EJ;
    Augmented_Bessel();
    Augmented_Bessel(const Augmented_Bessel& a);
    Augmented_Bessel(Augmented_Bessel&& a);
    Augmented_Bessel(const int n, const lm l, const double kappa, const GSL::Vector center, const Logarithmic_mesh mesh);

    ~Augmented_Bessel();

    Augmented_Bessel& operator=(const Augmented_Bessel& a);
    Augmented_Bessel& operator=(Augmented_Bessel&& a);

    void update(std::vector<double> v, const double en, const bool core);
};

namespace std {
	template<>
	struct hash<Augmented_Bessel>{
		size_t operator()(const Augmented_Bessel &f) const
		{
			return (std::hash<int>()((f.l.l >> 1) + f.l.m)^
			std::hash<int>()(f.n >> (f.l.l + f.l.m)));
		}
	};
}

#endif // AUGMENTED_FUN_H
