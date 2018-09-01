#ifndef AUGMENTED_FUN_H
#define AUGMENTED_FUN_H
#include "utils.h"
#include "log_mesh.h"
#include "atomic_quantity.h"
#include "../../GSL-lib/src/vector.h"

/***************************************************************************//**
* Superclass used to describe the general behaviour of an augmented function.\n
* Contains:\n
* __n__ - Principal quantum number of the augmented function.\n
* __l__ - Combined orbital and angular momentum of the augmented function, l = {l,m}.\n
* __kappa__ - Energy of waves outside augmentation region, used for matching at
* sphere boundary.\n
* __radius__ - Radius of augmentation sphere.\n
* __center__ - Position of the center of the augmentation sphere.\n
* __mesh__ - Logarithmic mesh of the augmentation sphere.\n
* __val__ - values of the augmented function at the mesh points.\n
*******************************************************************************/
class Augmented_function{
public:
    int n;
    lm l;
    double kappa, radius;
    GSL::Vector center;
    Logarithmic_mesh mesh;
    std::vector<double> val;

    //! Obtain value of augmented function at point r.
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

bool operator==(const Augmented_function &a, const Augmented_function &b);

bool operator!=(const Augmented_function &a, const Augmented_function &b);

/***************************************************************************//**
* Class used to represent the augmented Hankel functions.\n
* Contains:\n
* __EH__ - Hankel energy, obtained by solving the radial Schrödinger equaiton.\n
*******************************************************************************/
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


/***************************************************************************//**
* Class used to represent the augmented Bessel functions.\n
* Contains:\n
* __EJ__ - Bessel energy, obtained by solving the radial Schrödinger equaiton.\n
*******************************************************************************/
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
