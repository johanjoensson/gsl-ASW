#ifndef BLOCH_SUM_H
#define BLOCH_SUM_H
#include "spherical_fun.h"
#include "crystal.h"
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/complex.h"


class Bloch_sum{

    lm l;
    double kappa;
    Crystal c;
    double eta;
public:
    Bloch_sum(const lm l, const double kappa, const Crystal c);

    GSL::Complex calc_d1(GSL::Vector& tau, GSL::Vector& kp);
    GSL::Complex calc_d2(GSL::Vector& tau, GSL::Vector& kp);
    GSL::Complex calc_d3(GSL::Vector& tau);
};


#endif // BLOCH_SUM_H
