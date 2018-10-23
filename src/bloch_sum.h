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

    GSL::Complex calc_d1(const GSL::Vector& tau, const GSL::Vector& kp);
    GSL::Complex calc_d2(const GSL::Vector& tau, const GSL::Vector& kp);
    GSL::Complex calc_d3(const GSL::Vector& tau);
    GSL::Complex calc_d1_dot(const GSL::Vector& tau, const GSL::Vector& kp);
    GSL::Complex calc_d2_dot(const GSL::Vector& tau, const GSL::Vector& kp);
    GSL::Complex calc_d3_dot(const GSL::Vector& tau);
public:
    Bloch_sum(const lm l, const double kappa, const Crystal& c);

    GSL::Complex hankel_envelope(const GSL::Vector& tau, const GSL::Vector& kp);
    GSL::Complex hankel_envelope_dot(const GSL::Vector& tau, const GSL::Vector& kp);
};


#endif // BLOCH_SUM_H
