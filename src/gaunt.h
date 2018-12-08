#ifndef GAUNT_H
#define GAUNT_H
#include "spherical_fun.h"
#include "../../GSL-lib/src/special_functions.h"

GSL::Result gaunt(lm l1, lm l2, lm l3);
GSL::Result real_gaunt(lm l1, lm l2, lm l3);


#endif //GAUNT_H
