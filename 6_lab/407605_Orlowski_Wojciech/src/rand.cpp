#include "rand.hpp"
#include <cmath>

std::tuple<double,double> Normal(double mu, double sigma)
{
    double U1 = static_cast<double>(rand())/(static_cast<double>(RAND_MAX));
    double U2 = static_cast<double>(rand())/(static_cast<double>(RAND_MAX));
    
    double x = sigma*sqrt((-2. * log(1. - U1))) * cos(2.*M_PI*U2);
    double y = sigma*sqrt((-2. * log(1. - U1))) * sin(2.*M_PI*U2);
    
    return std::tuple<double,double>(x,y);
}
