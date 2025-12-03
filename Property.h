#ifndef PROPERTY_H
#define PROPERTY_H

#include "Grid.h"
#include <vector>

class Property {
private:
    double rho, viscosity, Re;
    double inletVel;
    std::vector<double> inlet_profile;
    Grid& grid;

public:
    Property(double density,double Re_flow,Grid& g);

    double getinflowProfile(int k) const { return inlet_profile[k]; }
    double getinlet_average() const { return inletVel; }
    double getrho() const { return rho; }
    double getviscosity() const { return viscosity; }
};

#endif