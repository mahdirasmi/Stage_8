#ifndef PROPERTY_H
#define PROPERTY_H

#include <vector>

class Property {
public:
Property( double density, double Re_flow, Grid& grid);
        
    double getinflowProfile(int k) const { return inlet_Profile_velocity[k]; } 
    double getinlet_average() const { return Inlet_average_velocity; } 
    double getrho() const { return rho; }
    double getviscosity() const { return viscosity; }

    public:
    std::double rho;
    std::double viscosity;
    std::double Re;
    std::double Inlet_average_velocity;
    std::vector<double> inlet_Profile_velocity;
    std::Grid& Domain_grid;
    
};