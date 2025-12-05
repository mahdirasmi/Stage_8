#ifndef CYLINDER_H
#define CYLINDER_H

#include "Grid.h"
#include "Property.h"

#include <vector>

class Cylinder {
private:
    Grid& grid;
public:
    double xsc, ysc, D;
    Cylinder(double xc,double yc,Grid& g);
    bool is_inside(double x,double y) const;
    double getXs() const { return xsc; }
    double getYs() const { return ysc; }
};

#endif
