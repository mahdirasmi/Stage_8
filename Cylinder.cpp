#include "Grid.h"
#include "Property.h"
#include "Cylinder.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace std;

  
Cylinder::Cylinder(double xc,double yc,Grid& g)
    : xsc(xc), ysc(yc), grid(g), D(1) {}

bool Cylinder::is_inside(double x,double y) const {
    return (x>=xsc-D/2 && x<=xsc+D/2 && y>=ysc-D/2 && y<=ysc+D/2);
}


