#include "Grid.h"
#include "Property.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace std;

Property::Property(double density,double Re_flow,Grid& g)
    : rho(density),Re(Re_flow), grid(g), inletVel(1.0),
      viscosity(density*1*1/Re_flow), inlet_profile(g.getNy()+1,0.0){}


