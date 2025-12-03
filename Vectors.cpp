#include "Grid.h"
#include "Property.h"
#include "Cylinder.h"
#include "Vectors.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace std;
Vectors::Vectors(Grid& g,Property& f,Cylinder& C)
    : grid(g), fluid(f), cyl(C)
{
    Nx=grid.getNx();
    Ny=grid.getNy();

    u.assign(Ny, vector<double>(Nx+1,0.0));
    u_next=uupp=u_old=sumu=u;

    v.assign(Ny+1, vector<double>(Nx,0.0));
    v_next=vvpp=v_old=v;

    P.assign(Ny, vector<double>(Nx,0.0));
    P_next=P;
}

void Vectors::initilization_u(){
    for(int i=0;i<=Nx;i++){
        for(int j=0;j<Ny;j++){
            if(cyl.is_inside(grid.getX(i), grid.getY(j)))
                u[j][i]=0;
            else
                u[j][i]=fluid.getinlet_average();   
        }
    }
}

