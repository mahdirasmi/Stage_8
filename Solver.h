#ifndef SOLVER_H
#define SOLVER_H

#include "Grid.h"
#include "Property.h"
#include "Cylinder.h"
#include "Vectors.h"
#include <vector>

class Solver {

private:
    Grid& Domain_grid;
    Property& fluid_property;
    Vectors& field_vectors;
    Cylinder& Square_cylinder;

public:
    double flowtime = 0.0, initialdt;
    int MaxIter = 30;
    double Max_error = 1e-5;
    double dtsum = 0;
    double total_error = 0.0;
    double report_error = 0;
    double omega = 1.0;
    double Cd = 0.0, Cdp = 0.0;
    double CFL_max=0.05;
   

    Solver(Grid& grid, Property& fuidprop, Vectors& vec_vec, Cylinder& cyl);

    double getinitialdt() const;

    void Boundary_condition(double dtt);
    void Interrior_Velocity(double timestep);
    void pressure_boundary();
    void Guessiedel(double dtt);
    void correction_velocity(double dtt);
    void solve(double ultimate_time);

    void compute_lift_drag(double& Cd, double& Cdp);
    void generate_report();
    void Monior(double t, double dt);
    
};

#endif