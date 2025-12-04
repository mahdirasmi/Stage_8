#include "Grid.h"
#include "Property.h"
#include "Cylinder.h"
#include "Vectors.h"
#include "Solver.h"
#include <iostream>

int main()
{

    Grid grid(
        3*15, 25, 3*25,       
        3*15, 25, 3*15,       
        26.0, 20.0);    

grid.mesh_generation();


  Property property(1, 100.0, grid);
  Cylinder cylinder(9, 10.0, grid);
 
   
    Vectors vectors(grid,property,cylinder);
    vectors.initilization_u();
    vectors.u_old = vectors.u;
    vectors.uupp = vectors.u; 
    vectors.u_next = vectors.u;
    vectors.v_old = vectors.v;
    vectors.vvpp = vectors.v;
    vectors.v_next = vectors.v;
   
  
  
   Solver solver(grid, property, vectors, cylinder);

   
   solver.solve(50);
  solver.generate_report();

   std:: cout << "\n✅ Simulation complete!\n";

 
    return 0;
}