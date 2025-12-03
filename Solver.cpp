#include "Solver.h"
#include "Vectors.h"
#include "Cylinder.h"
#include "Property.h"
#include "Grid.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

using namespace std;

Solver::Solver(Grid& grid, Property& fuidprop, Vectors& vec_vec, Cylinder& cyl)
    : Domain_grid(grid),
      fluid_property(fuidprop),
      field_vectors(vec_vec),
      Square_cylinder(cyl),
      initialdt(0.000001)
{
}

double Solver::getinitialdt() const {
    return initialdt;
}

/*****************  Boundary Condition  *****************/
void Solver::Boundary_condition(double dtt) {

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double dt = dtt;
    double Un = fluid_property.getinlet_average();

    // Top & bottom boundaries
    for (int i = 0; i <= Nx; ++i) {
        field_vectors.u[Ny-1][i] = field_vectors.u[Ny-2][i];
        field_vectors.u[0][i]   = field_vectors.u[1][i];
    }

    // Inlet/outlet for u
    for (int j = 0; j < Ny; ++j) {

        field_vectors.u[j][0] = Un;

        field_vectors.u_next[j][Nx] =
            field_vectors.u[j][Nx] -
            (field_vectors.u[j][Nx] - field_vectors.u[j][Nx-1]) *
            (Un * dt / Domain_grid.dx(Nx-1));
    }

    // v inlet / outlet
    for (int j = 0; j <= Ny; ++j) {
        field_vectors.v[j][0] = 0;

        field_vectors.v_next[j][Nx-1] =
            field_vectors.v[j][Nx-1] -
            (field_vectors.v[j][Nx-1] - field_vectors.v[j][Nx-2]) *
            (Un * dt / Domain_grid.dx(Nx-1));
    }

    // wall BCs for v
    for (int i = 0; i < Nx; ++i) {
        field_vectors.v[0][i]  = field_vectors.v[1][i];
        field_vectors.v[Ny][i] = field_vectors.v[Ny-1][i];
    }

    // cylinder no-slip
    for (int j = 1; j < Ny; ++j) {
        for (int i = 1; i < Nx; ++i) {

            if (Square_cylinder.is_inside(Domain_grid.getX(i), Domain_grid.getY(j))) {
                field_vectors.u[j][i] = 0;
                field_vectors.v[j][i] = 0;
                continue;
            }
        }
    }
}

/***************** Interior Velocity Update *****************/
void Solver::Interrior_Velocity(double timestep) {

    double dt = timestep;
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double rho = fluid_property.getrho();
    double mu  = fluid_property.getviscosity();

    /********** u predictor **********/
    for (int j = 1; j < Ny - 1; ++j) {
        for (int i = 1; i < Nx; ++i) {

            if (Square_cylinder.is_inside(Domain_grid.getX(i), Domain_grid.getY(j))) {
                field_vectors.u[j][i] = 0;
                field_vectors.uupp[j][i] = 0;
                field_vectors.v[j][i] = 0;
                field_vectors.vvpp[j][i] = 0;
                continue;
            }

            double vn = 0.5*(field_vectors.v[j+1][i-1] + field_vectors.v[j+1][i]);
            double un = 0.5*(field_vectors.u[j+1][i]   + field_vectors.u[j][i]);
            double vs = 0.5*(field_vectors.v[j][i-1]   + field_vectors.v[j][i]);
            double us = 0.5*(field_vectors.u[j-1][i]   + field_vectors.u[j][i]);
            double ue = 0.5*(field_vectors.u[j][i+1]   + field_vectors.u[j][i]);
            double uw = 0.5*(field_vectors.u[j][i-1]   + field_vectors.u[j][i]);

            double Ae = Domain_grid.dy(j);
            double Aw = Domain_grid.dy(j);
            double An = (Domain_grid.dx(i-1) + Domain_grid.dx(i)) / 2;
            double As = An;

            double m_e = rho * ue * Ae;
            double m_w = rho * uw * Aw;
            double m_n = rho * vn * An;
            double m_s = rho * vs * As;

            double U_N = field_vectors.u[j+1][i];
            double U_S = field_vectors.u[j-1][i];
            double U_E = field_vectors.u[j][i+1];
            double U_W = field_vectors.u[j][i-1];
            double U_P = field_vectors.u[j][i];

            double d_pe = Domain_grid.dx(i);
            double d_pw = Domain_grid.dx(i-1);
            double d_pn = (Domain_grid.dy(j) + Domain_grid.dy(j+1)) / 2;
            double d_ps = (Domain_grid.dy(j) + Domain_grid.dy(j-1)) / 2;

            double diff_e = mu*(U_E-U_P)*Ae/d_pe;
            double diff_w = mu*(U_P-U_W)*Aw/d_pw;
            double diff_n = mu*(U_N-U_P)*An/d_pn;
            double diff_s = mu*(U_P-U_S)*As/d_ps;

            double gamma = ((Domain_grid.dx(i-1)+Domain_grid.dx(i))/2) * Domain_grid.dy(j);

            double Ru_x = -(m_e*ue - m_w*uw + m_n*un - m_s*us)
                           + (diff_e - diff_w + diff_n - diff_s);

            field_vectors.uupp[j][i] = U_P + (dt/(rho*gamma))*Ru_x;
        }
    }

    /********** v predictor **********/
    for (int j = 1; j < Ny; ++j) {
        for (int i = 1; i < Nx-1; ++i) {

            if (Square_cylinder.is_inside(Domain_grid.getX(i), Domain_grid.getY(j))) {
                field_vectors.u[j][i] = 0;
                field_vectors.uupp[j][i] = 0;
                field_vectors.v[j][i] = 0;
                field_vectors.vvpp[j][i] = 0;
                continue;
            }

            double vn = 0.5*(field_vectors.v[j+1][i] + field_vectors.v[j][i]);
            double vs = 0.5*(field_vectors.v[j-1][i] + field_vectors.v[j][i]);
            double ue = 0.5*(field_vectors.u[j][i+1] + field_vectors.u[j-1][i+1]);
            double ve = 0.5*(field_vectors.v[j][i+1] + field_vectors.v[j][i]);
            double uw = 0.5*(field_vectors.u[j][i]   + field_vectors.u[j-1][i]);
            double vw = 0.5*(field_vectors.v[j][i-1] + field_vectors.v[j][i]);

            double Ae = (Domain_grid.dy(j)+Domain_grid.dy(j-1))/2;
            double Aw = Ae;
            double An = Domain_grid.dx(i);
            double As = An;

            double m_e = rho*ue*Ae;
            double m_w = rho*uw*Aw;
            double m_n = rho*vn*An;
            double m_s = rho*vs*As;

            double V_N = field_vectors.v[j+1][i];
            double V_S = field_vectors.v[j-1][i];
            double V_E = field_vectors.v[j][i+1];
            double V_W = field_vectors.v[j][i-1];
            double V_P = field_vectors.v[j][i];

            double d_pe = (Domain_grid.dx(i)+Domain_grid.dx(i+1))/2;
            double d_pw = (Domain_grid.dx(i)+Domain_grid.dx(i-1))/2;
            double d_pn = Domain_grid.dy(j);
            double d_ps = Domain_grid.dy(j-1);

            double diff_e = mu*(V_E-V_P)*Ae/d_pe;
            double diff_w = mu*(V_P-V_W)*Aw/d_pw;
            double diff_n = mu*(V_N-V_P)*An/d_pn;
            double diff_s = mu*(V_P-V_S)*As/d_ps;

            double gamma = Domain_grid.dx(i)*(Domain_grid.dy(j)+Domain_grid.dy(j-1))/2;

            double Ru_y = -(m_e*ve - m_w*vw + m_n*vn - m_s*vs)
                           + (diff_e - diff_w + diff_n - diff_s);

            field_vectors.vvpp[j][i] = V_P + (dt/(rho*gamma))*Ru_y;
        }
    }
}

/******************** Pressure Boundary Conditions ************************/
void Solver::pressure_boundary() {

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

            if (j == 0)
                field_vectors.P_next[j][i] = field_vectors.P_next[j+1][i];

            if (j == Ny-1)
                field_vectors.P_next[j][i] = field_vectors.P_next[j-1][i];

            if (i == 0)
                field_vectors.P_next[j][i] = field_vectors.P_next[j][i+1];

            if (i == Nx-1)
                field_vectors.P_next[j][i] = field_vectors.P_next[j][i-1];

            if (Square_cylinder.is_inside(Domain_grid.getX(i), Domain_grid.getY(j))) {

                double X = Domain_grid.getX(i);
                double Y = Domain_grid.getY(j);

                if (X == Square_cylinder.getXs() - Domain_grid.getD()/2 &&
                    Y >= Square_cylinder.getYs() - Domain_grid.getD()/2 &&
                    Y <= Square_cylinder.getYs() + Domain_grid.getD()/2)
                    field_vectors.P_next[j][i] = field_vectors.P_next[j][i-1];

                if (X == Square_cylinder.getXs() + Domain_grid.getD()/2 &&
                    Y >= Square_cylinder.getYs() - Domain_grid.getD()/2 &&
                    Y <= Square_cylinder.getYs() + Domain_grid.getD()/2)
                    field_vectors.P_next[j][i] = field_vectors.P_next[j][i+1];

                if (Y == Square_cylinder.getYs() - Domain_grid.getD()/2 &&
                    X >= Square_cylinder.getXs() - Domain_grid.getD()/2 &&
                    X <= Square_cylinder.getXs() + Domain_grid.getD()/2)
                    field_vectors.P_next[j][i] = field_vectors.P_next[j-1][i];

                if (Y == Square_cylinder.getYs() + Domain_grid.getD()/2 &&
                    X >= Square_cylinder.getXs() - Domain_grid.getD()/2 &&
                    X <= Square_cylinder.getXs() + Domain_grid.getD()/2)
                    field_vectors.P_next[j][i] = field_vectors.P_next[j+1][i];
            }
        }
    }
}

/********************** Gauss–Seidel ***************************/
void Solver::Guessiedel(double dtt) {

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double rho = fluid_property.getrho();
    double dt  = dtt;

    for (int iter = 0; iter < MaxIter; ++iter) {

        double maxerror = 0;

        for (int j = 1; j < Ny-1; ++j) {
            for (int i = 1; i < Nx-1; ++i) {

                if (Square_cylinder.is_inside(Domain_grid.getX(i), Domain_grid.getY(j)))
                    continue;

                double Ae = Domain_grid.dy(j);
                double Aw = Domain_grid.dy(j);
                double An = Domain_grid.dx(i);
                double As = Domain_grid.dx(i);

                double d_pe = (Domain_grid.dx(i)+Domain_grid.dx(i+1))/2;
                double d_pw = (Domain_grid.dx(i)+Domain_grid.dx(i-1))/2;
                double d_pn = (Domain_grid.dy(j)+Domain_grid.dy(j+1))/2;
                double d_ps = (Domain_grid.dy(j)+Domain_grid.dy(j-1))/2;

                double aE = Ae/d_pe;
                double aW = Aw/d_pw;
                double aN = An/d_pn;
                double aS = As/d_ps;

                double ap = aE + aW + aN + aS;

                double pN = field_vectors.P_next[j+1][i];
                double pS = field_vectors.P_next[j-1][i];
                double pE = field_vectors.P_next[j][i+1];
                double pW = field_vectors.P_next[j][i-1];

                double vpN = field_vectors.vvpp[j+1][i];
                double vpS = field_vectors.vvpp[j][i];
                double upE = field_vectors.uupp[j][i+1];
                double upW = field_vectors.uupp[j][i];

                double bp = -((rho*upE)*Ae - (rho*upW)*Aw + (rho*vpN)*An - (rho*vpS)*As)/dt;

                double PP = (aN*pN + aS*pS + aE*pE + aW*pW + bp)/ap;

                double newP = field_vectors.P[j][i]*(1-omega) + omega*PP;

                maxerror = max(maxerror, fabs(newP - field_vectors.P_next[j][i]));

                field_vectors.P_next[j][i] = newP;
            }
        }

        field_vectors.P = field_vectors.P_next;
        report_error = maxerror;

        if (maxerror < Max_error)
            break;
    }
}

/********************* Velocity Correction *************************/
void Solver::correction_velocity(double dtt) {

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double rho = fluid_property.getrho();
    double dt  = dtt;

    for (int j = 1; j < Ny-1; ++j) {
        for (int i = 1; i < Nx; ++i) {

            if (Square_cylinder.is_inside(Domain_grid.getX(i), Domain_grid.getY(j))) {
                field_vectors.u_next[j][i] = 0;
                continue;
            }

            field_vectors.u_next[j][i] =
                field_vectors.uupp[j][i] -
                dt * (field_vectors.P_next[j][i] -
                      field_vectors.P_next[j][i-1]) /
                (rho * (Domain_grid.dx(i)+Domain_grid.dx(i-1))/2);
        }
    }

    for (int j = 1; j < Ny; ++j) {
        for (int i = 1; i < Nx-1; ++i) {

            if (Square_cylinder.is_inside(Domain_grid.getX(i), Domain_grid.getY(j))) {
                field_vectors.v_next[j][i] = 0;
                continue;
            }

            field_vectors.v_next[j][i] =
                field_vectors.vvpp[j][i] -
                dt * (field_vectors.P_next[j][i] -
                      field_vectors.P_next[j-1][i]) /
                (rho * (Domain_grid.dy(j)+Domain_grid.dy(j-1))/2);
        }
    }
}

/********************* MAIN SOLVER LOOP ***************************/
void Solver::solve(double ultimate_time) {

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double rho = fluid_property.getrho();
    double dt  = initialdt;

    while (flowtime < ultimate_time) {

        field_vectors.u_old = field_vectors.u;
        field_vectors.v_old = field_vectors.v;

        field_vectors.u = field_vectors.u_next;
        field_vectors.v = field_vectors.v_next;

        Boundary_condition(dt);

        field_vectors.uupp = field_vectors.u;
        field_vectors.vvpp = field_vectors.v;

        Interrior_Velocity(dt);
        pressure_boundary();
        Guessiedel(dt);
        correction_velocity(dt);

        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i <= Nx; ++i)
                field_vectors.sumu[j][i] += field_vectors.u[j][i] * dt;

        dt = 0.00001;

        compute_lift_drag(Cd, Cdp);
        Monior(flowtime, dt);

        flowtime += dt;
    }

    cout << "=== Simulation Finished ===\n";
}


void Solver::compute_lift_drag(double& Cd, double& Cdp) {

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
 
    double rho = fluid_property.getrho();
    double mu  = fluid_property.getviscosity();
    double Uin = fluid_property.getinlet_average();
    double D   = Domain_grid.getD();
    

    double xC = Square_cylinder.getXs();
    double yC = Square_cylinder.getYs();
    double D2 = D / 2.0;

    double Fxtot = 0.0, Fx_p = 0.0, Fx_pplus=0,Fx_vplus=0;


  
    const double eps = 1e-3;

    for (int j = 1; j < Ny-1; ++j) {
        for (int i = 1; i < Nx-1; ++i) {

            double x = Domain_grid.getX(i);
            double y = Domain_grid.getY(j);

    
            bool left_face   = false;
            bool right_face  = false;
            bool bottom_face = false;
            bool top_face    = false;

            if (fabs(x - (xC - D2)) < eps && y >= yC - D2 - eps && y <= yC + D2 + eps)
                left_face = true;
            if (fabs(x - (xC + D2)) < eps && y >= yC - D2 - eps && y <= yC + D2 + eps)
                right_face = true;
            if (fabs(y - (yC - D2)) < eps && x >= xC - D2 - eps && x <= xC + D2 + eps)
                bottom_face = true;
            if (fabs(y - (yC + D2)) < eps && x >= xC - D2 - eps && x <= xC + D2 + eps)
                top_face = true;


            if (!(left_face || right_face || top_face || bottom_face))
                continue;

    
            double p_face = field_vectors.P[j][i];
            if (left_face) {
         
                p_face = 0.5 * (field_vectors.P[j][i] + field_vectors.P[j][i-1]);
            } else if (right_face) {
                p_face = 0.5 * (field_vectors.P[j][i] + field_vectors.P[j][i+1]);
            } else if (bottom_face) {
                p_face = 0.5 * (field_vectors.P[j][i] + field_vectors.P[j-1][i]);
            } else if (top_face) {
                p_face = 0.5 * (field_vectors.P[j][i] + field_vectors.P[j+1][i]);
            }

  
            double u = field_vectors.u[j][i];
            double v = field_vectors.v[j][i];

            
            double du_dx = 0.0, du_dy = 0.0, dv_dx = 0.0, dv_dy = 0.0;

            if (left_face) {
                du_dx = (field_vectors.u[j][i+1] - field_vectors.u[j][i]) / Domain_grid.dx(i);
            } else if (right_face) {
         
                du_dx = (field_vectors.u[j][i] - field_vectors.u[j][i-1]) / Domain_grid.dx(i-1);
            } else {
          
                du_dx = (field_vectors.u[j][i+1] - field_vectors.u[j][i-1]) / (2.0*Domain_grid.dx(i));
            }

            if (bottom_face) {
                dv_dy = (field_vectors.v[j+1][i] - field_vectors.v[j][i]) / Domain_grid.dy(j);
            } else if (top_face) {
                dv_dy = (field_vectors.v[j][i] - field_vectors.v[j-1][i]) / Domain_grid.dy(j-1);
            } else {
                dv_dy = (field_vectors.v[j+1][i] - field_vectors.v[j-1][i]) / (2.0*Domain_grid.dy(j));
            }

 // check here again
            if (j-1 >= 0 && j+1 < Ny) {
                du_dy = (field_vectors.u[j+1][i] - field_vectors.u[j-1][i]) / (2.0*Domain_grid.dy(j));
            } else if (j+1 < Ny) {
                du_dy = (field_vectors.u[j+1][i] - field_vectors.u[j][i]) / Domain_grid.dy(j);
            } else if (j-1 >= 0) {
                du_dy = (field_vectors.u[j][i] - field_vectors.u[j-1][i]) / Domain_grid.dy(j);
            }

            if (i-1 >= 0 && i+1 < Nx) {
                dv_dx = (field_vectors.v[j][i+1] - field_vectors.v[j][i-1]) / (2.0*Domain_grid.dx(i));
            } else if (i+1 < Nx) {
                dv_dx = (field_vectors.v[j][i+1] - field_vectors.v[j][i]) / Domain_grid.dx(i);
            } else if (i-1 >= 0) {
                dv_dx = (field_vectors.v[j][i] - field_vectors.v[j][i-1]) / Domain_grid.dx(i);
            }

            // Normal vector and face area
            double nx = 0.0, ny = 0.0, dS = 0.0;
            if (left_face)   { nx = -1.0; ny = 0.0; dS = Domain_grid.dy(j); }
            if (right_face)  { nx =  1.0; ny = 0.0; dS = Domain_grid.dy(j); }
            if (bottom_face) { nx =  0.0; ny = -1.0; dS = Domain_grid.dx(i); }
            if (top_face)    { nx =  0.0; ny =  1.0; dS = Domain_grid.dx(i); }


            double Fx_p_local = -p_face * nx * dS;
            Fx_pplus += Fx_p_local;


            double tau_xx = 2.0 * mu * du_dx;
            double tau_xy = mu * (du_dy + dv_dx);
            double tau_yy = 2.0 * mu * dv_dy;

          
            double Fx_v_local = (tau_xx * nx + tau_xy * ny) * dS;
            Fx_vplus += Fx_v_local;
        }
    }


    double Fxtot_all = Fx_pplus + Fx_vplus;


    cout << " Fx_p = " << Fx_pplus << "  Fx_v = " << Fx_vplus
         << "  Fxtot = " << Fxtot_all << endl;


    Cd = (2.0 * (Fx_vplus + Fx_pplus)) / (rho * Uin * Uin * D);
    Cdp = (2.0 * Fx_pplus) / (rho * Uin * Uin * D);

 

}

void Solver::generate_report() {
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();


    std::ofstream ContourV("Contour_V.txt");
    std::ofstream Contouru("Contour_u.txt");
    std::ofstream Contourv("Contour_v.txt");
    std::ofstream ContourP("Contour_P.txt");
    std::ofstream Recirculation("Recirculation2.txt");
     std::ofstream Recirculation_mean("Recirculation_mean.txt");

    if (!ContourV.is_open()) {
        std::cerr << "Error!\n";
        return;

        
    }

    double y_center=0, dy_center=0;
    // Contours
    for (int j = 0; j < Ny; ++j) {
        double y = Domain_grid.getY(j);
        for (int i = 0; i < Nx; ++i) {
            double x = Domain_grid.getX(i);
            y_center = Domain_grid.getY(j);
            if (y_center-10<1e-3)
            dy_center=Domain_grid.dy(j);

            double u = field_vectors.u[j][i];
            double v = field_vectors.v[j][i];
            double P = field_vectors.P[j][i];
            double V = std::sqrt(u*u + v*v);  // magnitude


            ContourV << x << "\t" << y << "\t" << V<< "\n";
            Contouru<< x << "\t" << y << "\t"<< u<< "\n";
            Contourv<< x << "\t" << y << "\t"<< v<< "\n";
            ContourP<< x << "\t" << y << "\t" << P<< "\n";
        }
        ContourV << "\n";
        Contouru << "\n";
        Contourv << "\n";
        ContourP << "\n";
    }
    int jj=static_cast<int> (10/dy_center);
    for (int i = 0; i <= Nx; ++i) {

        double uu = (field_vectors.u[jj][i]);
        double uuu = field_vectors.sumu[jj][i] / flowtime;
        double x = Domain_grid.getX(i);
        Recirculation << x << " " << uu << "\n";
       Recirculation_mean << x << " " << uuu << "\n";

    }
ContourV.close();
Contouru.close();
Contourv.close();
ContourP.close();
Recirculation.close();
Recirculation_mean.close();
    std::cout << "Contour data saved\n";
}

void Solver::Monior(double flowtime, double dt)
{
   
int Nx = Domain_grid.getNx();
int Ny = Domain_grid.getNy();
  
   
    double timestep=dt;


   cout << "time = " << setw(6) << flowtime<<'\t'<< "   dt= "<<timestep<<'\t'<<"   u = " << field_vectors.u[Ny/2][Nx/2] <<'\t'<<"Error=  "<< report_error<< endl;
 
ofstream forceFile("Forces.txt", ios::app);
forceFile << setw(12) << flowtime << setw(15) << Cd << setw(15) << Cdp << endl;
forceFile.close();
    
if (fmod(flowtime,100)<1e-6)
    {
cout << "time = " << setw(8) << flowtime<< "   Cd = " << setw(10) << Cd<< "   Cl = " << setw(10) << Cdp<<endl;
 forceFile << "time = " << setw(8) << flowtime<< "   Cd = " << setw(10) << Cdp << "   Cl = " << setw(10) << Cdp<<endl;
           
}
forceFile.close();

}
