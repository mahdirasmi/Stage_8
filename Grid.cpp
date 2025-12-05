#include "Grid.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace std;

Grid::Grid(int Nx1, int Nx2, int Nx3,
           int Ny1, int Ny2, int Ny3,
           double Ld, double Hd)
    : Nx{Nx1, Nx2, Nx3},
      Ny{Ny1, Ny2, Ny3},
      L(Ld), H(Hd),
      dxx(Nx1 + Nx2 + Nx3),
      dyy(Ny1 + Ny2 + Ny3),
      x_c(Nx1 + Nx2 + Nx3 + 1, 0.0),
      y_c(Ny1 + Ny2 + Ny3 + 1, 0.0),
      dl{ 0.326923 * L, 0.0384615 * L, 0.63461538 * L },
dH { 0.475 * H, 0.05 * H, 0.475 * H }

{
    compute_dxx();
    compute_dyy();
    build_coordinates();
}

double Grid::stretch(double s, double k) {
    return pow(0.5*s, k);
}

void Grid::compute_dxx() {
    int N_total = Nx[0] + Nx[1] + Nx[2];
    vector<double> w(N_total);

    for (int i = 0; i < Nx[0]; i++) {
        double s = 1.0 - double(i) / Nx[0];
        w[i] = stretch(s, 1.0);
    }

    for (int i = 0; i < Nx[1]; i++) {
        double xi_local = (double(i) + 0.5) / Nx[1];
        double symmetric_weight = 1-pow( 2*xi_local -1, 2.0);
     //   w[Nx[0] + i] = symmetric_weight;
      double m= (double (Nx[1])-1)/2;
        double A=2;
        double im=double (i)-m;
      w[Nx[0] + i] = w[Nx[0] -1]*(1+A*(1-pow(im,2)/pow(m,2)));
    }

    for (int i = 0; i < Nx[2]; i++) {
        double s = double(i + 1) / Nx[2];
        w[Nx[0] + Nx[1] + i] = stretch(s, 1.0);
    }

    double s[3] = {0.0};
    for (int i = 0; i < N_total; i++) {
        if (i < Nx[0]) s[0] += w[i];
        else if (i < Nx[0] + Nx[1]) s[1] += w[i];
        else s[2] += w[i];
    }

    for (int i = 0; i < Nx[0]; i++)
        dxx[i] = dl[0] / s[0] * w[i];
    for (int i = 0; i < Nx[1]; i++){
        dxx[Nx[0] + i] = dl[1] / s[1] * w[Nx[0] + i];
        double m= (double (Nx[1])-1)/2;
        double A=2;
        double im=double (i)-m;
        //dxx[Nx[0] + i] =  dxx[Nx[0] - 1]*(1+A*(1-pow(im,2)/pow(m,2)));
       // dxx[Nx[0] + i] =  dxx[Nx[0] - 1]*1;
    }
    for (int i = 0; i < Nx[2]; i++)
        dxx[Nx[0] + Nx[1] + i] = dl[2] / s[2] * w[Nx[0] + Nx[1] + i];
}

void Grid::compute_dyy() {
    int N_total = Ny[0] + Ny[1] + Ny[2];
    vector<double> w(N_total);

    for (int j = 0; j < Ny[0]; j++) {
        double s = 1.0 - double(j) / Ny[0];
        w[j] = stretch(s, 1.0);
    }

    for (int j = 0; j < Ny[1]; j++) {
        double xi_local = (double(j) + 0.5) / Ny[1];
        double symmetric_weight = 1.0 - pow(2.0 * xi_local - 1.0, 2.0);
        //w[Ny[0] + j] = 0.9 * symmetric_weight;
        w[Ny[0] + j] =  w[Ny[0] + -1];
    }

    for (int j = 0; j < Ny[2]; j++) {
        double s = double(j + 1) / Ny[2];
        w[Ny[0] + Ny[1] + j] = stretch(s, 1.0);
    }

    double s[3] = {0.0};
    for (int j = 0; j < N_total; j++) {
        if (j < Ny[0]) s[0] += w[j];
        else if (j < Ny[0] + Ny[1]) s[1] += w[j];
        else s[2] += w[j];
    }

    for (int j = 0; j < Ny[0]; j++)
        dyy[j] = dH[0] / s[0] * w[j];
    for (int j = 0; j < Ny[1]; j++)
    {
        dyy[Ny[0] + j] = dH[1] / s[1] * w[Ny[0] + j];
        double m= (double (Ny[1])-1)/2;
        double A=2;
        double im=double (j)-m;
        //dyy[Ny[0] + j] =  dyy[Ny[0] - 1]*(1+A*(1-pow(im,2)/pow(m,2)));
       // dyy[Ny[0] + j] =  dyy[Ny[0] - 1]*1;
    }
    for (int j = 0; j < Ny[2]; j++)
        dyy[Ny[0] + Ny[1] + j] = dH[2] / s[2] * w[Ny[0] + Ny[1] + j];
}

void Grid::build_coordinates() {
    for (int i = 1; i < x_c.size(); i++)
        x_c[i] = x_c[i - 1] + dxx[i - 1];

    for (int j = 1; j < y_c.size(); j++)
        y_c[j] = y_c[j - 1] + dyy[j - 1];
}



    void Grid::mesh_generation()
    {
           ofstream f("Mesh.dat");
        
  int Nx_total = Nx[0] + Nx[1] + Nx[2], Ny_total = Ny[0] + Ny[1] + Ny[2];

     cout<<"===== dx ====="<<endl;

    for(int j=0;j<=Ny_total;j++){
    for(int i=0;i<=Nx_total;i++)
        {
            f<< x_c[i]<<"\t"<<y_c[j] <<"\n";
         
        }
        f<<"\n";
        
    }
    f.close();
    for(int i=0;i<Nx_total;i++)
        {
            
            cout<<"i=  "<< i<<'\t'<<"dx=  "<<'\t'<<dxx[i]<<endl;
        }

    cout<<"===== dy ====="<<endl;

    for(int j=0;j<Ny_total;j++){
        cout<<"j=  "<< j<<'\t'<<"dy=  "<<'\t'<<dyy[j]<<endl;

    }

    }


double Grid::getX(int i) const { return x_c[i]; }
double Grid::getY(int j) const { return y_c[j]; }
double Grid::dx(int i) const { return dxx[i]; }
double Grid::dy(int j) const { return dyy[j]; }
int Grid::getNx() const { return x_c.size() - 1; }
int Grid::getNy() const { return y_c.size() - 1; }
double Grid::getD() const{return 1;};