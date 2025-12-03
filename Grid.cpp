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
      dl{8.5, 1.0, 16.5},
      dH{9.5, 1.0, 9.5}
{
    compute_dxx();
    compute_dyy();
    build_coordinates();
}

double Grid::stretch(double s, double k) {
    return pow(s, k);
}

void Grid::compute_dxx() {
    int N_total = Nx[0] + Nx[1] + Nx[2];
    vector<double> w(N_total);

    for (int i = 0; i < Nx[0]; i++) {
        double s = 1.0 - double(i) / Nx[0];
        w[i] = stretch(s, 2.0);
    }

    for (int i = 0; i < Nx[1]; i++) {
        double xi_local = (double(i) + 0.5) / Nx[1];
        double symmetric_weight = 1.0 - pow(2.0 * xi_local - 1.0, 2.0);
        w[Nx[0] + i] = 0.9 * symmetric_weight;
    }

    for (int i = 0; i < Nx[2]; i++) {
        double s = double(i + 1) / Nx[2];
        w[Nx[0] + Nx[1] + i] = stretch(s, 2.0);
    }

    double s[3] = {0.0};
    for (int i = 0; i < N_total; i++) {
        if (i < Nx[0]) s[0] += w[i];
        else if (i < Nx[0] + Nx[1]) s[1] += w[i];
        else s[2] += w[i];
    }

    for (int i = 0; i < Nx[0]; i++)
        dxx[i] = dl[0] / s[0] * w[i];
    for (int i = 0; i < Nx[1]; i++)
        dxx[Nx[0] + i] = dl[1] / s[1] * w[Nx[0] + i];
    for (int i = 0; i < Nx[2]; i++)
        dxx[Nx[0] + Nx[1] + i] = dl[2] / s[2] * w[Nx[0] + Nx[1] + i];
}

void Grid::compute_dyy() {
    int N_total = Ny[0] + Ny[1] + Ny[2];
    vector<double> w(N_total);

    for (int j = 0; j < Ny[0]; j++) {
        double s = 1.0 - double(j) / Ny[0];
        w[j] = stretch(s, 2.0);
    }

    for (int j = 0; j < Ny[1]; j++) {
        double xi_local = (double(j) + 0.5) / Ny[1];
        double symmetric_weight = 1.0 - pow(2.0 * xi_local - 1.0, 2.0);
        w[Ny[0] + j] = 0.9 * symmetric_weight;
    }

    for (int j = 0; j < Ny[2]; j++) {
        double s = double(j + 1) / Ny[2];
        w[Ny[0] + Ny[1] + j] = stretch(s, 2.0);
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
        dyy[Ny[0] + j] = dH[1] / s[1] * w[Ny[0] + j];
    for (int j = 0; j < Ny[2]; j++)
        dyy[Ny[0] + Ny[1] + j] = dH[2] / s[2] * w[Ny[0] + Ny[1] + j];
}

void Grid::build_coordinates() {
    for (size_t i = 1; i < x_c.size(); i++)
        x_c[i] = x_c[i - 1] + dxx[i - 1];

    for (size_t j = 1; j < y_c.size(); j++)
        y_c[j] = y_c[j - 1] + dyy[j - 1];
}

double Grid::getX(int i) const { return x_c[i]; }
double Grid::getY(int j) const { return y_c[j]; }
double Grid::dx(int i) const { return dxx[i]; }
double Grid::dy(int j) const { return dyy[j]; }
int Grid::getNx() const { return x_c.size() - 1; }
int Grid::getNy() const { return y_c.size() - 1; }