#ifndef GRID_H
#define GRID_H

#include <vector>

class Grid {
public:
    Grid(int Nx1, int Nx2, int Nx3,
         int Ny1, int Ny2, int Ny3,
         double Ld, double Hd);

    double getX(int i) const;
    double getY(int j) const;
    double dx(int i) const;
    double dy(int j) const;
    int getNx() const;
    int getNy() const;

private:
    enum BlockIndex { LEFT = 0, MIDDLE = 1, RIGHT = 2 };

    std::vector<int> Nx, Ny;
    double L, H;

    std::vector<double> x_c, y_c;
    std::vector<double> dxx, dyy;

    std::vector<double> dl;
    std::vector<double> dH;

    double stretch(double s, double k);

    void compute_dxx();
    void compute_dyy();
    void build_coordinates();
};


#endif