//
// Created by Nevsky1_2 on 22.07.2021.
//

#pragma once
#ifndef NAVIE_STOKS_FUNCTIONS_H
#define NAVIE_STOKS_FUNCTIONS_H
#include "Polygonish.h"
#include "Mesh.h"

struct parameters {
    //main
    double ro, p, h, H, e, E, T, u, w, v;

    //transport
    double mu, Pr, la;

    double Cp, Gam, Gm;

    //vectors
    double* U, * U1;
    double* V; //ro, u, v, h
};

struct changes {
    double* dU;
};

struct Wall {

    int vel; //0 - Free clip, 1 - no clip
    int temp; //1 - Tw is set, 2 - qw is set
    double value; //значение Tw или qw

};

struct Zone {

    int grantype;
    //1. Wall. 2. Supersonic Inlet. 3. Symmetry
    //4. Supersonic Outlet. 5. Free boundary
    //6. Subsonic Inlet. 7. Subsonic Outlet.
    Wall* wall;

};

void Init(parameters* (&p), int nCells, int Nm);
void Viscous(parameters* p, changes* (&du), Mesh mesh, Cell* cells, double dt);
void GetParams(parameters* (&p),int nCells,int Nm);
void Tecplot(parameters* p, Cell* cells, int Nx, int Ny, int nCells);
void Convect(parameters* p, changes* (&u), Mesh mesh, Cell* cell, int It, double dt);
void Yw(Mesh mesh, Cell* (&cells), int nCells);
double Dist(Point A, Point B, Point E);
void SetGran(Mesh& mesh);

#endif //NAVIE_STOKS_FUNCTIONS_H
