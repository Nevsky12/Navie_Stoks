//
// Created by Nevsky1_2 on 18.07.2021.
//

#pragma once
#ifndef NAVIE_STOKS_MESH_H
#define NAVIE_STOKS_MESH_H
#include "Cell.h"
#include "Functions.h"

struct Face {
    int nodes[2];           //индексы узлов грани (2)
    int cr;                 //индекс правой ячейки
    int cl;                 //индекс левой ячейки
    bool is_boundary;       //является ли эта грань граничной
    Point f_center;         //координаты центра
    double length;          //длина грани
    int zone;               //номер зоны, на которой устанавливается граничное условие

    void Print(int k) {
        cout << "Face index: " << k << endl;
        cout << "Node indexes: " << nodes[0] << ", " << nodes[1] << endl;
        cout << "Cr index: " << cr << ", Cl index: " << cl << endl;
        cout << "Center: " << f_center.x << ", " << f_center.y << endl;
        cout << "Length: " << length << endl;
    }
};

class Mesh {
public:

    Mesh();
    ~Mesh();

    [[nodiscard]] int Get_Nx() const {return Nx;}
    [[nodiscard]] int Get_Ny() const {return Ny;}
    [[nodiscard]] int Get_nNodes() const {return nNodes;}
    [[nodiscard]] int Get_nCells() const {return nCells;}
    [[nodiscard]] int Get_nFaces() const {return nFaces;}
    Point* Get_nodes() {return nodes;}
    Face* Get_faces() {return faces;}
    int Get_nZones() {return nZones;}
    Zone* Get_zones() {return zones;}
    void Set_Nx(int Nx_) {this->Nx = Nx_;}
    void Set_Ny(int Ny_) {this->Ny = Ny_;}
    void Set_nNodes(int nNodes_) {this->nNodes = nNodes_;}
    void Set_nCells(int nCells_) {this->nCells = nCells_;}
    void Set_nFaces(int nFaces_) {this->nFaces = nFaces_;}
    void Set_Zones() {this->nZones = 4; this->zones = new Zone[nZones];}

    void ReadingStruct(const string& filename);
    void CreateCell(Cell* (&cell));
    void CreateFaces();
    void CellFuncs(Cell* (&cells));

private:

    int nZones;
    Zone* zones;
    Face* faces;
    int Nx, Ny;                 //чисто точек по x,y соответственно в структурированной сетке
    int nNodes, nCells, nFaces; //число узлов, ячеек и граней соответственно
    Point* nodes;               //координаты узлов

};


#endif //NAVIE_STOKS_MESH_H
