//
// Created by Nevsky1_2 on 17.07.2021.
//

#pragma once
#ifndef NAVIE_STOKS_CELL_H
#define NAVIE_STOKS_CELL_H
#include "Polygonish.h"

class Cell {
private:

    Point c; //координаты центра тяжести

    double S; //площадь ячейки

    int nFaces; //число граней, окружающих ячейку

    int nNodes; //число узлов, окружающих ячейку

    int* faces; //номера граней, окружающих ячейку (nFaces)

    int* nodes; //номера узлов, окружающих ячейку (nNodes)

    int* fType; //типы граней: = 0 -> внутренняя грань (nFaces)
                //             = 1 -> граничная грань

    int* cells; //номера соседних ячеек (nFaces)
                //fType = 0 -> номер ячейки
                //fType = 1 -Ю номер грани

    double Yw; //расстояние от центра ячейки до стенки

public:
    Cell();
    ~Cell();

    Point Get_c() {return c;}
    void Set_c(Point c_);

    void Set_S(double S_) {this->S = S_;}
    [[nodiscard]] double Get_S() const {return S;}

    void Set_Yw(double Yw_) {this->Yw = Yw_;}
    [[nodiscard]] double Get_Yw() const {return Yw;}

    void Set_nFaces(int nf);
    [[nodiscard]] int Get_nFaces() const {return nFaces;}

    void Set_nNodes(int nn);
    [[nodiscard]] int Get_nNodes() const {return nNodes;}

    int Get_Node(int i) {return nodes[i];}
    int Get_Face(int i);
    void Set_Nodes(int* nodes_, int nn);

    void Set_Face(int iFace, int fIndex) {faces[iFace] = fIndex;}

    void Set_fType(int iFace, int ftype) {fType[iFace] = ftype;}

    void Set_Cells(int iFace, int cell) {cells[iFace] = cell;}

    void Print(int i);
};


#endif //NAVIE_STOKS_CELL_H
