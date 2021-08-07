//
// Created by Nevsky1_2 on 17.07.2021.
//

#include "Cell.h"

Cell::Cell() {
    this->c.x = this->c.y = 0.;
    this->S = 0.;
    this->nFaces = this->nNodes = 0;
    this->faces = new int [nFaces];
    this->nodes = new int [nNodes];
    this->fType = new int [nFaces];
    this->cells = new int [nFaces];
}

Cell::~Cell() {
    delete [] faces;
    delete [] nodes;
    delete [] fType;
    delete [] cells;
}

void Cell::Set_c(Point c_) {
    this->c.x = c_.x;
    this->c.y = c_.y;
}

void Cell::Set_nFaces(int nf) {
    this->nFaces = nf;
    this->cells = new int [nFaces];
    this->faces = new int [nFaces];
    this->fType = new int [nFaces];
}

void Cell::Set_nNodes(int nn) {
    this->nNodes = nn;
    this->nodes = new int [nNodes];
}

int Cell::Get_Face(int i) {
    if (i < nFaces) return faces[i];
    else {
        cout << "This number is out of range!" << endl;
        return -10;}
}

void Cell::Set_Nodes(int *nodes_, int nn) {
    this->nNodes = nn;
    this->nodes = new int [nNodes]; //номера узлов, окружающих ячейку
    for (int i = 0; i < nNodes; ++i) {
        this->nodes[i] = nodes_[i];
    }
}

void Cell::Print(int m) {
    //Создаем поток для записи
    string f = "cells.txt";
    ofstream record(f,ios::out | ios::app);

    if (record) {
        record << "Cell index: " << m << endl;
        record << "Cell center: " << c.x << " ," << c.y << endl;
        record << "Cell square: " << S << endl;
        record << "Cell number of faces: " << nFaces << endl;
        record << "Cell face indexes and types: " << endl;
        for (int i = 0; i < nFaces; ++i) {
            record << "index: " << faces[i] << ", type: " << fType[i] << endl;
        }
        record << "Cell number of nodes: " << nNodes << endl;
        record << "Node indexes: " << endl;
        for (int i = 0; i < nNodes; ++i) {
            record << nodes[i];
            if (i < nNodes - 1) cout << ", ";
        }
        record << endl;
        record << "Neighbours indexes: " << endl;
        for (int i = 0; i < nFaces; ++i) {
            record << cells[i];
            if (i < nFaces - 1) cout << ", ";
        }
        record << endl;

        record << "****************************************" << endl;
    } else cout << "This file wasn't opened!" << endl;
    record.close();
}
