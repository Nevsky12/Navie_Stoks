//
// Created by Nevsky1_2 on 18.07.2021.
//

#include "Mesh.h"

Mesh::Mesh() {
    this->Nx = this->Ny = this->nNodes = this->nFaces = this->nCells = 0;
    this->nodes = new Point [nNodes];
}

Mesh::~Mesh() {
    delete [] nodes;
}

void Mesh::ReadingStruct(const string& filename) {
    ifstream reading(filename);
    int temp;
    if (reading) {
        reading >> this->Nx >> this->Ny >> temp >> temp;

        this->nNodes = Nx * Ny;
        this->nCells = (Nx - 1) * (Ny - 1);
        this->nFaces = Nx * (Ny - 1) + Ny * (Nx - 1);
        this->nodes = new Point [nNodes];

        int k = 0;

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; i < Ny; ++j) {

                reading >> temp >> temp >> this->nodes[k].x >> this->nodes[k].y;
                ++k;
            }
        }
        cout << "Data from this file was read!"  << endl;
    } else {
        cout << "This file wasn't opened!" << endl;
    }

}

void Mesh::CreateCell(Cell *&cells) {
    int k = 0; //счётчик граней
    for (int i = 0; i < Nx - 1; ++i) {
        for (int j = 0; j < Ny - 1; ++j) {

            int nc = (Ny - 1) * i + j; //номер ячейки
            int nNodes_ = 4; //число узлов вокруг ячейки
            int* nodes_ = new int [nNodes_]; //создаём массив индексов узлов, расположенных вокруг ячейки, проходя узлы против часовой стрелки

            int n2 = Ny * i + j; //элемент под номером 0
            nodes_[0] = n2;
            n2 = Ny * (i + 1) + j; //элемент под номером 1
            nodes_[1] = n2;
            n2 = Ny * (i + 1) + (j + 1); //элемент под номером 2
            nodes_[2] = n2;
            n2 = Ny * i  + (j + 1); //элемент под номером 3
            nodes_[3] = n2;

            cells[nc].Set_Nodes(nodes_,nNodes_);
        }
    }
}

void Mesh::CreateFaces() {
    faces = new Face[nFaces];

    int k = 0; //счётчик граней
    //вертикальные грани
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny - 1; ++j) {

            int n1 = Ny * i + j;
            int n2 = Ny * i + ++j;
            faces[k].nodes[0] = n1;
            faces[k].nodes[1] = n2;

            if (i == 0) {

                faces[k].is_boundary = true;
                faces[k].cr = -1;
                faces[k].cl = (Ny - 1) * i + j;
                faces[k].zone = 1;

            } else if (i == Nx - 1) {

                faces[k].is_boundary  = true;
                faces[k].cl = -1;
                faces[k].cr = (Ny - 1) * (i - 1) + j;
                faces[k].zone = 3;

            } else {

                faces[k].is_boundary = false;
                faces[k].cl = (Ny - 1) * i + j;
                faces[k].cr = (Ny - 1) * (i - 1) + j;
                faces[k].zone = -1;

            }

            faces[k].f_center.x = 0.5 * (nodes[n1].x + nodes[n2].x);
            faces[k].f_center.y = 0.5 * (nodes[n1].y + nodes[n2].y);

            double dx = (nodes[n1].x - nodes[n2].x);
            double dy = (nodes[n1].y - nodes[n2].y);

            faces[k].length = sqrt(dx * dx + dy * dy);

            //faces[k].Print(k);

            ++k;
        }
    }
    //горизонтальные грани
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx - 1; ++i) {

            int n1 = Ny * i + j;
            int n2 = Ny * (i + 1) + ++j;
            faces[k].nodes[0] = n1;
            faces[k].nodes[1] = n2;

            if (j == 0) {

                faces[k].is_boundary = true;
                faces[k].cl = -1;
                faces[k].cr = (Ny - 1) * i + j;
                faces[k].zone = 2;

            } else if (j == Ny - 1) {

                faces[k].is_boundary  = true;
                faces[k].cr = -1;
                faces[k].cl = (Ny - 1) * (i - 1) + j;
                faces[k].zone = 4;

            } else {

                faces[k].is_boundary = false;
                faces[k].cl = (Ny - 1) * i + j;
                faces[k].cr = (Ny - 1) * i + (j - 1);
                faces[k].zone = -1;

            }

            faces[k].f_center.x = 0.5 * (nodes[n1].x + nodes[n2].x);
            faces[k].f_center.y = 0.5 * (nodes[n1].y + nodes[n2].y);

            double dx = (nodes[n1].x - nodes[n2].x);
            double dy = (nodes[n1].y - nodes[n2].y);

            faces[k].length = sqrt(dx * dx + dy * dy);

            //faces[k].Print(k);

            ++k;
        }
    }
}

void Mesh::CellFuncs(Cell *&cells) {
    for (int i = 0; i < nCells; ++i) {
        //массив точек - узлов вокруг ячейки
        int mNodes = cells[i].Get_nNodes();

        int* nds = new int [mNodes]; //массив номеров узлов
        auto* pnts = new Point [mNodes]; //массив координат узлов

        for (int j = 0; j < mNodes; ++j) {

            int n_ = cells[i].Get_Node(j);
            nds[j] = n_;

            //координаты узлов
            pnts[j].x = nodes[n_].x;
            pnts[j].y = nodes[n_].y;

        }
        //создаем многоугольник
        Polygonish pl(pnts, mNodes);

        Point center = pl.Center();
        double S = pl.Square();

        cells[i].Set_c(center);
        cells[i].Set_S(S);
        cells[i].Set_nFaces(mNodes);
    }

    for (int k = 0; k < nFaces; ++k) {

        int n1 = faces[k].nodes[0];
        int n2 = faces[k].nodes[1];

        int cl = faces[k].cl;
        int cr = faces[k].cr;

        int ftype = 0;
        if (faces[k].is_boundary) ftype = 1;

        if (cr >= 0) {

            //массив узлов ячейки cr
            int nn = cells[cr].Get_nNodes(); //размер массива узлов
            int* nncr = new int [nn];

            for (int  i = 0; i < nn; ++i) {
                nncr[i] = cells[cr].Get_Node(i);
            }

            int iFace = 0;

            for (int  i = 0; i < nn; ++i) {

                int j = ++i;
                if (j == nn) j = 0;
                if (nncr[i] == n1 && nncr[j] == n2) {

                    iFace = i;
                    cells[cr].Set_Face(iFace, k);
                    cells[cr].Set_fType(iFace, ftype);
                    cells[cr].Set_Cells(iFace, cl);

                }

            }
        }

        if (cl >= 0) {

            //массив узлов ячейки cl
            int nn = cells[cl].Get_nNodes(); //размер массива узлов
            int* nncl = new int [nn];

            for (int  i = 0; i < nn; ++i) {
                nncl[i] = cells[cl].Get_Node(i);
            }

            int iFace = 0;

            for (int  i = 0; i < nn; ++i) {

                int j = ++i;
                if (j == nn) j = 0;
                if (nncl[i] == n2 && nncl[j] == n1) {

                    iFace = i;
                    cells[cl].Set_Face(iFace, k);
                    cells[cl].Set_fType(iFace, ftype);
                    cells[cl].Set_Cells(iFace, cl);

                }

            }
        }
    }
}
