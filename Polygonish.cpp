//
// Created by Nevsky1_2 on 16.07.2021.
//

#include "Polygonish.h"

Polygonish::Polygonish() {
    this->n = 0;
    this->p = new Point[0];
}

Polygonish::~Polygonish() {
    delete [] p;
}

Polygonish::Polygonish(Point *p_other, int n_other) {
    this->n = n_other;
    this->p = new Point[n];
    for (int i = 0;i < n; ++i) {
        p[i] = p_other[i];
    }
}

void Polygonish::DataEntry(const string& filename) {
    ifstream reading(filename);
    if (reading) {
        reading >> n;
        p = new Point[n];
        for (int i = 0;i < n; ++i) {
            reading >> p[i].x >> p[i].y;
        }
    } else cout << "Error!" << endl;
}

Point Polygonish::Center() {
    int m = ++n;
    auto* xn = new Point[m];
    xn[n].x = xn[0].x;
    xn[n].y = xn[0].y;
    double A = 0;
    double Cx = 0., Cy = 0.;

    for (int i = 0;i < n; ++i) {
        A += xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y;
        Cx += (xn[i].x + xn[i+1].x) * (xn[i].x * xn[i+1].y - xn[i+1].x * xn[i].y);
        Cy += (xn[i].y + xn[i+1].y) * (xn[i].x * xn[i+1].y - xn[i+1].x * xn[i].y);
    }
    A = 0.5 * A;
    Cx = Cx / (6. * A);
    Cy = Cy / (6. * A);
    Point p1{};
    p1.x = Cx;
    p1.y = Cy;

    return p1;
}

double Polygonish::Square() {
    int m = ++n;
    auto* xn = new Point[m];
    xn[n].x = xn[0].x;
    xn[n].y = xn[0].y;
    double A = 0;

    for (int i = 0;i < n; ++i) {
        A += xn[i].x * xn[i + 1].y - xn[i + 1].x * xn[i].y;
    }
    A = 0.5 * fabs(A);
    return A;
}
