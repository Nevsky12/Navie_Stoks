//
// Created by Asus on 16.07.2021.
//

#pragma once
#ifndef NAVIE_STOKS_POLYGONISH_H
#define NAVIE_STOKS_POLYGONISH_H
#include <iostream>
#include <windows.h>
#include <string>
#include <fstream>
using namespace std;

struct Point {
    double x, y;
};

class Polygonish {
public:
    Polygonish();
    ~Polygonish();
    Polygonish(Point* p_other, int n_other);

    void DataEntry(const string& filename);

    [[nodiscard]] int Get_n() const {return n;}
    Point Get_p(int i) {return p[i];}

    Point Center();
    double Square();

private:
    Point* p;
    int n;
};


#endif //NAVIE_STOKS_POLYGONISH_H
