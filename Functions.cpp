//
// Created by Nevsky1_2 on 22.07.2021.
//

#include "Functions.h"

void Init(parameters* (&p), int nCells, int Nm) {

    double T0, Cp, la, P0, Gm, Gam, R, ro, mu;

    double U0;

    //Входные данные
    T0 = 293.0; //K
    la = 2.5658e-2; //W/(m K)
    mu = 17.863e-6;

    U0 = 0.1;

    P0 = 101325.0; //Pa
    Gam = 1.4;
    Gm = 28.97;

    //расчёт
    R = 8314.41 / Gm; //J/(kg K)
    ro = P0 / (R * T0);
    Cp = Gam / (Gam - 1.) * R;

    for (int i = 0; i < nCells; ++i) { //ro, p, h, H, E, T, mu, u, v, w

        p[i].p = P0;
        p[i].ro = ro;
        p[i].u = U0;
        p[i].v = 0.;
        p[i].w = 0.;
        p[i].T = T0;

        p[i].Cp = Cp;
        p[i].la = la;
        p[i].mu = mu;
        p[i].Gam = Gam;
        p[i].Gm = Gm;

        p[i].Pr = p[i].mu * p[i].Cp / p[i].la;

        p[i].h = Cp * T0;
        double q2 = 0.5 * (p[i].u * p[i].u + p[i].v * p[i].v);
        p[i].H = p[i].h + q2;
        p[i].E = p[i].H - p[i].p / p[i].ro;

        p[i].e = p[i].h / p[i].Gam;

        p[i].U = new double [Nm];
        p[i].U1 = new double [Nm];

        p[i].V = new double [Nm];

        p[i].U1[0] = p[i].ro  * p[i].E;

        p[i].V[0] = p[i].h;

    }

}

void Viscous(parameters *p, changes *&du, Mesh mesh, Cell *cells, double dt) {

    int nFaces = mesh.Get_nFaces();

    for (int i = 0; i < nFaces; ++i) {

        Face face = mesh.Get_faces()[i];

        int cr = face.cr;
        int cl = face.cl;

        if (face.is_boundary) {

            int c = max(cl, cr);

            //Координаты ц.т. ячейки
            Point xc = cells[c].Get_c();

            //Координаты ц. грани
            Point xf = face.f_center;

            //расстояние до стенки
            double dx = (xc.x - xf.x);
            double dy = (xc.y - xf.y);
            double dl = sqrt(dx * dx + dy * dy);

            double length = face.length;             //длина грани
            double S = cells[c].Get_S();            //площадь ячейки c

            int z = face.zone;


            double Tw;
            if (face.zone == 2) Tw = 500.;

            if (face.zone == 4) Tw = 200.;

            if (face.zone == 2 || face.zone == 4) {

                double hw = p[c].Cp * Tw;

                //dh_dn
                double dh_dn = (p[c].h - hw) / dl;

                //mu/Pr
                double mu_Pr = p[c].mu / p[c].Pr;

                //Fv - поток через грань
                double Fv = mu_Pr * dh_dn;

                du[c].dU[0] += - Fv * length / S * dt;

            }

            if (face.zone == 1) {
                //Fv - поток через грань
                double Fv = 100.;
                du[c].dU[0] += Fv * length / S * dt;
            }

            if (face.zone == 3) {
                //Fv - поток через грань
                double Fv = -10.;
                du[c].dU[0] += Fv * length / S * dt;
            }


        } else {

            //Координаты левой и правой ячеек
            Point xr = cells[cr].Get_c();
            Point xl = cells[cl].Get_c();

            //расстояние между двумя ячейками
            double dx = (xr.x - xl.x);
            double dy = (xr.y - xl.y);
            double dl = sqrt(dx * dx + dy * dy);

            //dh/dn
            double dh_dn = (p[cr].h - p[cl].h) / dl;

            //среднее mu/Pr
            double mu_Pr = 0.5 * ((p[cr].mu + p[cl].mu) / p[cl].Pr);

            //Fv - поток через грань
            double Fv = mu_Pr * dh_dn;

            //приращения
            //если p[cr].h > p[cl].h => dh_dn > 0 => Fv > 0
            //тепло втекает в левую ячейку => для неё в приращении Fv берется с плюсом
            //для правой - с минусом

            double length = face.length;              //длина грани
            double Sr = cells[cr].Get_S();            //площадь правой ячейки
            double Sl = cells[cl].Get_S();            //площадь левой ячейки

            du[cr].dU[0] += - Fv * length / Sr * dt;
            du[cl].dU[0] +=  Fv * length / Sl * dt;

        }

    }


}

void GetParams(parameters *&p, int nCells, int Nm) {

    for (int i = 0; i < nCells; ++i) {

        p[i].E = p[i].U1[0] / p[i].ro;

        double q2 = 0.5 * sqrt(p[i].u * p[i].u + p[i].v * p[i].v);

        p[i].e = p[i].E - q2;

        p[i].h = p[i].e * p[i].Gam;

        p[i].T = p[i].h / p[i].Cp;

        p[i].H = p[i].h + q2;

        double R = 8314.41 / p[i].Gm;

        p[i].V[0] = p[i].h;

    }

}

void Tecplot(parameters *p, Cell *cells, int Nx, int Ny, int nCells) {

    //Создаём поток для записи
    string f = "Tecplot/T.plt";
    ofstream record(f, ios::out);

    if (record) {

        record << "VARIABLES = \"X\", \"Y\", \"T,K\"" << endl;

        record << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;

        for (int i = 0; i < nCells; ++i) {
            record << cells[i].Get_c().x << " " << cells[i].Get_c().y << " " << p[i].T << endl;
        }

    }
    record.close();

}

void Convect(parameters *p, changes *&du, Mesh mesh, Cell *cells, int It, double dt) {

    int nFaces = mesh.Get_nFaces();

    for (int i = 0; i < nFaces; ++i) {

        Face face = mesh.Get_faces()[i];

        int cr = face.cr;
        int cl = face.cl;

        //nodes
        int n1 = face.nodes[0];
        int n2 = face.nodes[1];

        double dl = face.length; //длина грани

        Point x1 = mesh.Get_nodes()[n1];
        Point x2 = mesh.Get_nodes()[n2];

        double nx = -(x2.y - x1.y) / dl;
        double ny = (x2.x - x1.x) / dl;

        //при таком подходе вектор n является внешним по отношению к левой ячейки и внутренним по отношению к правой

        double Fc = 0.;

        if (face.is_boundary) {

            int c = max(cl, cr);
            double S = cells[c].Get_S(); //площадь ячейки c

            if (face.zone == 1) {

                //входные параметры
                double uInlet, vInlet, T_Inlet;
                uInlet = p[c].u;
                vInlet = p[c].v;
                T_Inlet = 800.;

                double un = uInlet * nx + vInlet * ny;

                double h = p[c].Cp * T_Inlet;
                double H = h + 0.5 * (uInlet * uInlet + vInlet * vInlet);

                double E = h / p[c].Gam + 0.5 * (uInlet * uInlet + vInlet * vInlet);

                //Fc - поток через грань
                Fc = p[c].ro * H * un;

                du[c].dU[0] += -Fc * dl / S * dt;
            }

            if (face.zone == 3) {

                double un = p[c].u * nx + p[c].v * ny;
                //Fc - поток через грань
                Fc = p[c].ro * p[c].H * un;
                du[c].dU[0] += Fc * dl / S * dt;

            }

            if (face.zone == 2 || face.zone == 4) {
                du[c].dU[0] += 0.;
            }
        } else { //internal face

            //средние значения
            double u_, v_, un_, H_, E_, ro_;

            u_ = 0.5 * (p[cr].u + p[cl].u);
            v_ = 0.5 * (p[cr].v + p[cl].v);
            un_ = u_ * nx + v_ * ny;
            H_ = 0.5 * (p[cr].H + p[cl].H);
            E_ = 0.5 * (p[cr].E + p[cl].E);
            ro_ = 0.5 * (p[cr].ro + p[cl].ro);

            double A = H_ / E_ * un_;
            double Apl, Amn;
            Apl = 0.5 * (A + abs(A));
            Amn = 0.5 * (A - abs(A));

            double UL, UR;
            UL = p[cl].U[0];
            UR = p[cr].U[0];

            //Fc - поток через грань
            Fc  = UL * Apl + UR * Amn;

            double Sr = cells[cr].Get_S();
            double Sl = cells[cl].Get_S();

            du[cr].dU[0] += +Fc * dl / Sr * dt;
            du[cl].dU[0] += +Fc * dl / Sl * dt;
            
        }
    }
}

void Yw(Mesh mesh, Cell* (&cells), int nCells) {

    int nFaces = mesh.Get_nFaces();

    for (int k = 0; k < nCells; ++k) {

        double z1 = 1.e10;
        Point E = cells[k].Get_c();

        for (int i = 0; i < nFaces; ++i) {
            if (mesh.Get_faces()[i].is_boundary) {

                int n1 = mesh.Get_faces()[i].nodes[0];
                int n2 = mesh.Get_faces()[i].nodes[1];
                Point A = mesh.Get_nodes()[n1];
                Point B = mesh.Get_nodes()[n2];

                double z2 = Dist(A, B, E);

                if (z1 > z2) z1 = z2;
            }
        }
        cells[k].Set_Yw(z1);

    }

    //Создаём поток для записи
    int Nx = mesh.Get_Nx();
    int Ny = mesh.Get_Ny();

    string f = "Tecplot/Yw.plt";
    ofstream record2(f, ios::out);
    if (record2) {

        record2 << "VARIABLES = \"X\", \"Y\", \"Yw,m\"" << endl;

        record2 << "ZONE I= " << Ny - 1 << ", J= " << Nx - 1 << ", DATAPACKING=POINT" << endl;

        for (int i = 0; i < nCells; ++i) {
            record2 << cells[i].Get_c().x << " " << cells[i].Get_c().y << " " << cells[i].Get_Yw() << endl;
        }

    }
    record2.close();

}

double Dist(Point A, Point B, Point E) {

    //Расстояние от точки E до грани AB

    double AB[2], BE[2], AE[2];

    double AB_BE, AB_AE, x, y;
    double x1, x2, y1, y2, mod;

    double dist;

    //vector AB
    AB[0] = B.x - A.x;
    AB[1] = B.y - A.y;

    //vector BE
    BE[0] = E.x - B.x;
    BE[1] = E.y - B.y;

    //vector AE
    AE[0] = E.x - A.x;
    AE[1] = E.y - A.y;

    //Вычисление скалярного произведения
    AB_BE = (AB[0] * BE[0] + AB[1] * BE[1]);
    AB_AE = (AB[0] * AE[0] + AB[1] * AE[1]);

    //Случай 1
    if (AB_BE > 0) {

        y = E.y - B.y;
        x = E.x - B.x;
        dist = sqrt(x * x + y * y);

    } else if (AB_AE < 0){ //Случай 2

        y = E.y - A.y;
        x = E.x - A.x;
        dist = sqrt(x * x + y * y);

    } else { //Случай 3

        //Поиск высоты
        x1 = AB[0];
        y1 = AB[1];
        x2 = AE[0];
        y2 = AE[1];
        mod = sqrt(x1 * x1 + y1 * y1);
        dist = abs(x1 * y2 - x2 * y1) / mod;

    }
    return dist;
}

void SetGran(Mesh& mesh) {

    mesh.Get_zones()[0].grantype = 1; //Left gran = Wall
    mesh.Get_zones()[1].grantype = 1; //Bottom gran = Wall
    mesh.Get_zones()[2].grantype = 1; //Right gran = Wall
    mesh.Get_zones()[3].grantype = 1; //Top gran = Wall

    int nZones = mesh.Get_nZones();
    for (int i = 0; i < nZones; ++i) {

        int tp = mesh.Get_zones()[i].grantype;
        if (tp == 1) {
            mesh.Get_zones()[i].wall = new Wall[1];
        }
        if (tp == 2) {
            //mesh.Get_zones()[i].inlet = new Wall[1];
        }

        mesh.Get_zones()[0].wall[0].vel = 1;
        mesh.Get_zones()[0].wall[0].temp = 2; //qw
        mesh.Get_zones()[0].wall[0].value = 0; //=qw

        mesh.Get_zones()[1].wall[0].vel = 1;
        mesh.Get_zones()[1].wall[0].temp = 1; //Tw
        mesh.Get_zones()[1].wall[0].value = 500.; //=Tw

        mesh.Get_zones()[2].wall[0].vel = 1;
        mesh.Get_zones()[2].wall[0].temp = 2; //qw
        mesh.Get_zones()[2].wall[0].value = 0; //=qw

        mesh.Get_zones()[1].wall[0].vel = 1;
        mesh.Get_zones()[1].wall[0].temp = 1; //Tw
        mesh.Get_zones()[1].wall[0].value = 200.; //=Tw
    }
}
