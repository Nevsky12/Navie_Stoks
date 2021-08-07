#include "Functions.h"
#include "Mesh.h"

int main() {

    //Русификация
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);

    Mesh mesh;

    string filename = "grid_head_fine.dat";

    mesh.ReadingStruct(filename);

    int nCells = mesh.Get_nCells();

    auto cells = new Cell [nCells];

    mesh.CreateCell(cells);

    mesh.CreateFaces();

    mesh.CellFuncs(cells);

    remove("cells.txt");

    //Ввод исходных данных
    auto* p = new parameters [nCells];
    auto* du = new changes [nCells];

    //Инициализация
    int Nm = 1;

    Init(p, nCells, Nm);

    mesh.Set_Zones();
    SetGran(mesh);

    Yw(mesh, cells, nCells);

    int ItMax = 50000;
    int It = 0;

    double resMin = 1.e-6;

    double res = 1.;

    for (int i = 0; i < nCells; ++i) {
        du[i].dU = new double [Nm];
        //du[i].dU = 0.;
    }

    while (It < ItMax && res > resMin) {

        It++;
        for (int i = 0; i < nCells; ++i) {
            for (int j = 0; j < Nm; ++j) {
                p[i].U[j] = p[i].U1[j]; //переприсвоение
                du[i].dU[j] = 0.;
            }
        }

        double dt = 2.e-1;

        //Приращение за счёт вязкости
        Convect(p, du, mesh, cells, It, dt);

        //Приращение за счёт вязкости
        Viscous(p, du, mesh, cells, dt);

        for (int i = 0; i < nCells; ++i) {
            for (int j = 0; j < Nm; ++j) {
                p[i].U1[j] = p[i].U[j] + du[i].dU[j];
            }
        }

        GetParams(p, nCells, Nm);

        double res = 0.;
        for (int i = 0; i < nCells; ++i) {
            double res_ = abs(du[i].dU[0] / p[i].U1[0]);
            if (res < res_) res = res_;
        }

        cout << "It= " << It << ", res= " << res <<endl;

        int Nx = mesh.Get_Nx();
        int Ny = mesh.Get_Ny();

        Tecplot(p, cells, Nx, Ny, nCells);

    }
    return 0;
}
