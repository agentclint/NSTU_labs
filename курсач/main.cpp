#include "Solver.h"
#include <fstream>
#include <iostream>
#include <vector>
struct node
{
    double r;
    double z;
};
struct material
{
    double lambda;
    int gamma_id;
};
struct element
{
    std::vector<int> node_loc;
    int mater;
    int f_id;
};
std::vector<element> all_elems;
std::vector<material> all_materials;
std::vector<std::vector<int>> S1;
double gamma(double r, double z, int gam_id)
{
    switch (gam_id)
    {
    case 0:
        return 1;
    default:
        break;
    }
}
double func_f(double r, double z, int f_id)
{
    switch (f_id)
    {
    case 0:
        return r * z - z / r;
    default:
        break;
    }
}

double func_S1(double r, double z, int s1_id)
{
    switch (s1_id)
    {
    case 0:
        return r * z;
    default:
        break;
    }
}

int Input()
{
    int N, Nmat, Kel, NS1;
    std::ifstream in;
    in.open("info.txt");
    in >> N >> Nmat >> Kel >> NS1;
    in.close();
    in.open("rz.txt");
    all_nodes.resize(N);
    for (int i = 0; i < N; i++)
    {
        in >> all_nodes[i].r >> all_nodes[i].z;
    }
    in.close();
    in.open("S1.txt");
    S1.resize(NS1);
    for (int i = 0; i < NS1; i++)
    {
        int size;
        in >> size;
        S1[i].resize(size);
        for (int j = 0; j < size; j++)
        {
        }
        in >> S1[i][j];
    }
    in.close();
    in.open("material.txt");
    all_materials.resize(Nmat);
    for (int i = 0; i < Nmat; i++)
    {
        in >> all_materials[i].lambda >> all_materials[i].gamma_id;
    }
    in.close();
    in.open("elem.txt");
    all_elems.resize(Kel);
    for (int i = 0; i < Kel; i++)
    {
        all_elems[i].node_loc.resize(4);
        in >> all_elems[i].node_loc[0] >> all_elems[i].node_loc[1] >> all_elems[i].node_loc[2] >> all_elems[i].node_loc[3] >> all_elems[i].mater >> all_elems[i].f_id;
    }
    in.close();
    return 0;
}

double GetG_Loc(double rp, double lambda, double hr, double hz,
                std::vector<std::vector<double>> &G_loc)
{
    double a1 = (lambda * hz * rp) / (6 * hr),
           a2 = (lambda * hz) / (12),
           a3 = (lambda * hr * rp) / (6 * hz),
           a4 = (lambda * hr * hr) / (12 * hz);
    G_loc[0][0] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
    G_loc[0][1] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
    G_loc[0][2] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
    G_loc[0][3] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
    G_loc[1][0] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
    G_loc[1][1] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
    G_loc[1][2] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
    G_loc[1][3] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;
    G_loc[2][0] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
    G_loc[2][1] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
    G_loc[2][2] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
    G_loc[2][3] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
    G_loc[3][0] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
    G_loc[3][1] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;
    G_loc[3][2] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
    G_loc[3][3] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
    return 0;
}

double GetM_Loc(double rp, double zs, int gam, double hr, double hz, std::vector<std::vector<double>> &M_loc)
{
    double g1 = gamma(rp, zs, gam),
           g2 = gamma(rp + hr, zs, gam), g3 = gamma(rp, zs + hz, gam),
           g4 = gamma(rp + hr, zs + hz, gam);
    M_loc[0][0] += hr * (g1 * (rp / 4 + hr / 20) * hz / 4 +
                         g2 * (rp / 12 + hr / 30) * hz / 4 + g3 * (rp / 4 + hr / 20) * hz / 12 +
                         g4 * (rp / 12 + hr / 30) * hz / 12);
    M_loc[0][1] += hr * (g1 * (rp / 12 + hr / 30) * hz / 4 + g2 * (rp / 12 + hr / 20) * hz / 4 +
                         g3 * (rp / 12 + hr / 30) * hz / 12 + g4 * (rp / 12 + hr / 20) * hz / 12);
    M_loc[0][2] += hr * (g1 * (rp / 4 + hr / 20) * hz / 12 +
                         g2 * (rp / 12 + hr / 30) * hz / 12 +
                         g3 * (rp / 4 + hr / 20) * hz / 12 +
                         g4 * (rp / 12 + hr / 30) * hz / 12);
    M_loc[0][3] += hr * (g1 * (rp / 12 + hr / 30) * hz / 12 + g2 * (rp / 12 + hr / 20) * hz / 12 + g3 * (rp / 12 + hr / 30) * hz / 12 + g4 * (rp / 12 + hr / 20) * hz / 12);
    M_loc[1][0] += hr * (g1 * (rp / 12 + hr / 30) * hz / 4 + g2 * (rp / 12 + hr / 20) * hz / 4 +
                         g3 * (rp / 12 + hr / 30) * hz / 12 + g4 * (rp / 12 + hr / 20) * hz / 12);
    M_loc[1][1] += hr * (g1 * (rp / 12 + hr / 20) * hz / 4 +
                         g2 * (rp / 4 + hr / 5) * hz / 4 +
                         g3 * (rp / 12 + hr / 20) * hz / 12 +
                         g4 * (rp / 4 + hr / 5) * hz / 12);
    M_loc[1][2] += hr * (g1 * (rp / 12 + hr / 30) * hz / 12 + g2 * (rp / 12 + hr / 20) * hz / 12 + g3 * (rp / 12 + hr / 30) * hz / 12 + g4 * (rp / 12 + hr / 20) * hz / 12);
    M_loc[1][3] += hr * (g1 * (rp / 12 + hr / 20) * hz / 12 +
                         g2 * (rp / 4 + hr / 5) * hz / 12 +
                         g3 * (rp / 12 + hr / 20) * hz / 12 +
                         g4 * (rp / 4 + hr / 5) * hz / 12);
    M_loc[2][0] += hr * (g1 * (rp / 4 + hr / 20) * hz / 12 +
                         g2 * (rp / 12 + hr / 30) * hz / 12 +
                         g3 * (rp / 4 + hr / 20) * hz / 12 +
                         g4 * (rp / 12 + hr / 30) * hz / 12);
    M_loc[2][1] += hr * (g1 * (rp / 12 + hr / 30) * hz / 12 + g2 * (rp / 12 + hr / 20) * hz / 12 + g3 * (rp / 12 + hr / 30) * hz / 12 + g4 * (rp / 12 + hr / 20) * hz / 12);
    M_loc[2][2] += hr * (g1 * (rp / 4 + hr / 20) * hz / 12 +
                         g2 * (rp / 12 + hr / 30) * hz / 12 +
                         g3 * (rp / 4 + hr / 20) * hz / 4 +
                         g4 * (rp / 12 + hr / 30) * hz / 4);
    M_loc[2][3] += hr * (g1 * (rp / 12 + hr / 30) * hz / 12 + g2 * (rp / 12 + hr / 20) * hz / 12 +
                         g3 * (rp / 12 + hr / 30) * hz / 4 + g4 * (rp / 12 + hr / 20) * hz / 4);
    M_loc[3][0] += hr * (g1 * (rp / 12 + hr / 30) * hz / 12 + g2 * (rp / 12 + hr / 20) * hz / 12 + g3 * (rp / 12 + hr / 30) * hz / 12 + g4 * (rp / 12 + hr / 20) * hz / 12);
    M_loc[3][1] += hr * (g1 * (rp / 12 + hr / 20) * hz / 12 +
                         g2 * (rp / 4 + hr / 5) * hz / 12 +
                         g3 * (rp / 12 + hr / 20) * hz / 12 +
                         g4 * (rp / 4 + hr / 5) * hz / 12);
    M_loc[3][2] += hr * (g1 * (rp / 12 + hr / 30) * hz / 12 + g2 * (rp / 12 + hr / 20) * hz / 12 +
                         g3 * (rp / 12 + hr / 30) * hz / 4 + g4 * (rp / 12 + hr / 20) * hz / 4);
    M_loc[3][3] += hr * (g1 * (rp / 12 + hr / 20) * hz / 12 +
                         g2 * (rp / 4 + hr / 5) * hz / 12 +
                         g3 * (rp / 12 + hr / 20) * hz / 4 +
                         g4 * (rp / 4 + hr / 5) * hz / 4);
    return 0;
}

int Getb_Loc(double rp, double zs, double hr, double hz,
             std::vector<double> &b_loc, int f_id)
{
    double f1 = func_f(rp, zs, f_id),
           f2 = func_f(rp + hr, zs, f_id), f3 = func_f(rp, zs + hz, f_id),
           f4 = func_f(rp + hr, zs + hz, f_id);
    b_loc[0] =
        f1 * (hr * hz / 3 * (rp / 3 + hr / 12)) + f2 * (hr * hz / 3 * (rp / 6 + hr / 12)) + f3 * (hr * hz / 6 * (rp / 3 + hr / 12)) +
        f4 * (hr * hz / 6 * (rp / 6 + hr / 12));
    b_loc[1] =
        f1 * (hr * hz / 3 * (rp / 6 + hr / 12)) +
        f2 * (hr * hz / 3 * (rp / 3 + hr / 4)) +
        f3 * (hr * hz / 6 * (rp / 6 + hr / 12)) +
        f4 * (hr * hz / 6 * (rp / 3 + hr / 4));
    b_loc[2] =
        f1 * (hr * hz / 6 * (rp / 3 + hr / 12)) + f2 * (hr * hz / 6 * (rp / 6 + hr / 12)) + f3 * (hr * hz / 3 * (rp / 3 + hr / 12)) +
        f4 * (hr * hz / 3 * (rp / 6 + hr / 12));
    b_loc[3] =
        f1 * (hr * hz / 6 * (rp / 6 + hr / 12)) +
        f2 * (hr * hz / 6 * (rp / 3 + hr / 4)) +
        f3 * (hr * hz / 3 * (rp / 6 + hr / 12)) +
        f4 * (hr * hz / 3 * (rp / 3 + hr / 4));
    return 0;
}

int Get_Loc(std::vector<std::vector<double>> &A_loc, std::vector<double> &b_loc,
            int el_id)
{
    element el = all_elems[el_id];
    double hr = all_nodes[el.node_loc[1]].r - all_nodes[el.node_loc[0]].r,
           hz = all_nodes[el.node_loc[2]].z - all_nodes[el.node_loc[0]].z;
    GetG_Loc(all_nodes[el.node_loc[0]].r, all_materials[el.mater].lambda, hr, hz,
             A_loc); // A_loc = G_loc
    GetM_Loc(all_nodes[el.node_loc[0]].r, all_nodes[el.node_loc[0]].z,
             all_materials[el.mater].gamma_id, hr, hz, A_loc);
    Getb_Loc(all_nodes[el.node_loc[0]].r, all_nodes[el.node_loc[0]].z, hr, hz, b_loc,
             el.f_id);
    return 0;
}

int GeneratePortrait(std::vector<int> &ia, std::vector<int> &ja, int N, int Kel)
{
    ia.resize(N + 1);
    ja.resize(16 * Kel);
    std::vector<int> temp_list1(16 * Kel),
        temp_list2(16 * Kel);
    std::vector<int> listbeg(N);
    int listsize = 0;
    for (int i = 0; i < N; i++)
    {
        listbeg[i] = 0;
    }
    for (int ielem = 0; ielem < Kel; ielem++)
    {
        for (int i = 0; i < 4; i++)
        {
            int k = all_elems[ielem].node_loc[i];
            for (int j = i + 1; j < 4; j++)
            {
                int ind1 = k;
                int ind2 = all_elems[ielem].node_loc[j];
                if (ind2 < ind1)
                {
                    ind1 = ind2;
                    ind2 = k;
                }
                int iaddr = listbeg[ind2];
                if (iaddr == 0)
                {
                    listsize++;
                    listbeg[ind2] = listsize;
                    temp_list1[listsize] = ind1;
                    temp_list2[listsize] = 0;
                }
                else
                {
                    while (temp_list1[iaddr] < ind1 && temp_list2[iaddr] > 0)
                    {
                        iaddr = temp_list2[iaddr];
                    }
                    if (temp_list1[iaddr] > ind1)
                    {
                        listsize++;
                        temp_list1[listsize] = temp_list1[iaddr];
                        temp_list2[listsize] = temp_list2[iaddr];
                        temp_list1[iaddr] = ind1;
                        temp_list2[iaddr] = listsize;
                    }
                    else if (temp_list1[iaddr] < ind1)
                    {
                        listsize++;
                        temp_list2[iaddr] = listsize;
                        temp_list1[listsize] = ind1;
                        temp_list2[listsize] = 0;
                    }
                }
            }
        }
    }
    ia[0] = 0;
    for (int i = 0; i < N; i++)
    {
        ia[i + 1] = ia[i];
        int iaddr = listbeg[i];
        while (iaddr != 0)
        {
            ja[ia[i + 1]] = temp_list1[iaddr];
            ia[i + 1]++;
            iaddr = temp_list2[iaddr];
        }
    }
    ja.resize(ia[N]);
    return 0;
}

int AddLocal(std::vector<int> &ia, std::vector<int> &ja, std::vector<double> &di,
             std::vector<double> &al, std::vector<double> &au,
             std::vector<std::vector<double>> &A_loc,
             std::vector<double> &b, std::vector<double> &b_loc, int el_id)
    std::vector<int> L = all_elems[el_id].node_loc;
int k = all_elems[el_id].node_loc.size(); // размерность локальной матрицы
for (int i = 0; i < k; i++)
{
    di[L[i]] += A_loc[i][i];
    b[L[i]] += b_loc[i];
}
for (int i = 0; i < 4; i++)
{
    int temp = ia[L[i]];
    for (int j = 0; j < i; j++)
    {
        for (int k = temp; k < ia[L[i] + 1]; k++)
        {
            if (ja[k] == L[j])
            {
                al[k] += A_loc[i][j];
                au[k] += A_loc[j][i];
                k++;
            }
            break;
        }
    }
}
return 0;
}