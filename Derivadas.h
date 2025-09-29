#ifndef DERIVADAS_H
#define DERIVADAS_H

#include <stdio.h>
#include <stdlib.h>
#include "Malha.h"

#define QUAD(x)( ((x)*(x)) )

double Lx(double **U, int i, int j, double t, MALHA M);

double Ly(double **V, int i, int j, MALHA M);

double Cx(double **U, double **V, double **W, int i, int j, double t, MALHA M);

double Cy(double **U, double **V, int i, int j, MALHA M);

double Ct(double **U, double **V, double **W, int i, int j, MALHA M);

double GradPy(double **P, int i, int j, MALHA M);

double GradPx(double **P, int i, int j, MALHA M);

double DivT_u(double **Txx, double **Txy, double **Ttt, int i, int j, MALHA M);

double DivT_v(double **Txy, double **Tyy, int i, int j, MALHA M);

// The equivalent of the "divergent" part of the stress tensor term for generalized newtonian models
double DivT_u_NewtGen(MALHA M, double **U, double **V, double **ViscosityFunction, int i, int j);

double DivT_v_NewtGen(MALHA M, double **U, double **V, double **ViscosityFunction, int i, int j);

double DivT_t(double **Txt, double **Tyt, int i, int j, MALHA M);

double DivT_u_knownpressure(double **Txx, double **Txy, double **Ttt, int i, int j, MALHA M);

double Conv_Tensor(double **T, double **U, double **V, int i, int j, const char *Tipo, MALHA M);

double TensorParede(MALHA M, int i, int j, const char *Tipo, char Face);

double FuncaoCubista(double PhiR, double PhiU, double PhiD, double Xf, double Xr, double Xu, double Xd);

double FuncaoConvCentral(double PhiU, double PhiD, double Xf, double Xu, double Xd);

double FuncaoUpwind(double PhiU);

double Interpolacao(double h1, double h2, double Phi1, double Phi2);

double ValorInflowTensor(int i, int j, double **U, const char *Tipo, MALHA M, char Direcao);

void CalculaAutovalores_Jacobi(double a[4][4], int n, double d[4], double b[4], double z[4], double v[4][4], int *nrot);

#endif // DERIVADAS_H
