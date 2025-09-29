#ifndef POISSON_H
#define POISSON_H

#include <petscksp.h>
#include "Malha.h"
#include "ContornoSuperficie.h"

PetscErrorCode InicializaMatrizPoisson(Mat A, int NumLinhas, MALHA M, int iStart, int iEnd);

PetscErrorCode MontaMatrizPoisson(Mat A, MALHA M, double **ViscosityNewtGen);

PetscErrorCode VetorIndependentePoisson(Vec VetVetor, double **U, double **V, double **PVelho, double **PNovo, double **Txx, double **Txy, double **Tyy, double **Psi, double **ViscosityNewtGen, MALHA M, LISTA_CELULAS_SURFACE *ListaSurf);

PetscErrorCode LinhaCondicaoContornoPoisson(Mat A, MALHA *M, double **ViscosityNewtGen, int i, int j, int Linha);

double LadoDireitoContornoPoisson(MALHA M, int i, int j, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, double **ViscosityNewtGen);

void DeterminaVetorNormalDeCadaCelula(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, double RaioParticulas);

void DeterminaCurvaturaDeCadaCelula(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, double RaioParticulas);

void DeterminaCurvaturaDeCadaCelulaVersao2(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, int QtdParticulas);

void DeterminaCurvaturaDeCadaCelulaVersaoSplines(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, int QtdParticulas);

void DeterminaCurvaturaDeCadaCelulaVersaoBezier(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, int QtdParticulas);

double CalculaCurvaturaDaParticulaSpline(MALHA *M, int QtdParticulas, PONTO *PontoEscolhido);

double CalculaCurvaturaBezier(MALHA *M, int QtdParticulas, PONTO *PontoEscolhido);

void AlocaLinhaCondicaoContornoPoisson(MALHA *M, int i, int j, int indiceVetor, int iStart, int iEnd, int *d_nnz, int *o_nnz, int *qtd);

void LaplacianoSupExplicito(MALHA *M, LISTA_CELULAS_SURFACE *CelulasSurf, double **U, double **V, double **P, double **Psi);

void LaplacianoSupImplicitoMatriz(Mat A, MALHA *M, int i, int j, int Linha);

void LaplacianoSupImplicitoRHS(Vec Vetor, MALHA *M, int i, int j, int Linha, double **U, double **V, double **P);

#endif // POISSON_H
