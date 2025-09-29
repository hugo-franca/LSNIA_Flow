#ifndef EQ_MOVIMENTO_H
#define EQ_MOVIMENTO_H

#include <petscksp.h>
#include "Malha.h"
#include "Derivadas.h"
#include "ContornoSuperficie.h"

PetscErrorCode InicializaMatrizCoeficientesU(Mat A, int NumLinhas, MALHA M, int iStart, int iEnd);

PetscErrorCode MontaMatrizCoeficientesU(Mat A, MALHA M, double **ViscosityNewtGen);

PetscErrorCode VetorIndependenteU(Vec VetVetor,
                                  double **UVelho, double **VVelho, double **WVelho,
                                  double **UNovo, double **VNovo, double **WNovo,
                                  double **PVelho,
                                  double **Txx, double **Txy, double **Ttt,
                                  double **ViscosityFunction,
                                  int n, MALHA M);

PetscErrorCode InicializaMatrizCoeficientesV(Mat A, int NumLinhas, MALHA M, int iStart, int iEnd);

PetscErrorCode MontaMatrizCoeficientesV(Mat A, MALHA M, double **ViscosityNewtGen);

PetscErrorCode VetorIndependenteV(Vec VetVetor, 
                                    double **UVelho, double **VVelho, double **UNovo, double **VNovo,
                                    double **PVelho, double **Txy, double **Tyy, 
                                    double **ViscosityFunction,
                                    MALHA M);

PetscErrorCode InicializaMatrizCoeficientesW(Mat A, int NumLinhas, MALHA M, int iStart, int iEnd);

PetscErrorCode MontaMatrizCoeficientesW(Mat A, MALHA M);

PetscErrorCode VetorIndependenteW(Vec VetVetor, double **UVelho, double **VVelho, double **WVelho, double **PVelho,
                                  double **Txt, double **Tyt, MALHA M);

void CalculaNormalCurvaturaEquacaoMovimento(MALHA *M, double X, double Y, double *ResultNx, double *ResultNy, double *ResultCurv);

void ConfereDirecaoVetorNormal(MALHA *M, PONTO *P, CURVA *C, double *NormalX, double *NormalY);

void CalculaDirecaoExternaSurfaceFull(MALHA *M, int i, int j, double *DirecaoX, double *DirecaoY);

void CalculaResiduoMomentoU(MALHA *M, double **UVelho, double **UMeio, double **UNovo,
                            double **VVelho, double **VMeio, double **VNovo,
                            double **PVelho, double **PNovo,
                            double **TxxVelho, double **TxxNovo, double **TxyVelho, double **TxyNovo, double **TttVelho, double **TttNovo,
                            double **Residuo);



#endif
