#ifndef METODO_PROJECAO_H
#define METODO_PROJECAO_H

#include <time.h>
#include <math.h>
#include <petscksp.h>
#include <petsctime.h>
#include "FuncoesAuxiliares.h"
#include "Malha.h"
#include "ContornoSuperficie.h"
#include "Derivadas.h"
#include "EquacaoMovimento.h"
#include "EquacaoPoisson.h"
#include "EquacaoConstitutiva.h"
#include "EVPT_Equacoes.h"
#include "Visualizacao.h"

//#define U(i, j) U[ (i)*Ny + j ]
//#define V(i, j) V[ (i)*(Ny+1) + j ]
//#define P(i, j) P[ (i)*(Ny) + j ]



typedef struct
{
    Mat A_u, A_v, A_w, A_p, A_s; //Matrizes de sistemas lineares
    Vec vetU, vetV, vetW, vetP, vetS; //Vetores independente de cada sistema
    Vec solU, solV, solW, solPsi, solS; //Solucoes dos sistemas lineares (U, V e Psi)
    KSP kspU, kspV, kspW, kspP, kspS; //Solvers do petsc

    double **UVelho, **VVelho, **WVelho, **PVelho;
    double **UMeio, **VMeio, **WMeio, **Psi;
    double **UNovo, **VNovo, **WNovo, **PNovo;
    double **TxxVelho, **TyyVelho, **TxyVelho, **TxtVelho, **TytVelho, **TttVelho; //Os tres ultimos sao apenas pro axissimetrico, coordenadas cilindricas
    double **TxxNovo, **TyyNovo, **TxyNovo, **TxtNovo, **TytNovo, **TttNovo; //Os tres ultimos sao apenas pro axissimetrico, coordenadas cilindricas
    double **lambdaVelho, **muVelho, **nuVelho;
    double **lambdaNovo, **muNovo, **nuNovo;
    double **evptStructureVelho, **evptStructureNovo;
    double **MuNewtGenVelho, **MuNewtGenNovo; //Apparent viscosity for generalized newtonian models (always 1 for purely Newtonian/viscoelastic)
    double **norm_tau_dev; // Norm of tau_dev for the EVP model. Used to determine if yielded/unyielded

    /// Some residuals to plot occasionally for debugging
    double **residualFaceU, **residualFaceV, **residualCell;

    /// === o metodo log conformation usa um monte de variavel auxiliar. coloquei tudo em uma struct
    /// === Estrutura definida no arquivo EquacaoConstitutiva.h
    VARIAVEIS_LOG_CONFORMATION variaveis_log_conf;

    //Auxiliares pra passsar pro termo convectivo da formulacao NSF
    double **lambdaAux, **nuAux;

    // Algumas coisas a mais apenas pra imprimir no arquivo VTK e visualizar
    double **termoDominanteXX, **termoDominanteXY, **termoDominanteYY;

    char sistemasDestruidos;

} SOLVER_PROJECAO;

double FuncaoParabolica(double PontoMin, double PontoMax, double ValorMax, double Ponto, double t, MALHA M);

void CondicaoInicialProjecao(SOLVER_PROJECAO Solver, MALHA M, double T);

void CondicaoContornoProjecao(double **U, double **V, double **W, MALHA M, double T);

PetscErrorCode InicializaEstruturasProjecao(SOLVER_PROJECAO *Solver, MALHA M);

PetscErrorCode ExecutaPassoTemporalProjecao(SOLVER_PROJECAO *Solver, MALHA *M, int n, int rank, LISTA_CELULAS_SURFACE *ListaSurf, double *Residuo);

void ProcuraParticulaNaN(MALHA *M, const char *Mensagem);

//Calcula UNovo e VNovo depois que UMeio, VMeio e Psi ja sao conhecidos
void CalculaVelocidadeNova(double **UNovo, double **VNovo, double **WNovo, double **UMeio, double **VMeio, double **WMeio, double **Psi, MALHA M);

void CalculaPressaoNova(double **PNovo, double **P, double **Psi, MALHA M);

PetscErrorCode AdvectaInterface(MALHA *M, double **UVelho, double **VVelho, double **UNovo, double **VNovo);

PetscErrorCode VerificaInterface(MALHA *M, double **UVelho, double **VVelho, double **UNovo, double **VNovo);

/* Funcoes auxiliares */
PetscErrorCode CopiaVecUPraMatriz(Vec Vetor, double **Matriz, MALHA M, int rank);

PetscErrorCode CopiaVecVPraMatriz(Vec Vetor, double **Matriz, MALHA M, int rank);

PetscErrorCode CopiaVecPsiPraMatriz(Vec Vetor, double **Matriz, MALHA M, int rank);

double NormaEntreSolucoesU(SOLVER_PROJECAO Solver, MALHA M);

void CopiaSolucaoNovaParaVelha(SOLVER_PROJECAO *Solver, MALHA M);

PetscErrorCode DestroiSistemasLineares(SOLVER_PROJECAO *Solver, MALHA M);

void FinalizaPrograma(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, SOLVER_PROJECAO *Solver);



/// Energias
double CalculaEnergiaCinetica(MALHA *M, double **U, double **V);

double CalculaEnergiaGrav(MALHA *M);

double CalculaEnergiaSuperficie(MALHA *M);

double CalculaEnergiaDissipativaVisc(MALHA *M, double **U, double **V, double EnergiaAnterior);

void CalculaTodasEnergias(MALHA *M, double **U, double **V, double **P, int PassoTemporal);

double CalculaFlowType(MALHA *M, double **U, double **V, int i, int j);

double CalculaMagnitudeVelocidade(double **U, double **V, int i, int j);

double CalculaEnergiaDissipativaViscCelula(MALHA *M, double **U, double **V, int i, int j, double EnergiaAnterior);

#endif // METODO_PROJECAO_H
