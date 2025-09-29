#ifndef CONTORNO_SUPERFICIE_H
#define CONTORNO_SUPERFICIE_H

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <petscksp.h>

#include "FuncoesAuxiliares.h"
#include "Malha.h"
#include "InterfaceFrontTracking.h"
#include "Derivadas.h"

#define CELULA_FLUIDO(i, j) (M->celulas[i][j]==FULL || M->celulas[i][j]==SURFACE)

#define LADO_VERTICAL_DIREITA(i, j) (LADO_VERTICAL_DIREITA_func(i, j, M))

#define LADO_VERTICAL_ESQUERDA(i, j) (LADO_VERTICAL_ESQUERDA_func(i, j, M))

#define LADO_HORIZONTAL_BAIXO(i, j) (LADO_HORIZONTAL_BAIXO_func(i, j, M))

#define LADO_HORIZONTAL_CIMA(i, j) (LADO_HORIZONTAL_CIMA_func(i, j, M))

#define QUINA_BAIXO_DIREITA(i, j) (QUINA_BAIXO_DIREITA_func(i, j, M))

#define QUINA_BAIXO_ESQUERDA(i, j) (QUINA_BAIXO_ESQUERDA_func(i, j, M))

#define QUINA_CIMA_DIREITA(i, j) (QUINA_CIMA_DIREITA_func(i, j, M))

#define QUINA_CIMA_ESQUERDA(i, j) (QUINA_CIMA_ESQUERDA_func(i, j, M))

#define DEGENERADO_CIMA_BAIXO(i, j) (DEGENERADO_CIMA_BAIXO_func(i, j, M))

#define DEGENERADO_ESQ_DIR(i, j) (DEGENERADO_ESQ_DIR_func(i, j, M))

#define DEGENERADO_ESQ_CIMA_BAIXO(i, j) (DEGENERADO_ESQ_CIMA_BAIXO_func(i, j, M))

typedef enum {SEM_EXTRAPOLACAO, ESQ_CIMA, ESQ_BAIXO, DIR_CIMA, DIR_BAIXO} DIRECAO_EXTRAPOLACAO;

typedef struct CelulaSurface {
    int i, j;
    int passoCelula;
    int novaCelulaSurface; // Indica se essa celula acabou de sair de EMPTY e virar SURFACE neste passo
    double **matrizCurvatura;
    struct CelulaSurface *prox, *ant;

    // Pra usar na funcao CalculaCurvaturaVersao2
    PONTO *pontoEscolhido;
    double menorDistancia;

} CELULA_SURFACE;

typedef struct {
    CELULA_SURFACE *prim;
    CELULA_SURFACE *ult;
} LISTA_CELULAS_SURFACE;

LISTA_CELULAS_SURFACE ClassificaCelulasEPontos(MALHA *M, CURVA *C, double **U, double **V, double **P);

void ReclassificaCelulasEPontos(MALHA *M, LISTA_CELULAS_SURFACE *ListaCelulas, int PassoTemporal);

void InicializaNovasPropriedadesSURFACE(MALHA *M, LISTA_CELULAS_SURFACE *ListaCelulas, double **EVPT_Structure);

void TratandoCasosDegenerados(MALHA *M, LISTA_CELULAS_SURFACE *ListaCelulas, int PassoTemporal);

int EncontraGotasPequenas_Surface(MALHA *M, int i, int j, char **CelulasVisitadas);

void RemoveGotaPequena_Surface(MALHA *M, int i, int j, char **CelulasVisitadas, LISTA_CELULAS_SURFACE *ListaSurf);

PetscErrorCode CondicaoContornoTangencial(int n, MALHA *M, LISTA_CELULAS_SURFACE listaCelulas, double **U, double **V, double **UNovo, double **VNovo, double **PNovo, double **Txx, double **Txy, double **Tyy, KSP *Ksp, Mat *A, Vec *Vet, Vec *Sol);

void NovaIncognitaSuperficie(MALHA *M, int i, int j, char Incognita, int *QtdIncognitas, char Descarta, int *TamanhoBuffer);

int NovaEquacaoSuperficie(MALHA M, Mat A, Vec Vet, double **U, double **V, double **Txx, double **Txy, double **Tyy, int i, int j, char TipoEquacao, char DirecaoExtrapolacao, int *QtdEquacoes, int LinhaMatriz);

double ValorPressaoCelulaSurface(MALHA M, int i, int j, char TipoEquacao, char Direcao, double **UVelho, double **VVelho, double **UNovo, double **VNovo);

void VetorNormalForaDaInterface(MALHA M, double PontoX, double PontoY, double *NormalX, double *NormalY);

void VelocidadeEmUmaParticula(MALHA M, PONTO *P, double **U, double **V, double *VelU, double *VelV);

double FuncaoDistribuicao(double X, double Y, PONTO *P, double dx, double dy);

void EncontraCelulaSeedFloodFill(MALHA M, LISTA_CELULAS_SURFACE L, int *iSeed, int *jSeed);

void EncontraVelPressaoSeedFloodFill(MALHA M, int *iSeedU, int *jSeedU, int *iSeedV, int *jSeedV, int *iSeedP, int *jSeedP);

//Manutencao da lista de celulas
int AdicionaCelulaSurface(MALHA *M, LISTA_CELULAS_SURFACE *L, int i, int j, int Passo);

//void ModificaCelulaSurface(LISTA_CELULAS_SURFACE *L, CELULA_SURFACE *Celula, int i, int j);

void RemoveCelulaSurface(LISTA_CELULAS_SURFACE *L, CELULA_SURFACE *Celula, MALHA *M, TIPO_CELULA Tipo);

CELULA_SURFACE *EncontraCelulaSurface(LISTA_CELULAS_SURFACE *L, int i, int j);

LISTA_CELULAS_SURFACE UniaoListasCelulaSurface(MALHA *M, LISTA_CELULAS_SURFACE *L, int QtdListas);

void DestroiListaCelulasSurface(LISTA_CELULAS_SURFACE *L);




//Funcoes para detectar casos degenerados
int CasoDegeneradoCima(MALHA M, int i, int j);

void DeletaCurvasFULL(MALHA *M);

void ArrumaMudancaTopologicaAxissimetrica(MALHA *M);

int DeletaCurvasPequenasCompletamenteSURFACE(MALHA *M);

int RecursaoDeletaCurvasSURFACE(MALHA *M, CURVA *curva, int i, int j);

int is_point_in_path(double x, double y, CURVA *curva);




void SalvarEstado(MALHA M, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, LISTA_CELULAS_SURFACE ListaSurf, int n);

void SalvarEstadoBinario(MALHA M, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, double **Lambda, double **Mu, double **Nu, LISTA_CELULAS_SURFACE ListaSurf, int n);

void CarregarEstado(MALHA *M, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, LISTA_CELULAS_SURFACE *ListaSurf, int *n);

void CarregarEstadoBinario(MALHA *M, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, double **Lambda, double **Mu, double **Nu, LISTA_CELULAS_SURFACE *ListaSurf, int *n);


///CRITERIOS PARA CLASSIFICACAO DE CELULAS
int LADO_VERTICAL_DIREITA_func(int i, int j, MALHA *M);

int LADO_VERTICAL_ESQUERDA_func(int i, int j, MALHA *M);

int LADO_HORIZONTAL_BAIXO_func(int i, int j, MALHA *M);

int LADO_HORIZONTAL_CIMA_func(int i, int j, MALHA *M);

int QUINA_BAIXO_DIREITA_func(int i, int j, MALHA *M);

int QUINA_BAIXO_ESQUERDA_func(int i, int j, MALHA *M);

int QUINA_CIMA_DIREITA_func(int i, int j, MALHA *M);

int QUINA_CIMA_ESQUERDA_func(int i, int j, MALHA *M);

int DEGENERADO_CIMA_BAIXO_func(int i, int j, MALHA *M);

int DEGENERADO_ESQ_DIR_func(int i, int j, MALHA *M);

int DEGENERADO_ESQ_CIMA_BAIXO_func(int i, int j, MALHA *M);

#endif // CONTORNO_SUPERFICIE_H
