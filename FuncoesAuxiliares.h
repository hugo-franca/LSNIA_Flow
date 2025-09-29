#ifndef FUNCOES_AUXILIARES_H
#define FUNCOES_AUXILIARES_H

#include <petscksp.h>
#include <stdlib.h>
#include <time.h>

typedef struct ItemPilhaGeral {
    void *valor;
    struct ItemPilhaGeral *prox;
}ITEM_PILHA_GERAL;

typedef struct PilhaGeral {
    ITEM_PILHA_GERAL *prim;
} PILHA_GERAL;

typedef struct ParametroOpcional {
    char nomeParametro[100];
    double valorParametro;
} PARAMETRO_OPCIONAL;

void PushPilha(PILHA_GERAL *P, void *Valor);

void *PopPilha(PILHA_GERAL *P);

void InsereParametroOpcional(PILHA_GERAL *P, char *StringParametro);

void **AlocaMatriz(int Linhas, int Colunas, size_t TamanhoElemento, void *ValorInicial);

void DesalocaMatriz(void **Matriz, int Linhas, int Colunas);

void TrocaPonteirosMatrizes(double ***Matriz1, double ***Matriz2);

void InicializaSolverKSP(KSP *Ksp, Mat A);

void ImprimeMatrizPraMatlab(Mat M, const char *NomeArquivo);

void ImprimeVetorPraMatlab(Vec V, const char *NomeArquivo);

void BroadcastMatrizMPI(double **Matriz, int Linhas, int Colunas);

double DeuProblema(const char* format, ...);

double FinalizouNormalmente(const char* format, ...);

void PrintDebug(const char* format, ...);

void Check_Read(int QtdLida, int QtdEsperada);

void IniciaTimer();

double FinalizaTimer();

#endif // FUNCOES_AUXILIARES_H
