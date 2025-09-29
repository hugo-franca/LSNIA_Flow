#ifndef INTER_FRONTTRACKING
#define INTER_FRONTTRACKING

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "FuncoesAuxiliares.h"

//Indica se a regiao eh fisica ou nao fisica
typedef enum StatusRegiao {INDEFINIDA, FISICA, NAO_FISICA, REMOVIDA} STATUS_REGIAO;

//Indica se a regiao eh externa ou interna (usado apenas no calculo de vetor normal, nao eh necessario na mudanca topologica
typedef enum TipoLocalRegiao {LOCAL_INDEFINIDO, EXTERNA, INTERNA} TIPO_LOCAL_REGIAO;

//Indica se uma curva esta saindo ou entrando em um node NH
typedef enum TipoCurva {TIPO_INDEFINIDO, ENTRANDO, SAINDO} TIPO_CURVA;

//Indica o formato-tipo da curva inserido na condicao inicial
typedef enum TipoRegiaoFree {TIPO_FREE_NENHUM, TIPO_FREE_ELIPSE, TIPO_FREE_IMAGEM, TIPO_FREE_RETANGULO, TIPO_FREE_INFLOW} TIPO_REGIAO_FREE;

//Representa um ponto na lista da curva (Nao eh um NODE, eh apenas um ponto qualquer no meio da curva)
typedef struct Ponto {
    double x, y; //Coordenadas no plano
    struct Ponto *prox, *ant; //Ponteiros pro proximo e anterior (Eh uma lista encadeada)
    char isNode; //Indica se eh um node ou nao

    //Se for node, indica se ele eh NH ou nao
    char isNH;

    //Se for NH vai usar essas 4 curvas no algoritmo de desembaracas
    void *curvas[4];
    char curvaEntrando[4]; //Se cada uma das 4 curvas ta entrando ou saindo deste node
    struct Ponto *seg1P1, *seg2P1; //Usados apenas no desembaracamento. Sao os pontos iniciais de cada segmento que embaracam

    //Vetor normal a interface neste ponto. Eh calculado e usado na funcao que gera a funcao marcadora
    double vetNormalX, vetNormalY;
    double comprimento; //Comprimento do segmento usado pra calcular a normal (soh eh usado no calculo da marker)

    int fixo; //Indica se este ponto vai ficar fico (sem ser advectado)

    int outflow; //Indica se este ponto passou pelo outflow (die-swell, por exemplo)

    //Usado apenas temporariamente no momento da adveccao
    double novoX, novoY;

    // Usado apenas pra um teste na simulacao FaradayWaves do Bernardo
    int wave;

    //Indica em qual celula da malha este ponto esta
    int i, j;
} PONTO;

typedef struct Regiao {
    int numero; //Apenas um identificador
    STATUS_REGIAO status;
    TIPO_LOCAL_REGIAO localRegiao;
} REGIAO;

typedef struct Curva {
    PONTO *pontoInicial;

    REGIAO *regiaoEsq; //Regiao da esquerda
    REGIAO *regiaoDir; //Regiao da direita

    double medioX, medioY; // Ponto medio desta curva, usado apenas para direcionar o vetor normal no novo algoritmo SURFACE_FULL
    double qtdPontosCurva; //Calculado e usado APENAS no algoritmo SURFACE_FULL

    struct Curva *proxCurva; //Ponteiro pra proxima curva da interface
    struct Curva *curvaAnt;
} CURVA;

//A interface eh uma lista de curvas
typedef struct Interface {
    CURVA *curvaInicial;
    int qtdCurvas;
} INTERFACE;

#define LOOP_CURVA(c) \
    for( struct { int primeiro_ponto; PONTO *pInicial; PONTO *ponto; } loop_stuff = { 1, (c)->pontoInicial, (c)->pontoInicial }; \
    (loop_stuff.ponto!=loop_stuff.pInicial) || (loop_stuff.primeiro_ponto); \
    loop_stuff.ponto=loop_stuff.ponto->prox, loop_stuff.primeiro_ponto=0 )

// Usado pra fazer a parametrizacao da interface por splines cubicas
typedef struct PolinomioGrau3
{
    double a;
    double b;
    double c;
    double d;
} POL3; //Estrutura para representar um polinomio de grau 3 por seus coeficientes

///FUNCOES DE ADVECCAO DE PARTICULAS
void AdvectaParticula(PONTO *Ponto, double U, double V, double dt, double *NovoX, double *NovoY);

void VetorNormal(CURVA *Curva, PONTO *Ponto, double *X, double *Y);

void VetorNormalMinimosQuadrados(INTERFACE *Inter, double Xc, double Yc, double Raio, double *ResultX, double *ResultY);

double CalculaCurvatura(INTERFACE *Inter, double Xc, double Yc, double RaioCurvatura, double RaioNormal, double NxCelulas, double NyCelulas);

double CalculaCurvaturaML(INTERFACE *Interf, double Xc, double Yc);

double AnguloVetor(double vx, double vy);

double Bezier(double *Beta, double T, int N);

double Dif_Bezier(double *Beta, double T, int N);

double Dif2_Bezier(double *Beta, double T, int N);

double CalculaCurvaturaPeloAnguloDeContato(INTERFACE *Interf, double Xc, double Yc, double RaioNormal);


double AreaTriangulo(PONTO *pI, PONTO *pJ, PONTO *pK);

int SOLVESYSTEM (int n, double **a, double *x);

char CriaEmbaracamentoForcado(INTERFACE *Interface, double DistanciaMinima);


///FUNCOES PARA LEITURA E CRIACAO DA INTERFACE A PARTIR DO ARQUIVO
CURVA *AdicionaCurvaElipse(INTERFACE *Interface, double XCentro, double YCentro, double RaioX, double RaioY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, double Angulo, int QtdSegmentos);

CURVA *AdicionaCurvaElipseAxissimetrico(INTERFACE *Interface, double XCentro, double YCentro, double RaioX, double RaioY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, double Angulo, int QtdSegmentos);

CURVA *AdicionaCurvaImagemAxissimetrico(INTERFACE *Interface, char *filename, REGIAO *RegiaoEsq, REGIAO *RegiaoDir);

CURVA *AdicionaCurvaSemiElipseAxissimetrico(INTERFACE *Interface, double XCentro, double YCentro, double RaioX, double RaioY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, double Angulo, int QtdSegmentos);

CURVA *AdicionaCurvaSemiElipse(INTERFACE *Interface, double XCentro, double YCentro, double RaioX, double RaioY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, int QtdSegmentos);

CURVA *AdicionaCurvaSemiElipseCoalescencia(INTERFACE *Interface, double R0, double ContactAngle, double TranslateX, double TranslateY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, int QtdSegmentos);

CURVA *AdicionaCurvaRetangulo(INTERFACE *Interface, double xMin, double xMax, double yMin, double yMax, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, double MaxDistancia);

CURVA *AdicionaCurvaFaradayWaves(INTERFACE *Interface, double Alfa, double Epsilon, REGIAO *RegiaoEsq, REGIAO *RegiaoDir);

CURVA *AdicionaCurvaFilamento(INTERFACE *Interface, double Rp, double R0, double Lz, REGIAO *RegiaoEsq, REGIAO *RegiaoDir);

CURVA *AdicionaCurvaCaberAxi(INTERFACE *Interface, double Rp, double R0, double Lz, REGIAO *RegiaoEsq, REGIAO *RegiaoDir);

CURVA *AdicionaCurvaSplinesUniforme(INTERFACE *Interface, CURVA *CurvaOriginal);

void Thomas(int NumEquacoes, double *Diag, double *DiagInf, double *DiagSup, double *Ind, double *Solucao);

void SplineCubica(double *X, double *Y, POL3 *S, int N);

void SplineCubicaNotAKnot(double *X, double *Y, POL3 *S, int N);

double ValorSpline(double t, double *T, POL3 *S, int N);

double ValorSplineDerivada(double t, double *T, POL3 *S, int N);

double ValorSplineDerivada2(double t, double *T, POL3 *S, int N);

CURVA *AdicionaCurvaMaziSpreading(INTERFACE *Interface, double h_inf, double R0, double max_r, REGIAO *RegiaoEsq, REGIAO *RegiaoDir);

CURVA *AdicionaCurvaBezierUniforme(INTERFACE *Interface, CURVA *CurvaOriginal);

CURVA *AdicionaCurvaBezierUniforme2(INTERFACE *Interface, CURVA *CurvaOriginal);

void AdicionaCurvaSessilePendingDrops(INTERFACE *Interface, double RaioCima, double RaioBaixo, double MeioX, double MeioY, double MinY, double MaxY, int BaseFixa);

double AlgoritmoBisseccaoBezier(double *BetaX, double *BetaY, int N, double T1, double T2, double X_Ant, double Y_Ant, double h);

double AlgoritmoBisseccaoSpline(POL3 *S_x, POL3 *S_y, double *T_vec, int N, double T1, double T2, double X_Ant, double Y_Ant, double h);

void AjustaExtremosFilamento(double *X, double *Y, double Lz);

CURVA *AdicionaCurvaInflow(INTERFACE *Interface, double Posicao, double ComprimentoInicial, double CoordMin, double CoordMax, double TamanhoBorda, char Direcao, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, int QtdSegmentos);

CURVA *AdicionaCurvaDoArquivo(INTERFACE *Interface, const char *NomeArquivo, REGIAO *RegiaoEsq, REGIAO *RegiaoDir);

double FuncaoFilamento(double z, double R0, double Lz);

double FuncaoZero_ColisaoGotas(double Theta, double C_x, double C_y, double b);

double Bisseccao_ColisaoGotas(double T1, double T2, double C_x, double C_y, double b);

///FUNCOES PRA UTILIZAR A ESTRUTURA DE DADOS DA INTERFACE
void AdicionaCurvaNaInterface(INTERFACE *Interface, CURVA *Curva);

REGIAO *CriaNovaRegiao(INTERFACE *Interface);

void SetaTodasRegioesComoIndefinidas(INTERFACE *Interface);

void InverteDirecaoCurva(CURVA *Curva);

void JuntaCurvasComNodesIguais(INTERFACE *Interface);

void FinalizaCurvasCiclicas(INTERFACE *Interface);

void AtualizaLocalDasRegioes(INTERFACE *Interface);

void MudaNodeDaCurva(INTERFACE *Interface, double X, double Y);

void AndaComONodeDaCurva(CURVA *C, PONTO *P, int Passos);

void RemovePontoDaCurva(CURVA *C, PONTO *P);

PONTO *UltimoPonto(PONTO *P);

char PontosIguais(PONTO *P1, PONTO *P2);

char MenorAntiHorario(PONTO *A, PONTO *B, PONTO *Centro);

char CurvaEntrando(CURVA *Curva, PONTO *Node);

char IsNodeNH(PONTO *Node, PONTO *NodesNH);

PONTO *CopiaPonto(PONTO *P);



///FUNCOES DE DESEMBARACAMENTO
void RealizaDesembaracamento(INTERFACE *Interface, PONTO *PontosNH);

char PossuiInterseccaoSegmentos(PONTO *Seg1Ponto1, PONTO *Seg1Ponto2, PONTO *Seg2Ponto1, PONTO *Seg2Ponto2, double *X, double *Y);

PONTO *PossuiEmbaracamentoInterface(INTERFACE *Interface);

void DivideSegmentosEmbaracados(INTERFACE *Interface, PONTO *NodesNH);

void OrdenaCurvasAntiHorario(PONTO *NodeH);

void DefineComponentesCurvas(INTERFACE *Interface, PONTO *NodesNH);

char PreencheUmaNovaRegiao(INTERFACE *Interface);

char ProcuraRegioesFisicas(INTERFACE *Interface, PONTO *NodesNH);

char ProcuraRegioesNaoFisicas(INTERFACE *Interface, PONTO *NodesNH);

void RemoveRegioesNaoFisicas(INTERFACE *Interface, PONTO *NodesNH);

void LiberaMemoriaCurva(CURVA *C);

#endif
