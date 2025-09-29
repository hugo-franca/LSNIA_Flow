#ifndef MALHA_H
#define MALHA_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "Python.h"

#include "InterfaceFrontTracking.h"
#include "FuncoesAuxiliares.h"
//#include "DropFunctions.h"

#define EPS_MALHA 1e-05
#define ZERO(x) ( fabs(x)<EPS_MALHA ? 1 : 0 )

///Macros para verificar se um ponto de velocidade possui um vizinho a direita, esquerda, cima ou baixo
#define DIREITA_V(i, j) ( (i==M.Nx-1) || ( (j==M.Ny || M.pontosU[(i)+1][(j)].tipo==EMPTY || M.pontosU[(i)+1][(j)].tipo==BOUNDARY) && (j==0 || M.pontosU[(i)+1][(j)-1].tipo==EMPTY || M.pontosU[(i)+1][(j)-1].tipo==BOUNDARY) ) )
#define ESQUERDA_V(i, j) ( (i==0) || ( (j==M.Ny || M.pontosU[(i)][(j)].tipo==EMPTY || M.pontosU[(i)][(j)].tipo==BOUNDARY) && (j==0 || M.pontosU[(i)][(j)-1].tipo==EMPTY || M.pontosU[(i)][(j)-1].tipo==BOUNDARY) ) )

#define BAIXO_U(i, j) ( (i==M.Nx || M.pontosV[i][j].tipo==EMPTY || M.pontosV[i][j].tipo==BOUNDARY) && (i==0 || M.pontosV[i-1][j].tipo==BOUNDARY || M.pontosV[i-1][j].tipo==EMPTY) )
#define CIMA_U(i, j) ( (i==M.Nx || M.pontosV[i][j+1].tipo==EMPTY || M.pontosV[i][j+1].tipo==BOUNDARY) && (i==0 || M.pontosV[i-1][j+1].tipo==BOUNDARY || M.pontosV[i-1][j+1].tipo==EMPTY) )

typedef enum TipoCelula {EMPTY, FULL, BOUNDARY, PERMA_EMPTY, SURFACE, INCOGNITA, MUDANCA_TOPOL, ERRO_CELULA} TIPO_CELULA;

typedef enum { NAO_BOUNDARY, INFLOW, NEUMANN, NOSLIP, SLIP, EIXO_AXISSIMETRIA, SIMETRIA, PERIODICIDADE, KNOWN_PRESSURE, ERRO_BOUNDARY } TIPO_BOUNDARY;

typedef enum { NAO_EXTREMO, ESQUERDA, DIREITA, CIMA, BAIXO } EXTREMO_DOMINIO;

typedef enum { ERRO_DIRECAO, DIRECAO_X, DIRECAO_Y, DIRECAO_XY } DIRECAO;

typedef enum { SEM_MOVIMENTO, MOVIMENTO_X, MOVIMENTO_Y, MOVIMENTO_THETA } MOVIMENTO_PAREDE;

typedef enum { NAO_INFLOW, PARABOLICO, PARABOLICO_SIMETRICO_INICIO, PARABOLICO_SIMETRICO_FIM, RETO, ERRO_PERFIL } PERFIL_INFLOW;

typedef enum { CARTESIANO, AXI_CILINDRICO } TIPO_COORDENADAS;

typedef enum { MODELO_VE, MODELO_VP, MODELO_EVP, MODELO_EVPT, MODELO_EVP_SARAMITO } TIPO_MODELO;

typedef enum { VISCO_CART, VISCO_NSF, VISCO_HIBRIDO, VISCO_LOG } FORMULACAO_VISCO;

typedef struct Posicao {
    int i, j;

    char incognita; //soh vai ser usado no sistema linear da superficie livre. ignorado no resto do codigo
    char descarta;
} POSICAO;

typedef struct InfoPonto {
    TIPO_CELULA tipo;
    TIPO_BOUNDARY tipoBoundary;
    MOVIMENTO_PAREDE movimento;
    EXTREMO_DOMINIO extremo;

    double valorDirichlet; //Usado como velocidade na condicao INFLOW ou pressao na condicao KNOWN_PRESSURE
    double valorDirichlet2; //Usado como gradiente da pressao na condicao GRAD_PRESSURE
    // double **tensorParedeVelho, **tensorParedeNovo;
} INFO_PONTO;

typedef struct Bloco {
    TIPO_CELULA tipoU, tipoV;
    TIPO_BOUNDARY tipoBoundaryU, tipoBoundaryV;
    DIRECAO direcaoNormal;
    MOVIMENTO_PAREDE movimentoParede;

    PERFIL_INFLOW perfilInflow;

    double valorDirichlet, valorDirichlet2;
    int iMin, iMax, jMin, jMax;

    //Extremos desta boundary
    double x1, y1, x2, y2;

    //Extremos e caracteristicas se for uma boundary free surface
    TIPO_REGIAO_FREE tipoRegiaoFree;
    double centroX, centroY, raioX, raioY; //Usados em caso de elipse
    double xMin, yMin, xMax, yMax; // Usados em caso de retangulo
    double freeVelX, freeVelY;

    int id;
    struct Bloco *prox; //Proximo bloco da lista
} BLOCO;

typedef struct BlocoParaview {
    int numFaces;
    double **matrizPontos;

    struct BlocoParaview *prox;
} BLOCO_PARAVIEW;

typedef struct PontoCelulaQuad {
    PONTO *ponto;

    struct PontoCelulaQuad *prox;
} PONTO_CELULA_QUAD;

typedef struct CelulaQuad {
    double xMin, xMax, yMin, yMax;
    PONTO_CELULA_QUAD *ponto_celula;
    int qtdPontos;

    struct CelulaQuad *cima_esq, *cima_dir, *baixo_esq, *baixo_dir;
    struct CelulaQuad *pai;
} CELULA_QUAD;

typedef struct Malha {
    INTERFACE interface; //Interface front-tracking

    int Nx, Ny, Nt; //Quantidade de celulas e passos temporais
    double *dx, *dy; //Espacamentos, malha com stretching
    double minDxDy; //O menor espacamento entre todos os dx e dy
    double dt;
    TIPO_CELULA **celulas;

    /// === Pra usar na malha quad_tree
    CELULA_QUAD *quad_celula;
    int qtdCelulas;

    TIPO_CELULA **celulasParede; //Gambiarra que eu adicionei depois pra evitar particulas de entrarem em um bloco que deveria ser de parede

    TIPO_COORDENADAS tipoCoord;
    double *r; //Valores das coordenadas r da malha pra usar no axissimetrico. Vou deixar tudo igual a 1 no cartesiano
    double *x, *y; //Coordenadas x e y da malha. Vou precisar disso pros convectivos no caso nao-uniforme

    INFO_PONTO **pontosU;
    INFO_PONTO **pontosV;
    INFO_PONTO **pontosW;
    INFO_PONTO **pontosP;

    BLOCO *primBloco;
    BLOCO *blocoLivre; //Blocos de regiao free-surface. Usado apenas na condicao inicial, mais nada

    double xMin, xMax;
    double yMin, yMax;
    double tMax;
    double Re, Fr; //Reynolds e Froude
    double We, beta; //Parametros viscoelasticos Oldroid-B, PTT e Geisekus
    double epsilon; //Parametro viscoelastico PTT
    double alpha; //Parametro viscoelastico Geisekus
    double weber; //Parametro da tensao superficial
    double gravX, gravY; //Gravidade adimensional (geralmente deve ser [0, -1])
    double tolNSF, tolNSFconversao; //Tolerancia usada apenas na formulacao NSF e HIBRIDA pra parte viscoelastica
    int grauHibrido; //Grau da formulacao hibrida NSF-CSF. Valor eh ignorado se nao usar hibrido
    double velInicialX, velInicialY; //Velocidade como condicao inicial
    FORMULACAO_VISCO formulacaoVisco; //Formulacao usada na parte viscoelastica
    
    // Parametros do modelo EVPT
    double eta_inf; // viscosity for the unstructured material
    double eta_0; // viscosity for the fully structured material
    double yield_stress; // yield stress...
    double K; //consistency parameter from that Herschel-Buckley flow curve
    double power_n; //Power law exponente from the Herschel-Buckley flow curve
    double time_eq; //Thixotropy equilibrium time
    double evpt_m; // Parameter that goes in the elastic modulus function (G(lambda))

    // Parametros do modelo EVP-Saramito
    double Bi;

    // Duas opcoes: MODELO_VE (viscoelastico) ou MODELO_EVPT (elasto-visco-plastico-thixotropico)
    TIPO_MODELO tipo_modelo;

    // Um parametro de teste para usar na simulacao colisao de duas gotas
    double merge_time;

    int tensaoSuperficial; //Indica se a simulacao vai usar tensao superficial ou nao
    int freqTSUR; //A frequencia com a qual vai usar o tsur (ZERO = nao vai usar tsur)
    double **curvaturaCelulas; //Curvatura de cada celula i, j
    PONTO** *pontoCelula; // Matriz de pontos, um ponto pra cada celula. Este eh o ponto da interface escolhido para cada celula
    double **vetNx, **vetNy; //Vetor normal de cada celula usado apenas pra calcular a curvatura

    POSICAO *indicesU, *indicesV, *indicesW, *indicesP;
    int **matIndicesU, **matIndicesV, **matIndicesW, **matIndicesP;
    int qtdIncognitasU, qtdIncognitasV, qtdIncognitasW, qtdIncognitasP; //Quantidade de incognitas

    //Usados apenas no sistema linear da superficie livre
    POSICAO *indicesInterface;
    char **pontosContinuidade, **pontosSupLivre;

    //Nome escolhido pelo usuario para os arquivos vtk de visualizacao
    char nomeArquivo[200];
    int intervaloVTK_surf, intervaloVTK_prop, intervaloEstado, intervaloPrint;

    // Passo temporal atual pra usar em alguns algoritmos que dependem do tempo...
    int passoTemporal;

    //Parametros opcionais passados no final do arquivo modelo
    PILHA_GERAL parametrosOpcionais;

    int **malhaRecursoes;

    int tensaoSup_full;

    //double **PontosBlocoParaview;
    double x1, x2, x3, x4;
    double y1, y2, y3, y4;
} MALHA;


///Essas estruturas de pilha sao apenas usadas no FloodFill que inicializa o dominio
typedef struct ItemPilha {
    int i, j;
    struct ItemPilha *prox;
} ITEM_PILHA;

typedef struct Pilha {
    ITEM_PILHA *topo;
} PILHA;


///Essas estruturas sao usadas apenas no algoritmo que faz a mudanca topologica forcada
typedef struct PontoMudancaTopol
{
    PONTO *ponto;
    double direcaoMovX, direcaoMovY;
    struct PontoMudancaTopol *prox, *ant;
} PONTO_MUDANCA_TOPOL;

typedef struct ListaMudancaTopol
{
    PONTO_MUDANCA_TOPOL *pontoInicial;
    PONTO_MUDANCA_TOPOL *ultimoPonto;
} LISTA_MUDANCA_TOPOL;

MALHA LerModeloDoArquivo(const char *NomeArquivo);

void InsereBloco(MALHA *M, int id, char *TipoCelula, char *DirecaoNormal, char *TipoBoundary, double ValorDirichlet, double ValorDirichlet2, char *PerfilInflow,
                 char *Movimento, int iMin, int iMax, int jMin, int jMax,
                 double X1, double Y1, double X2, double Y2);

void CriaMalha(MALHA *Malha, const char *TipoStretchingX, int QtdVertX, double *VerticesX, int *Nx, double *ParametrosX, const char *TipoStretchingY, int QtdVertY, double *VerticesY, int *Ny, double *ParametrosY);

void CriaMalhaUniforme(MALHA *Malha, int Nx, int Ny);

void CriaMalhaStretchingTrygvasson(MALHA *Malha, int QtdVertX, double *VerticesX, int *Nx, double *ParametrosX, int QtdVertY, double *VerticesY, int *Ny, double *ParametrosY);

void CriaMalhaStretchingGeometrico(MALHA *Malha, int QtdVertX, double *VerticesX, int *Nx, double *ParametrosX, int QtdVertY, double *VerticesY, int *Ny, double *ParametrosY);


void InicializaMalhaComEmpty(MALHA *M);

void InicializaIndiceIncognitas(MALHA *M);

void InicializaExtremosContorno(MALHA *M);

void LeituraMalhaStretching(FILE *Arq, char *Tipo, int *QtdVertices, double **Vertices, int **Espacamentos, double **Parametros);

TIPO_CELULA ReconheceTipoDaStringCelula(char *Cadeia);

TIPO_BOUNDARY ReconheceTipoDaStringBoundary(char *Cadeia);

DIRECAO ReconheceTipoDaStringDirecao(char *Cadeia);

MOVIMENTO_PAREDE ReconheceTipoDaStringMovimento(char *Cadeia);

PERFIL_INFLOW ReconheceTipoDaStringPerfilInflow(char *Cadeia);

void PreencheRegiao(MALHA M, int NumFaces, int *Faces);

void PreencheRegiaoParede(MALHA M, int NumFaces, int *Faces);

void PreencheRegiaoPoligono(MALHA M, CURVA *C);

BLOCO *BuscaBloco(MALHA M, int id);

void InsereBlocoFreeSurface(MALHA *M, TIPO_REGIAO_FREE TipoBloco, double VelX, double VelY, double xMin, double yMin, double xMax, double yMax, double RaioX, double RaioY, double CentroX, double CentroY);

BLOCO *BuscaBlocoFreeSurface(MALHA *M, double X, double Y);

void AlgoritmoFloodFillU(MALHA M, int i, int j, TIPO_CELULA Tipo);

void AlgoritmoFloodFillV(MALHA M, int i, int j, TIPO_CELULA Tipo);

void AlgoritmoFloodFillP(MALHA M, int i, int j, TIPO_CELULA Tipo);

void AlgoritmoFloodFillUIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo);

void AlgoritmoFloodFillVIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo);

void AlgoritmoFloodFillPIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo);

void AlgoritmoFloodFillCelulaParedeIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo);

void AlgoritmoFloodFillWIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo);

void AlgoritmoFloodFillCelulaIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo);

double AlgoritmoBisseccao(int Grau, double *Polinomio, double X1, double X2);

double ValorPolinomio(double X, int Grau, double *Pol);

void EncontraCelulaDaParticula(MALHA *M, double X, double Y, int *iResult, int *jResult);


/// Funcoes para forcar o embaracamento da interface dependendo da classificacao das celulas na malha
int EmbaracamentoForcado_Raycast(MALHA *M, double Distancia);

void EmbaracamentoForcado_DuasListas(MALHA *M, double Distancia, PONTO *p1C1, PONTO *p1C2, CURVA *c1, CURVA *c2);

int RealizaEmbaracamentoForcado(MALHA *M);

int RealizaMudancaTopologica_Uniao(MALHA *M);

int RealizaEmbaracamentoForcado_Quebra2(MALHA *M, double DistanciaMax, double DistanciaMin, PONTO **pontoRemover, int (*CriterioOpcional)(MALHA*, PONTO*, PONTO*, int) );

int CriterioOpcional_Quebra(MALHA *M, PONTO *p1, PONTO *p2, int Opcao);

void RemoveCurvaDestePonto(MALHA *M, PONTO *pRemover);

int RealizaEmbaracamentoForcado_Quebra(MALHA *M);

void AtualizaMinimoMaximo(int i, int j, int *MinI, int *MinJ, int *MaxI, int *MaxJ);

void EmbaracamentoCasoDuasListas(MALHA *M, LISTA_MUDANCA_TOPOL *L);

void EmbaracamentoCasoUmaLista(MALHA *M, LISTA_MUDANCA_TOPOL Lista);

void RemoveCelulasIndesejadasMudancaTopologicaPasso1(MALHA *M, int *QtdCelulas);

void RemoveCelulasIndesejadasMudancaTopologicaPasso2(MALHA *M, int *QtdCelulas);

void AdicionaNaListaDeMudancaTopol(LISTA_MUDANCA_TOPOL *L, PONTO *Ponto);

void RemoveParticulasFULL(MALHA *M);

void RemoveParticulasNaSimetriaSURFACE(MALHA *M);

void RemoveCurvasMuitoPequenas(MALHA *M, int MinimoParticulas);

void CheckUpInterface(MALHA *M, const char *NomeTeste);

void SuavizaOndulacoes(MALHA *M);

void SuavizaOndulacoesAoRedorDoPonto(MALHA *M, double X, double Y, int qtdPontos);

/// FUNCOES USADAS PRA CRIAR A FUNCAO QUE DA O VALOR PRA UMA CONDICAO DE CONTORNO NOSLIP EM MOVIMENTO
double FuncaoValorNoSlipHorizontal(MALHA M, int i, int j, double X, double T, double **U);

double FuncaoValorNoSlipVertical(MALHA M, int i, int j, double X, double T, double **V);


void AdicionaERemoveParticulas(MALHA *M, double TamanhoMaximo, int PermiteCelulasEmpty);

void DivideSegmentoSeNecessario(MALHA *M, PONTO *Ponto, double TamanhoMaximo, int PermiteCelulasEmpty);

void DivideSegmentoSpline(MALHA *M, PONTO *PontoEscolhido, double *NovoX, double *NovoY);

void EncontraCelulaDeTodasAsParticulas(MALHA *M);

void SuavizaInterfaceTSUR(MALHA *M);

void TSUR(MALHA *M, PONTO *P1);



///FUNCOES DE VISUALIZACAO
void DesenhaMalhaGnuplot(MALHA Malha);

void DesenhaMalhaVTK(MALHA Malha, int n);


///Funcoes auxiliares pra usar a pilha do flood fill
void PilhaPush(PILHA *P, int i, int j);

void PilhaPop(PILHA *P, int *i, int *j);

int SubstituiString(char *Str, const char *Orig, char *Sub);

void StringMaiuscula(char *String);

void ArquivoMaiusculoTemporario(FILE *ArqOriginal, FILE *ArqNovo);

void ArmazenaQuinasBlocoParaview(MALHA *M, BLOCO_PARAVIEW *PrimBloco, int NumFaces);

void AdicionaBlocoParaview(MALHA *M, BLOCO_PARAVIEW **PrimBloco, int NumFaces, int *Faces);

void EscreveArquivoBlocosParaview(BLOCO_PARAVIEW *PrimBloco, MALHA M);

double CalculaVelocidadeDaGota(MALHA *M, double **U, double **V, int Passo);



/// == Quad tree
void CriaMalhaQuadTree(MALHA *M);

void DivideCelulaQuad(MALHA *M, CELULA_QUAD *C);

CELULA_QUAD *NovaCelulaQuad(CELULA_QUAD *Pai, double xMin, double xMax, double yMin, double yMax);

void DesenhaMalhaQuadVTK(MALHA M, int n);

void ImprimePontosQuadTree(FILE *Arq, CELULA_QUAD *C);

void MovePontoQuadTree(CELULA_QUAD *C, PONTO_CELULA_QUAD *P);

void AdicionaPontoQuadTree(MALHA *M, CELULA_QUAD *C, PONTO *PontoInt, int MaxPontos);

double BuscaParametro(MALHA *M, const char *NomeParametro);

float calcAreaFluido(MALHA malha);

float calcSpread(MALHA malha);

float calcAltura(MALHA malha);

void fprintfVetorFormatado(FILE *arq, float *vet, int tam);

void pontoEsquerda(MALHA malha, int *i, int *j);

void pontoDireita(MALHA malha, int *i, int *j);

void pontoMeio(MALHA malha, int *i, int *j);

#endif // MALHA_H
