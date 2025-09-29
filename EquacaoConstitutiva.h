#ifndef EQUACAO_CONSTITUTIVA
#define EQUACAO_CONSTITUTIVA

#include <stdlib.h>
#include <stdio.h>
#include "Malha.h"
#include "Derivadas.h"


typedef struct
{
    double **PsiVelho_xx, **PsiVelho_xy, **PsiVelho_yy, **PsiVelho_xt, **PsiVelho_yt, **PsiVelho_tt;
    double **PsiNovo_xx, **PsiNovo_xy, **PsiNovo_yy, **PsiNovo_xt, **PsiNovo_yt, **PsiNovo_tt;
    double **O_11, **O_12, **O_13, **O_21, **O_22, **O_23, **O_31, **O_32, **O_33;
    double **lambda1, **lambda2, **lambda3;
} VARIAVEIS_LOG_CONFORMATION;

double GeneralizedNewtonian_Viscosity(MALHA M, int n, double **U, double **V, double **MuNewtGen);

void GeneralizedNewtonian_Tensor(MALHA M, int n, double **U, double **V, double **MuNewtGen,
                                double **Txx, double **Txy, double **Tyy);

void EquacoesConstitutivas(double **LambdaVelho, double **MuVelho, double **NuVelho, double **LambdaNovo,  double **MuNovo, double **NuNovo, double **LambdaAux, double **NuAux,
                            double **TxxVelho, double **TxyVelho, double **TyyVelho, double **TxtVelho, double **TytVelho, double **TttVelho,
                            double **TxxNovo, double **TxyNovo, double **TyyNovo, double **TxtNovo, double **TytNovo, double **TttNovo,
                            double **EVPT_Structure, double **Norm_tau_dev,
                            VARIAVEIS_LOG_CONFORMATION *Variaveis_log_conf,
                            double **UVelho, double **VVelho, double **WVelho, double **UNovo, double **VNovo, double **WNovo, double **P,
                            double **TermoDominanteXX, double **TermoDominanteXY, double **TermoDominanteYY,
                            int n, MALHA M);

void EquacoesConstitutivasCartesiano(int i, int j,
                                     double **Txx_Velho, double **Txy_Velho, double **Tyy_Velho, double **Txt_Velho, double **Tyt_Velho, double **Ttt_Velho,
                                     double **Txx_Novo, double **Txy_Novo, double **Tyy_Novo, double **Txt_Novo, double **Tyt_Novo, double **Ttt_Novo,
                                     double **EVPT_Structure, double **Norm_tau_dev,
                                     double **U, double **V, double **W, double **P, int n, MALHA M);

void EquacoesConstitutivasNSF(int i, int j, double **LambdaVelho, double **MuVelho, double **NuVelho, double **LambdaNovo,  double **MuNovo, double **NuNovo, double **LambdaAux, double **NuAux,
                              double **UVelho, double **VVelho, double **UNovo, double **VNovo, int n, MALHA M);

void EquacoesConstitutivasLogConformation(int i, int j,
                                          double **PsiVelho_xx, double **PsiVelho_xy, double **PsiVelho_yy, double **PsiVelho_xt, double **PsiVelho_yt, double **PsiVelho_tt,
                                          double **PsiNovo_xx, double **PsiNovo_xy, double **PsiNovo_yy, double **PsiNovo_xt, double **PsiNovo_yt, double **PsiNovo_tt,
                                          double **O_11, double **O_12, double **O_13, double **O_21, double **O_22, double **O_23, double **O_31, double **O_32, double **O_33,
                                          double **Lambda1, double **Lambda2, double **Lambda3,
                                          double **U, double **V, double **W, int n, MALHA M);

void LOG_Termo_Wi(double o11, double o12, double o13, double o21, double o22, double o23, double o31, double o32, double o33,
                  double lambda1, double lambda2, double lambda3, double alfa, int dimensao, double **Result);

double EVPT_CalculaSigmaDev(MALHA M, int i, int j, 
                            double **U, double **V, double **P, 
                            double **Txx, double **Txy, double **Tyy, double **Ttt);

double DerivadaDUDX(double **U, int i, int j, MALHA M);

double DerivadaDUDY(double **U, int i, int j, double t, MALHA M);

double DerivadaDVDX(double **V, int i, int j, MALHA M);

double DerivadaDVDY(double **V, int i, int j, MALHA M);

double DerivadaDWDX(double **W, int i, int j, MALHA M);

double DerivadaDWDY(double **W, int i, int j, MALHA M);

double DerivadaDUDT(double **UVelho, double **UNovo, int i, int j, MALHA M);

double DerivadaDVDT(double **VVelho, double **VNovo, int i, int j, MALHA M);

double EVPT_ElasticModulus(MALHA M, double Lambda);

double EVPT_Viscosity(MALHA M, double Lambda);

int ProximoDaParede_Hibrido(int i, int j, MALHA M);



double CalculaDetMinA(MALHA M, double **Txx, double **Txy, double **Tyy);

void CalculaResiduosEqsConstitutivas(MALHA *M, double **UVelho, double **UNovo , double **VVelho, double **VNovo,
                          double **TxxVelho, double **TxxNovo, double **TxyVelho, double **TxyNovo, double **TyyVelho, double **TyyNovo,
                          double **LambdaVelho, double **LambdaNovo, double **MuVelho, double **MuNovo, double **NuVelho, double **NuNovo,
                          double **LambdaAux, double **NuAux,
                          double **PsiVelho_XX, double **PsiNovo_XX, double **PsiVelho_XY, double **PsiNovo_XY, double **PsiVelho_YY, double **PsiNovo_YY,
                          double **ResTxx, double **ResTxy, double **ResTyy,
                          double **ResLambda, double **ResMu, double **ResNu,
                          double **ResPsiXX, double **ResPsiXY, double **ResPsiYY);

void CalculaTermoDominanteEqsConstituvas(int i, int j,
                                         double **Txx_Velho, double **Txy_Velho, double **Tyy_Velho, double **Txt_Velho, double **Tyt_Velho, double **Ttt_Velho,
                                         double **Txx_Novo, double **Txy_Novo, double **Tyy_Novo, double **Txt_Novo, double **Tyt_Novo, double **Ttt_Novo,
                                         double **U, double **V, double **W,
                                         double **TermoDominanteXX, double **TermoDominanteXY, double **TermoDominanteYY,
                                         int n, MALHA M);

#endif // EQUACAO_CONSTITUTIVA
