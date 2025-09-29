#ifndef VISUALIZACAO_H
#define VISUALIZACAO_H

#include "Malha.h"
#include "InterfaceFrontTracking.h"

void ImprimeArquivosUVPGnuplot(double **U, double **V, double **W, double **P, MALHA Malha, int n);

void ImprimeArquivoVTK(double **U, double **V, double **W, double **P,
                       double **Txx, double **Txy, double **Tyy, double **Txt, double **Tyt, double **Ttt,
                       double **Lambda, double **Mu, double **Nu,
                       double **EVPT_Structure, double **Norm_tau_dev,
                       double **MuNewtGen,
                       double **TermoDominanteXX, double **TermoDominanteXY, double **TermoDominanteYY,
                       double **ResidualFaceX, double **ResidualFaceY, double **ResidualCell,
                       MALHA Malha, int n);

void ImprimeInterfaceVTK(MALHA Malha, int n);

void ImprimeArquivoBinario(double **U, double **V, MALHA Malha, int n);

void ImprimeCorte(double **U, double **V, double **W, int Indice, const char *Direcao, const char *Variavel, const char *NomeArquivo, MALHA M);

void PlotU(MALHA Malha, int n);

void PlotV(MALHA Malha, int n);

void PlotP(MALHA Malha, int n);

void PlotW(MALHA Malha, int n);



///Superficie livre
void ImprimePontosCurva(CURVA Curva, FILE *Arquivo);

void ImprimeCurvasEPontos(INTERFACE Interface);

void PlotInterfaceGnuplot(MALHA Malha, char DesenhaMalha, char DesenhaNormais, int n);

void ImprimeArquivoColetaVTK(float **flowType, float **magnitudeVelocidade, float **energiaDissipativaCelulas, MALHA Malha, int n);

#endif // VISUALIZACAO_H
