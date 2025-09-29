#ifndef EVPT_EQUACOES_H
#define EVPT_EQUACOES_H

#include "Malha.h"
#include "EquacaoConstitutiva.h"

double FlowCurve_MetodoDeNewton(MALHA M, double GammaAtual, double Stress);

void Evolucao_StructureParameter(MALHA M, double **U, double **V, double **P,
                                double **Txx, double **Txy, double **Tyy,
                                double **EVPT_Lambda_Velho, double **EVPT_Lambda_Novo);

#endif