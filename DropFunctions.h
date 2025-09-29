#ifndef DROP_FUNCTIONS_H
#define DROP_FUNCTIONS_H

#include "Malha.h"
#include "InterfaceFrontTracking.h"
#include "ContornoSuperficie.h"
#include <math.h>

void DropDiameterAxisymmetric(MALHA M, double t);
void EstudoCodigo(MALHA M, LISTA_CELULAS_SURFACE *ListaCelulas, int PassoTemporal);
void ClassificaCelulasBlocoParaview(MALHA *M);

#endif
