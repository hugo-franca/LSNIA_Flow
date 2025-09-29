#include "DropFunctions.h"

void DropDiameterAxisymmetric(MALHA M, double t){

	double diameterX=0.0, xMax=0.0, xMin=0.0, yMax=0.0, yMin=0.0, diameterY=0.0;
	CURVA *curva;
	PONTO *ponto;

	//FILE *arqInterface;
	FILE *arqDiameter;

	//arqInterface = fopen("./ArquivosPlot/Interface/DropInterface.txt", "wt");
	for (curva = M.interface.curvaInicial; curva != NULL; curva = curva->proxCurva) {
	//for (curva = M.interface.curvaInicial; curva != NULL; curva = curva->proxCurva) {
        	ponto = curva->pontoInicial;
        	do {
			//printf("x=%f\n", ponto->x);
			//Diametro x:
			if (ponto->x > xMax){
				xMax = ponto->x;
			}else if (ponto->x < xMin){
				xMin = ponto->x;
			}

			//Menor ponto em y
			if (ponto->y < yMin){
				yMin = ponto->y;
			}else if (ponto->y > yMax){
				yMax = ponto->y;
			}
			 ponto = ponto->prox;
		} while (ponto != curva->pontoInicial);
		//fprintf(arqInterface, "%lf \t %lf\n", ponto->x, ponto->y);		
	}
	//fclose(arqInterface);

	//printf("xMax=%f \t xMin=%f\n", xMax, xMin);getchar();
	diameterY = yMax - yMin;
	diameterX = 2.0*(xMax - xMin);
	arqDiameter = fopen("./ArquivosPlot/Interface/DropDiameter_Tan_Nor.txt", "a");
	fprintf(arqDiameter, "%lf \t %lf \t %lf \n", t, diameterX, diameterY);
	fclose(arqDiameter);
}

void EstudoCodigo(MALHA M, LISTA_CELULAS_SURFACE *ListaCelulas, int PassoTemporal){

	CURVA *curva;
	PONTO *ponto;
	int i, j, cont=0;

	/*
	//Percorrer a curva que define a superficie livre
	for( curva=M.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
		cont++;
		//printf("cont=%d\n", cont);getchar();
		ponto = curva->pontoInicial;
		do{
			
			EncontraCelulaDaParticula(&M, ponto->x, ponto->y, &i, &j);
			printf("ponto->x=%f, ponto->y=%f\n", ponto->x, ponto->y);
			printf("i=%d, j=%d\n", i, j);
			printf("M.x[%d]=%f, M.y[%d]=%f\n", i, M.x[i], j, M.y[j]);getchar();


			ponto=ponto->prox;
		}while(ponto != curva->pontoInicial);

	}
	*/
	/*
	for (i=0; i<=M.Nx; i++){
		for (j=0; j<=M.Ny; j++){
			printf("M.celulas[i][j] = %d\n", M.celulas[i][j]);
			printf("M.pontosU[i][j].tipo=%d, M.pontosU[i+1][j].tipo=%d, M.pontosV[i][j].tipo=%d, M.pontosV[i][j+1].tipo=%d\n", M.pontosU[i][j].tipo, M.pontosU[i+1][j].tipo, M.pontosV[i][j].tipo, M.pontosV[i][j+1].tipo);
			printf("M.x[%d]=%f, M.y[%d]=%f\n", i, M.x[i], j, M.y[j]);getchar();
		}
	}
	*/

}


void ClassificaCelulasBlocoParaview(MALHA *M){


}


