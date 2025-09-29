#include "InterfaceFrontTracking.h"


///FUNCOES USADAS PRA FAZER ADVECCAO DAS PARTICULAS
void AdvectaParticula(PONTO *Ponto, double U, double V, double dt, double *NovoX, double *NovoY)
{
    *NovoX = Ponto->x + dt*U;
    *NovoY = Ponto->y + dt*V;
    return;
}

void VetorNormal(CURVA *Curva, PONTO *Ponto, double *X, double *Y)
{
    static PONTO *pontoAnt, *pontoProx, pontoTeste;
    static double vetX, vetY, normalX, normalY, norma;
    static double sentido;

    pontoAnt = Ponto->ant;
    pontoProx = Ponto->prox;

    vetX = (pontoProx->x - pontoAnt->x);
    vetY = (pontoProx->y - pontoAnt->y);
    normalX = -vetY;
    normalY = vetX;

    pontoTeste.x = 0.5*(pontoAnt->x + pontoProx->x) + normalX;
    pontoTeste.y = 0.5*(pontoAnt->y + pontoProx->y) + normalY;


    sentido = (pontoProx->x - pontoAnt->x)*(pontoTeste.y - pontoAnt->y) - (pontoTeste.x - pontoAnt->x)*(pontoProx->y - pontoAnt->y);

    //Sentido<0, a normal ta apontando pro lado direito da curva
    //Sentido>0, a normal ta apontando pro lado esquerdo

//    if( Ponto->x<0.012 ) {
//        PrintDebug("%lf %lf %lf %lf %d %d %lf\n", Ponto->x, Ponto->y, normalX, normalY, Curva->regiaoEsq->localRegiao, Curva->regiaoDir->localRegiao, sentido);
//        getchar();
//    }

    if( sentido<0 && Curva->regiaoEsq->localRegiao==EXTERNA) {
        normalX = -normalX;
        normalY = -normalY;
	}
	else if( sentido>0 && Curva->regiaoDir->localRegiao==EXTERNA ) {
        normalX = -normalX;
        normalY = -normalY;
//        PrintDebug("AFAFAIFAIJIAJF\n");
	}

	//Normalizando o vetor
    norma = sqrt( normalX*normalX + normalY*normalY );
    *X = normalX/norma;
    *Y = normalY/norma;

//    if( Ponto->x<0.012 ) {
//        PrintDebug("%lf %lf %lf %lf %d\n", Ponto->x, Ponto->y, *X, *Y, Curva->regiaoEsq->localRegiao);
//        getchar();
//    }

    return;
}

void VetorNormalMinimosQuadrados(INTERFACE *Inter, double Xc, double Yc, double Raio, double *ResultX, double *ResultY)
{
    double A[2][2] = {{0.0, 0.0}, {0.0, 0.0}}; //Matriz do sistema (minimos quadrados)
    double rhs[2] = {0.0, 0.0}; //Lado direito do sistema linear (minimos quadrados)
    double a; //Coeficientes da reta y=ax+b
    double det, distancia;
    PONTO *p;
    CURVA *c;

    for( c=Inter->curvaInicial; c!=NULL; c=c->proxCurva ) {
        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            distancia = sqrt( ((Xc - p->x)*(Xc - p->x)) + ((Yc - p->y)*(Yc - p->y)) );

            if( distancia > Raio )
                continue;

            ///Tirando os pontos fixos no inflow
            if( p->fixo )
                continue;

            A[0][0] += (p->x)*(p->x);
            A[0][1] += p->x;
            A[1][0] += p->x;
            A[1][1] += 1.0;
            rhs[0] += (p->x)*(p->y);
            rhs[1] += p->y;
        }
    }

    det = A[0][0]*A[1][1] - A[0][1]*A[1][0]; //Determinante
	if( fabs(det)<1e-9 ) {
		*ResultX = 1.0;
		*ResultY = 0.0;
		//DeuProblema("PROB: VETOR NORMAL MINIMOS QUADRADOS %lf %lf\n", Xc, Yc);
		return;
	}

	a = (rhs[0]*A[1][1] - rhs[1]*A[0][1])/det;
	//b = (rhs[1]*A[0][0] - rhs[0]*A[1][0])/det;


    double norma = sqrt(a*a + 1.0);
    *ResultX = -a/norma;
    *ResultY = 1.0/norma;

    //Freeflow
//    *ResultX = a/norma;
//    *ResultY = -1.0/norma;

    return;
}

double CalculaCurvatura(INTERFACE *Interf, double Xc, double Yc, double RaioCurvatura, double RaioNormal, double NxCelulas, double NyCelulas)
{
    double curvatura = 0.0;
    double **A; //Matriz do sistema linear (minimos quadrados)
    //double rhs[3] = {0.0, 0.0, 0.0}; //Lado direito do sistema linear (minimos quadrados)
    //long double x[3] = {0.0, 0.0, 0.0};
    double d/*, e, f*/; //Coeficientes da parabola dxx + ex + f
    double xi, eta;
    double nx, ny;
    double det, distancia, valor = 0.0;
    PONTO *p;
    CURVA *c;

    A = (double **)AlocaMatriz(3, 4, sizeof(double), &valor);


    //Calculando o vetor normal
    VetorNormalMinimosQuadrados(Interf, Xc, Yc, RaioNormal, &nx, &ny);

    for( c=Interf->curvaInicial; c!=NULL; c=c->proxCurva ) {

        LOOP_CURVA(c) {
            p = loop_stuff.ponto;
            distancia = sqrt( ((Xc - p->x)*(Xc - p->x)) + ((Yc - p->y)*(Yc - p->y)) );

            if( distancia > RaioCurvatura )
                continue;

            ///Tirando os pontos fixos no inflow
            if( p->fixo )
                continue;

    //        xi = ((Xc - p->x)*ny) + ((p->y - Yc)*nx);
    //        eta = ((p->x - Xc)*nx) + ((p->y - Yc)*ny);

            //Freeflow
            xi = (p->x-Xc)*ny-(p->y-Yc)*nx;
            eta = (p->x-Xc)*nx+(p->y-Yc)*ny;

            //Elementos da matriz
            A[0][0] += xi*xi*xi*xi;
            A[0][1] += xi*xi*xi;
            A[0][2] += xi*xi;
            A[1][0] += xi*xi*xi;
            A[1][1] += xi*xi;
            A[1][2] += xi;
            A[2][0] += xi*xi;
            A[2][1] += xi;
            A[2][2] += 1.0;

            //Vetor do lado direito
            A[0][3] += xi*xi*eta;
            A[1][3] += xi*eta;
            A[2][3] += eta;

    //        if( Yc==0.01 ) {
    //            cont++;
    //            PrintDebug("Particula %d: (%lf, %lf)\n", cont, p->x, p->y);
    //        }
    //
    //        fprintf(arq, "%.20lf %.20lf\n", p->x, p->y);
        }
    }

    //fclose(arq);

    /// == Resolvendo por cramer
    det = A[0][0]*A[1][1]*A[2][2] + A[2][0]*A[0][1]*A[1][2] + A[1][0]*A[2][1]*A[0][2]
            - A[2][0]*A[1][1]*A[0][2] - A[0][0]*A[2][1]*A[1][2] - A[2][2]*A[1][0]*A[0][1]; //Determinante

    //PrintDebug("%e ", det);
	if( fabs(det)<1e-9 ) {
        //PrintDebug("SISTEMA SINGULAR\n");
		return 0.0;
	}

    //Usando regra de cramer
	d = A[0][3]*A[1][1]*A[2][2] + A[2][3]*A[0][1]*A[1][2] + A[1][3]*A[2][1]*A[0][2]
            - A[2][3]*A[1][1]*A[0][2] - A[0][3]*A[2][1]*A[1][2] - A[2][2]*A[1][3]*A[0][1]; //Determinante
    d = d/det;

    curvatura = -2.0*d;

    //if( Yc==0.01 ) {
//        PrintDebug("(%lf %lf): Normal %lf %lf %lf\n", Xc, Yc, nx, ny, curvatura);
//        PrintDebug("MATRIZ:\n");
//        PrintDebug("A=[%.15e %.15e %.15e;\n", A[0][0], A[0][1], A[0][2]);
//        PrintDebug("%.15e %.15e %.15e;\n", A[1][0], A[1][1], A[1][2]);
//        PrintDebug("%.15e %.15e %.15e]\n", A[2][0], A[2][1], A[2][2]);
//        PrintDebug("b=[%.15e %.15e %.15e]\n", A[0][3], A[1][3], A[2][3]);
//        PrintDebug("\n\n");
   // }

    /// Resolvendo pelo freeflow
//    if( SOLVESYSTEM(2, A, x) ) {
//        curvatura = 0.0;
//    }
//    else {
//        curvatura = -2.0*x[0];
//    }


    if( (nx*NxCelulas+ny*NyCelulas) < 0.0 )
        curvatura = -curvatura;

    DesalocaMatriz((void **)A, 3, 4);

    //return -curvatura;

    return curvatura;
}

double AnguloVetor(double vx, double vy)
{
	double modulo = sqrt(vx*vx + vy*vy);

	if( vx>=0 && vy>=0 ) // Primeiro quadrante
		return asin(vy/modulo);
	else if( vx<0 && vy>=0 ) // Segundo quadrante
		return M_PI - asin(vy/modulo);
	else if( vx>=0 && vy<0 ) // Terceiro quadrante
		return asin(vy/modulo);
	else if( vx<0 && vy<0 ) // quarto quadrante
		return - M_PI - asin(vy/modulo);
	else
		DeuProblema("Caso nao-definido no AnguloVetor.\n");

    return 0.0;
}

double Bezier(double *Beta, double T, int N)
{
    double *betaAux, result;
    int i, j;

    betaAux = (double *)malloc( (N+1)*sizeof(double) );

    for( i=0; i<=N; i++ )
        betaAux[i] = Beta[i];

    for( j=1; j<=N; j++ ) {
        for( i=0;i<=N-j; i++ ) {
            betaAux[i] = betaAux[i]*(1.0 - T) + betaAux[i+1]*T;
        }
    }

    result = betaAux[0];
    free(betaAux);
    return result;
}

double Dif_Bezier(double *Beta, double T, int N)
{
    double *betaAux, result;
    int i, j;

    betaAux = (double *)malloc( (N+1)*sizeof(double) );

    for( i=0; i<=N; i++ )
        betaAux[i] = Beta[i];

    N = N-1;

    for( i=0; i<=N; i++ )
        betaAux[i] = betaAux[i+1] - betaAux[i];

    for( j=1; j<=N; j++ ) {
        for( i=0;i<=N-j; i++ ) {
            betaAux[i] = betaAux[i]*(1.0 - T) + betaAux[i+1]*T;
        }
    }

    result = (N+1)*betaAux[0];
    free(betaAux);
    return result;
}

double Dif2_Bezier(double *Beta, double T, int N)
{
    double *betaAux, result;
    int i, j;

    betaAux = (double *)malloc( (N+1)*sizeof(double) );

    for( i=0; i<=N; i++ )
        betaAux[i] = Beta[i];

    N = N-2;

    for( i=0; i<=N; i++ )
        betaAux[i] = betaAux[i+2] - 2*betaAux[i+1] + betaAux[i];

    for( j=1; j<=N; j++ ) {
        for( i=0;i<=N-j; i++ ) {
            betaAux[i] = betaAux[i]*(1.0 - T) + betaAux[i+1]*T;
        }
    }

    result = betaAux[0];
    free(betaAux);
    return result;
}

double CalculaCurvaturaPeloAnguloDeContato(INTERFACE *Interf, double Xc, double Yc, double RaioNormal)
{
    double ang_s = 50.0;
    double nx1, ny1, nx2, ny2, nx3, ny3;
    double ang, curvatura, beta, pi;


    ang = 0.0;
    curvatura = 0.0;
    beta =0.0;
    pi = 3.1415927;
    VetorNormalMinimosQuadrados(Interf, Xc, Yc, RaioNormal, &nx1, &ny1);
    beta = 2*Yc;

    ang = ang_s * (pi / 180.0); // convert
    nx2 = 0.0;
    ny2 =  cos(ang);

    // Calculando a normal da parede  //
    nx3 = 0.0;
    ny3 = 1.0;

    curvatura = ((nx1*nx3+ny1*ny3) /beta) - ((nx2*nx3+ny2*ny3) /beta);


    return curvatura;
}



double AreaTriangulo(PONTO *pI, PONTO *pJ, PONTO *pK)
{
    return 0.5*( (pJ->x - pI->x)*(pK->y - pI->y) - (pK->x - pI->x)*(pJ->y - pI->y) );
}

/// Funcao do freeflow usada na curvatura pro sistema linear. testando
int SOLVESYSTEM (int n, double **a, double *x)
{
  int   i,j,k;
  double  maxi, maxj, aux;
  double *y;

  /* Transformando a matriz em uma matriz triangular */
  /* superior                                        */

  /* looping do pivo */
  for (k = 0; k < n; k++) {

    /* determinado o maior elemento em valor absoluto */
    /* na coluna k a partir da linha k                */
    j    = k;
    maxj = fabs(a[k][k]);

    for (i = k+1; i <= n; i++) {

      maxi = fabs(a[i][k]);

      if (maxi > maxj) {

	j    = i;
	maxj = maxi;

      } /* if */

    } /* for i */

    /* verificando se o sistema e nao inversivel */
    if (maxj < 1.0e-10) {

      /* sistema nao inversivel */
      return 1;

    } /* if */

    /* escolhendo o pivo */
    if (j != k) {

      y    = a[j];
      a[j] = a[k];
      a[k] = y;

    } /* if */

    /* zerando os elementos na coluna k abaixo da linha k */
    for (i = k+1; i <= n; i++) {

      aux = a[i][k]/a[k][k];

      for (j = k+1; j <= n+1; j++) {

	a[i][j] -= aux*a[k][j];

      } /* for j */

    } /* for i */

  } /* for k */

  /* verificando se o sistema e nao inversivel */
  if (fabs(a[n][n]) < 1.0e-10) {

    /* sistema nao inversivel */
    return 1;

  } /* if */

  /* resolvendo o sistema triangular superior */
  x[n] = a[n][n+1]/a[n][n];

  for (i = n-1; i >= 0; i--) {

    for (j = i+1, aux = 0.0; j <= n; j++) aux += a[i][j]*x[j];

    x[i] = (a[i][n+1]-aux)/a[i][i];

  } /* for i */

  /* retornano que houve sucesso */
  return 0;

} /* solvesystem */

char CriaEmbaracamentoForcado(INTERFACE *Interface, double DistanciaMinima)
{
    CURVA *curva1, *curva2;
    PONTO *p1, *p2, *p1Escolhido, *p2Escolhido;
    double distancia, distanciaPraEmbaracar, norma, embaracou;
    double vetX, vetY;

    distanciaPraEmbaracar = DistanciaMinima;
    embaracou = 0;

    //No momento o algoritmo eh capaz de fazer 1 quebra por curva, depois vou tentar generalizar
    for( curva1=Interface->curvaInicial; curva1!=NULL; curva1=curva1->proxCurva ) {

        DistanciaMinima = distanciaPraEmbaracar;
        p1Escolhido = NULL;
        p2Escolhido = NULL;
        //printf("\n CURVA \n");

        LOOP_CURVA(curva1) {
            p1 = loop_stuff.ponto;
            for( curva2=Interface->curvaInicial; curva2!=NULL; curva2=curva2->proxCurva ) {
                LOOP_CURVA(curva2) {
                    p2 = loop_stuff.ponto;

                    //Ignora se os dois pontos forem consecutivos da mesma curva (ou se for o mesmo ponto tambem)
                    if( p1==p2 || p1->prox==p2 || p1->ant==p2 )
                        continue;
                    else if( PontosIguais(p1, p2) )
                        continue;
                    else if( p1->prox!=NULL && PontosIguais(p1->prox, p2) )
                        continue;
                    else if( p1->ant!=NULL && PontosIguais(p1->ant, p2) ) //Apenas caso tente comparar o ultimo segmento de uma curva ciclica com o primeiro
                        continue;

                    //Nao pode fazer embaracamento com node. Verifica se p1 e p2 sao nodes
                    if( p1->ant==NULL || p1->prox==NULL || p2->ant==NULL || p2->prox==NULL )
                        continue;

                    distancia = sqrt( (p2->x - p1->x)*(p2->x - p1->x) + (p2->y - p1->y)*(p2->y - p1->y) );
                    //printf("Distancia: %lf %lf\n", distancia, DistanciaMinima);
                    if( distancia<DistanciaMinima ) {
                        DistanciaMinima = distancia;
                        p1Escolhido = p1;
                        p2Escolhido = p2;
                    }
                }
            }
        }


        //Realiza o embaracamento usando o p1 e p2 escolhidos
        if( p1Escolhido!=NULL && p2Escolhido!=NULL ) {
            vetX = p2Escolhido->x - p1Escolhido->x;
            vetY = p2Escolhido->y - p1Escolhido->y;
            norma = sqrt( vetX*vetX + vetY*vetY );
            vetX = vetX/norma;
            vetY = vetY/norma;

            printf("NODE DESSA CURVA [%lf %lf]\n", curva1->pontoInicial->x, curva1->pontoInicial->y);
            printf("PONTOS ESCOLHIDOS [%lf %lf] [%lf %lf]\n", p1Escolhido->x, p1Escolhido->y, p2Escolhido->x, p2Escolhido->y);

            distanciaPraEmbaracar = DistanciaMinima + 0.00000001;
            p1Escolhido->x = p1Escolhido->x + distanciaPraEmbaracar*vetX;
            p1Escolhido->y = p1Escolhido->y + distanciaPraEmbaracar*vetY;

            p2Escolhido->x = p2Escolhido->x - distanciaPraEmbaracar*vetX;
            p2Escolhido->y = p2Escolhido->y - distanciaPraEmbaracar*vetY;
            embaracou = 1;


        }

    }


    return embaracou;
}


///FUNCOES PARA LEITURA E CRIACAO DA INTERFACE A PARTIR DO ARQUIVO
CURVA *AdicionaCurvaElipse(INTERFACE *Interface, double XCentro, double YCentro, double RaioX, double RaioY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, double Angulo, int QtdSegmentos)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    double deltaAngulo = 2*M_PI/(double)QtdSegmentos;
    int i;

    /// Diminuindo um residuo no raio
    RaioX -= 1e-08;
    RaioY -= 1e-08;

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    for( i=0; i<=QtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        //Setando as coordenadas da circunferencia
        ponto->x = XCentro + RaioX*cos(Angulo);
        ponto->y = YCentro + RaioY*sin(Angulo);
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        Angulo += deltaAngulo;
    }

    //Ultimo ponto adicionado tb eh um node
    pontoAnt->isNode = 1;
    pontoAnt->prox = NULL;

    free(proxPonto);
    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}

CURVA *AdicionaCurvaElipseAxissimetrico(INTERFACE *Interface, double XCentro, double YCentro, double RaioX, double RaioY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, double Angulo, int QtdSegmentos)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    int i;

    double anguloMin = 1.5*M_PI + 1e-8;
    double anguloMax = 2.5*M_PI - 1e-8;
    double deltaAngulo = (anguloMax - anguloMin)/(double)QtdSegmentos;

    /// Diminuindo um residuo no raio
    RaioX -= 1e-08;
    RaioY -= 1e-08;

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    Angulo = anguloMin;
    for( i=0; i<=QtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        //Setando as coordenadas da circunferencia
//        double anguloTemp = Angulo - anguloMin;
//        if( anguloTemp>0.5*M_PI )
//            anguloTemp = M_PI - anguloTemp;
////        double anguloTemp = Angulo;
//        RaioX = RaioY = 1.0 + 0.05*( 3.0*cos(anguloTemp)*cos(anguloTemp) - 1.0 );
        ponto->x = XCentro + RaioX*cos(Angulo);
        ponto->y = YCentro + RaioY*sin(Angulo);
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        Angulo += deltaAngulo;
    }

    // Imendando a curva pra virar um loop
    pontoAnt->prox = curva->pontoInicial;
    curva->pontoInicial->ant = pontoAnt;

    AdicionaCurvaNaInterface(Interface, curva);
    // FinalizaCurvasCiclicas(Interface);

    free(ponto);
    return curva;
}

CURVA *AdicionaCurvaImagemAxissimetrico(INTERFACE *Interface, char *filename, REGIAO *RegiaoEsq, REGIAO *RegiaoDir)
{
    double x, y;
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;

    curva = (CURVA *) calloc(1, sizeof(CURVA));

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *) calloc(1, sizeof(PONTO));
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    //Coletando os pontos do arquivo
    FILE *arq = fopen(filename, "rt");
    while (!feof(arq)) {
        proxPonto = (PONTO *) calloc(1, sizeof(PONTO));
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        //Setando as coordenadas
        fscanf(arq, "%lf %lf\n", &x, &y);
        ponto->x = x;
        ponto->y = y;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
    }

    // Fecha o arquivo temporário
    fclose(arq);  

    //Ultimo ponto adicionado tb eh um node
    pontoAnt->prox = curva->pontoInicial;
    curva->pontoInicial->ant = pontoAnt;

    AdicionaCurvaNaInterface(Interface, curva);
    //FinalizaCurvasCiclicas(Interface);

    free(proxPonto);
    return curva;
}

double FunctionMaziInitialShape(double r, double h_inf, double R0)
{
    double max_term = 1.0 - (r/R0)*(r/R0);
    double h = ((max_term) > 0.0) ? max_term : 0.0;
    return h_inf + R0*h;
}

CURVA *AdicionaCurvaMaziSpreading(INTERFACE *Interface, double h_inf, double R0, double max_r, REGIAO *RegiaoEsq, REGIAO *RegiaoDir)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    int i, qtdSegmentos = 400;
    double r, delta_r;

    max_r -= 1e-08;

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;
    curva->pontoInicial = ponto;



    delta_r = max_r/qtdSegmentos;
    r = 0.0;
    for( i=0; i<=qtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        ponto->x = r;
        ponto->y = FunctionMaziInitialShape(r, h_inf, R0);
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        r += delta_r;
    }
    pontoAnt->fixo = 1;

    /// Colocando o ponto pra ir pra parte de baixo-direita da camada de fluido h_inf
    proxPonto = (PONTO *)malloc( sizeof(PONTO) );
    proxPonto->isNode = 0;
    proxPonto->fixo = 0;
    proxPonto->outflow = 0;
    ponto->x = max_r;
    ponto->y = 0.0 + 1e-8;
    ponto->prox = proxPonto;
    ponto->ant = pontoAnt;
    ponto->fixo = 1;
    ponto->outflow = 0;
    pontoAnt =  ponto;
    ponto = proxPonto;

    /// Colocando o ponto pra ir pra parte de baixo-esquerda da camada de fluido h_inf
    proxPonto = (PONTO *)malloc( sizeof(PONTO) );
    proxPonto->isNode = 0;
    proxPonto->fixo = 0;
    proxPonto->outflow = 0;
    ponto->x = 0.0;
    ponto->y = 0.0 + 1e-8;
    ponto->prox = proxPonto;
    ponto->ant = pontoAnt;
    ponto->fixo = 1;
    ponto->outflow = 0;
    pontoAnt =  ponto;
    ponto = proxPonto;


    
    pontoAnt->prox = curva->pontoInicial;
    curva->pontoInicial->ant = pontoAnt;

    AdicionaCurvaNaInterface(Interface, curva);

    free(ponto);
    return curva;
}

CURVA *AdicionaCurvaSemiElipseAxissimetrico(INTERFACE *Interface, double XCentro, double YCentro, double RaioX, double RaioY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, double Angulo, int QtdSegmentos)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    int i;

    double anguloMin = 2.0*M_PI + 1e-8;
    double anguloMax = 2.5*M_PI - 1e-8;
    double deltaAngulo = (anguloMax - anguloMin)/(double)QtdSegmentos;

    /// Diminuindo um residuo no raio
    RaioX -= 1e-08;
    RaioY -= 1e-08;

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    Angulo = anguloMin;
    for( i=0; i<=QtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        //Setando as coordenadas da circunferencia
//        double anguloTemp = Angulo - anguloMin;
//        if( anguloTemp>0.5*M_PI )
//            anguloTemp = M_PI - anguloTemp;
////        double anguloTemp = Angulo;
//        RaioX = RaioY = 1.0 + 0.05*( 3.0*cos(anguloTemp)*cos(anguloTemp) - 1.0 );
        ponto->x = XCentro + RaioX*cos(Angulo);
        ponto->y = YCentro + RaioY*sin(Angulo);
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        Angulo += deltaAngulo;
    }

    // Adicionando um ponto na origem
    proxPonto = (PONTO *)malloc( sizeof(PONTO) );
    proxPonto->isNode = 0;
    proxPonto->fixo = 0;
    proxPonto->outflow = 0;
    ponto->x = XCentro;
    ponto->y = YCentro;
    ponto->prox = proxPonto;
    ponto->ant = pontoAnt;
    pontoAnt =  ponto;
    ponto = proxPonto;

    // Finalizando a curva com um ponto igual ao ponto inicial
    ponto->x = XCentro + RaioX*cos(anguloMin);
    ponto->y = YCentro + RaioY*sin(anguloMin);
    ponto->prox = NULL;
    ponto->ant = pontoAnt;
    ponto->isNode = 1;

    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}

CURVA *AdicionaCurvaSemiElipse(INTERFACE *Interface, double XCentro, double YCentro, double RaioX, double RaioY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, int QtdSegmentos)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    double deltaAngulo = M_PI/(double)QtdSegmentos;
    double angulo;
    int i;

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    angulo = 0.0;
    for( i=0; i<=QtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        //Setando as coordenadas da circunferencia
        ponto->x = XCentro + RaioX*cos(angulo);
        ponto->y = YCentro + RaioY*sin(angulo);
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        angulo += deltaAngulo;
    }

    // Ajustando os ponteiros do primeiro e ultimo ponto
    pontoAnt->prox = curva->pontoInicial;
    curva->pontoInicial->ant = pontoAnt;

    AdicionaCurvaNaInterface(Interface, curva);
    // FinalizaCurvasCiclicas(Interface);

    free(proxPonto);
    return curva;
}

CURVA *AdicionaCurvaSemiElipseCoalescencia(INTERFACE *Interface, double R0, double ContactAngle, double TranslateX, double TranslateY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, int QtdSegmentos)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    int i;

    double gap_angle = 0.0021;

    ContactAngle *= M_PI/180.0;


    // Calculating the radius and center of the full circle
    double tan_squared = tan(ContactAngle)*tan(ContactAngle);
    double R = sqrt( R0*R0*(1.0 + 1.0/tan_squared) );
    double center_x = 0.0;
    double center_y = - sqrt( R*R - R0*R0 );

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    // Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    double anguloInicial = atan(-center_y/R0) + 1e-9;
    double anguloFinal = (M_PI - anguloInicial - 1e-9) - gap_angle;
    TranslateX = R0 + 0.0001;
    
    double angulo = anguloInicial;
    double deltaAngulo = (anguloFinal - anguloInicial)/(double)QtdSegmentos;
    for( i=0; i<=QtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->outflow = 0;
        proxPonto->fixo = 0;

        //Setando as coordenadas da circunferencia
        ponto->x = center_x + R*cos(angulo) + TranslateX;
        ponto->y = center_y + R*sin(angulo) + TranslateY;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        angulo += deltaAngulo;
    }

    // Indo para a segunda gota
    TranslateX = - R0 - 0.0001;
    angulo = anguloInicial + gap_angle;
    deltaAngulo = (anguloFinal - anguloInicial)/(double)QtdSegmentos;
    for( i=0; i<=QtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->outflow = 0;
        proxPonto->fixo = 0;

        //Setando as coordenadas da circunferencia
        ponto->x = center_x + R*cos(angulo) + TranslateX;
        ponto->y = center_y + R*sin(angulo) + TranslateY;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        angulo += deltaAngulo;
    }

    //Adicionando o ultimo ponto em cima do primeiro
    ponto->x = center_x + R*cos(anguloInicial) + R0 + 0.0001;
    ponto->y = center_y + R*sin(anguloInicial) + TranslateY;
    ponto->prox = NULL;
    ponto->ant = pontoAnt;
    ponto->isNode = 1;

    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}

CURVA *AdicionaCurvaSemiElipseCoalescencia2(INTERFACE *Interface, double R0, double ContactAngle, double TranslateX, double TranslateY, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, int QtdSegmentos)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    int i;

    ContactAngle *= M_PI/180.0;


    // Calculating the radius and center of the full circle
    double tan_squared = tan(ContactAngle)*tan(ContactAngle);
    double R = sqrt( R0*R0*(1.0 + 1.0/tan_squared) );
    double center_x = 0.0;
    double center_y = - sqrt( R*R - R0*R0 );

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    // Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    double anguloInicial = atan(-center_y/R0) + 1e-9;
    double anguloFinal = (M_PI - anguloInicial - 1e-9);
    
    double angulo = anguloInicial;
    double deltaAngulo = (anguloFinal - anguloInicial)/(double)QtdSegmentos;
    for( i=0; i<=QtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->outflow = 0;
        proxPonto->fixo = 0;

        //Setando as coordenadas da circunferencia
        ponto->x = center_x + R*cos(angulo) + TranslateX;
        ponto->y = center_y + R*sin(angulo) + TranslateY;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        angulo += deltaAngulo;
    }

    //Adicionando o ultimo ponto em cima do primeiro
    ponto->x = center_x + R*cos(anguloInicial) + TranslateX;
    ponto->y = center_y + R*sin(anguloInicial) + TranslateY;
    ponto->prox = NULL;
    ponto->ant = pontoAnt;
    ponto->isNode = 1;

    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}

CURVA *AdicionaCurvaRetangulo(INTERFACE *Interface, double xMin, double xMax, double yMin, double yMax, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, double MaxDistancia)
{
    CURVA *curva;
    PONTO *p1, *p2, *p3, *p4, *p5;

    xMin += 1e-8;
    xMax -= 1e-8;
    yMin += 1e-8;
    yMax -= 1e-8;

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    p1 = (PONTO *)malloc( sizeof(PONTO) );
    p2 = (PONTO *)malloc( sizeof(PONTO) );
    p3 = (PONTO *)malloc( sizeof(PONTO) );
    p4 = (PONTO *)malloc( sizeof(PONTO) );
    p5 = (PONTO *)malloc( sizeof(PONTO) );

    curva->pontoInicial = p1;

    p1->x = xMin;
    p1->y = yMin;
    p1->isNode = 1;
    p1->fixo = 0;
    p1->ant = NULL;
    p1->prox = p2;
    p1->outflow = 0;

    p2->x = xMax;
    p2->y = yMin;
    p2->isNode = 0;
    p2->fixo = 0;
    p2->ant = p1;
    p2->prox = p3;
    p2->outflow = 0;

    p3->x = xMax;
    p3->y = yMax;
    p3->isNode = 0;
    p3->fixo = 0;
    p3->ant = p2;
    p3->prox = p4;
    p3->outflow = 0;

    p4->x = xMin;
    p4->y = yMax;
    p4->isNode = 0;
    p4->fixo = 0;
    p4->ant = p3;
    p4->prox = p5;
    p4->outflow = 0;

    p5->x = xMin;
    p5->y = yMin;
    p5->isNode = 1;
    p5->fixo = 0;
    p5->ant = p4;
    p5->prox = NULL;
    p5->outflow = 0;


    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}

CURVA *AdicionaCurvaFaradayWaves(INTERFACE *Interface, double Alfa, double Epsilon, REGIAO *RegiaoEsq, REGIAO *RegiaoDir)
{
    CURVA *curva;
    PONTO *p1, *p2, *ponto;
    int i;

    double xMin = 0.0 + 1e-8;
    double xMax = 17.62 - 1e-8;
    double yMin = 0.0 + 1e-8;

    // Alocando memoria para a curva front-tracking
    curva = (CURVA *)malloc( sizeof(CURVA) );

    // Inicializando as regioes (usado apenas em simulacoes que tenham mudanca topologica...)
    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    // Alocando memoria para dois pontos
    p1 = (PONTO *)malloc( sizeof(PONTO) );
    p2 = (PONTO *)malloc( sizeof(PONTO) );

    curva->pontoInicial = p1;

    // Criando o ponto p1: cantinho inferior-esquerdo da piscina
    p1->x = xMin;
    p1->y = yMin;
    p1->isNode = 1;
    p1->fixo = 0;
    p1->ant = NULL;
    p1->prox = p2;
    p1->outflow = 0;

    // Criando o ponto p2: cantinho inferior-direito da piscina
    p2->x = xMax;
    p2->y = yMin;
    p2->isNode = 0;
    p2->fixo = 0;
    p2->ant = p1;
    p2->prox = NULL;
    p2->outflow = 0;

    // Agora vamos criar os pontos da onda
    double qtdSegmentos = 10000;
    double delta_x = (xMax - xMin)/qtdSegmentos;
    double x_atual;
    PONTO *pontoAnterior = p2;
    PONTO *proximoPonto = NULL;

    //Criando o node inicial da onda
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;
    p2->prox = ponto;

    double h0 = 0.0013;
    double h1 = 0.00002;
    double L = 0.072;
    double L_scale = 0.00408407045;

    x_atual = xMax;
    for( i=0; i<=qtdSegmentos; i++ ) {
        proximoPonto = (PONTO *)malloc( sizeof(PONTO) );
        proximoPonto->isNode = 0;
        proximoPonto->fixo = 0;
        proximoPonto->outflow = 0;

        // Colocando as coordenadas e os ponteiros da lista encadeada de pontos
        ponto->x = x_atual;
//        ponto->y = (Alfa/M_PI)*(1.0 + Epsilon*sin(M_PI*x_atual - 0.5*M_PI));
        double x_dim = x_atual*L_scale;
        ponto->y = h0 - h1 + 0.0001*exp( -(10.0/L)*(10.0/L)*(x_dim - 0.5*L)*(x_dim - 0.5*L) );
        ponto->y /= L_scale; //Adimensionalizando
        ponto->prox = proximoPonto;
        ponto->ant = pontoAnterior;

        //Preparando para a proxima iteracao
        pontoAnterior =  ponto;
        ponto = proximoPonto;
        x_atual -= delta_x;
    }

    // Fechando o poligono: repetindo o ponto p1 (cantinho inferior-esquerdo)
    ponto->x = xMin;
    ponto->y = yMin;
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;
    ponto->ant = pontoAnterior;
    ponto->prox = NULL;

    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}

CURVA *AdicionaCurvaInflow(INTERFACE *Interface, double Posicao, double ComprimentoInicial, double CoordMin, double CoordMax, double TamanhoBorda, char Direcao, REGIAO *RegiaoEsq, REGIAO *RegiaoDir, int QtdSegmentos)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    int i, metadeSegmentos;
    double delta;

    curva = (CURVA *)malloc( sizeof(CURVA) );
    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    metadeSegmentos = QtdSegmentos/2;
    CoordMin += TamanhoBorda;
    CoordMax -= TamanhoBorda;
    delta = (CoordMax - CoordMin)/metadeSegmentos;

    //Colocando a metade dos segmentos que vao ficar livres
    for( i=0; i<=metadeSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;

        if( Direcao=='v' ) {
            ponto->x = Posicao + ComprimentoInicial;
            ponto->y = CoordMin + delta*i;
        }
        else if( Direcao=='h' ) {
            ponto->x = CoordMin + delta*i;
            ponto->y = Posicao - ComprimentoInicial - 1e-8;
        }
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        ponto->fixo = 0;
        ponto->outflow = 0;

        /// Mantendo os cantinhos fixos
        if( i==0 || (i==metadeSegmentos) )
            ponto->fixo = 1;

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
    }

    ///TIRAR ISSO
    //pontoAnt->fixo = 1;

    //Colocando a metade dos segmentos que vao ficar presos
    for( i=metadeSegmentos; i>=0; i-- ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;

        if( Direcao=='v' ) {
            ponto->x = Posicao;
            ponto->y = CoordMin + delta*i;
        }
        else if( Direcao=='h' ) {
            ponto->x = CoordMin + delta*i;
            ponto->y = Posicao - 1e-8;
        }
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;
        ponto->fixo = 1;
        ponto->outflow = 0;

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
    }

    //Adicionando mais um pra fechar a curva ciclica, caso ComprimentoInicial>0
    if( fabs(ComprimentoInicial)>1e-8 ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        if( Direcao=='v' ) {
            ponto->x = Posicao + ComprimentoInicial;
            ponto->y = CoordMin;
        }
        else if( Direcao=='h' ) {
            ponto->x = CoordMin;
            ponto->y = Posicao - ComprimentoInicial - 1e-8;
        }
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        ///VOLTAR ISSO PRA 0
        ponto->fixo = 1;
        pontoAnt =  ponto;
        ponto = proxPonto;
    }

    //Ultimo ponto adicionado tb eh um node
    pontoAnt->isNode = 1;
    pontoAnt->prox = NULL;
    free(proxPonto);


    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);

    ///Tirar isso
    //Interface->curvaInicial->nodeInicial->fixo = 1;

    return curva;
}

CURVA *AdicionaCurvaFilamento(INTERFACE *Interface, double Rp, double R0, double Lz, REGIAO *RegiaoEsq, REGIAO *RegiaoDir)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    int i, QtdSegmentosVertical = 100;
    double z, deltaZ;

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 1;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    z = -0.5*Lz;
    deltaZ = Lz/QtdSegmentosVertical;
    for( i=0; i<=QtdSegmentosVertical; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        ponto->x = FuncaoFilamento(z, R0, Lz);
        ponto->y = z;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        AjustaExtremosFilamento(&(ponto->x), &(ponto->y), Lz);

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        z += deltaZ;
    }
    pontoAnt->fixo = 1;

    /// Colocando o ponto pra ir pro outro lado
    proxPonto = (PONTO *)malloc( sizeof(PONTO) );
    proxPonto->isNode = 0;
    proxPonto->fixo = 0;
    proxPonto->outflow = 0;

    ponto->x = - FuncaoFilamento(0.5*Lz, R0, Lz);
    ponto->y = 0.5*Lz;
    ponto->prox = proxPonto;
    ponto->ant = pontoAnt;
    ponto->fixo = 1;
    ponto->outflow = 0;

    AjustaExtremosFilamento(&(ponto->x), &(ponto->y), Lz);

    //Preparando pra proxima iteracao
    pontoAnt =  ponto;
    ponto = proxPonto;

    z = 0.5*Lz;
    z -= deltaZ;
    for( i=1; i<=QtdSegmentosVertical; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        ponto->x = - FuncaoFilamento(z, R0, Lz);
        ponto->y = z;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        AjustaExtremosFilamento(&(ponto->x), &(ponto->y), Lz);


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        z -= deltaZ;
    }
    pontoAnt->fixo = 1;

    /// Colocando o ponto pra ir pro outro lado
    proxPonto = (PONTO *)malloc( sizeof(PONTO) );
    proxPonto->isNode = 0;
    proxPonto->fixo = 0;
    proxPonto->fixo = 0;

    ponto->x = FuncaoFilamento(-0.5*Lz, R0, Lz);
    ponto->y = -0.5*Lz;
    ponto->prox = proxPonto;
    ponto->ant = pontoAnt;
    ponto->fixo = 1;

    AjustaExtremosFilamento(&(ponto->x), &(ponto->y), Lz);

    //Preparando pra proxima iteracao
    pontoAnt =  ponto;
    ponto = proxPonto;


    //Ultimo ponto adicionado tb eh um node
    pontoAnt->isNode = 1;
    pontoAnt->prox = NULL;

    free(proxPonto);
    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}

CURVA *AdicionaCurvaCaberAxi(INTERFACE *Interface, double Rp, double R0, double Lz, REGIAO *RegiaoEsq, REGIAO *RegiaoDir)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    int i, QtdSegmentosVertical = 400;
    double z, deltaZ;

    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 1;
    ponto->outflow = 0;
    ponto->x = 0.0;
    ponto->y = -0.5*Lz;
    ponto->ant = NULL;

    curva->pontoInicial = ponto;

    pontoAnt = ponto;
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->fixo = 1;
    ponto->isNode = 0;
    ponto->outflow = 0;

    pontoAnt->prox = ponto;


    z = -0.5*Lz;
    deltaZ = Lz/QtdSegmentosVertical;
    for( i=0; i<=QtdSegmentosVertical; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        ponto->x = FuncaoFilamento(z, R0, Lz);
        ponto->y = z;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        AjustaExtremosFilamento(&(ponto->x), &(ponto->y), Lz);

        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        z += deltaZ;
    }
    pontoAnt->fixo = 1;

    /// Colocando o ponto pra ir pro outro lado
    proxPonto = (PONTO *)malloc( sizeof(PONTO) );
    proxPonto->isNode = 0;
    proxPonto->fixo = 0;
    proxPonto->outflow = 0;

    ponto->x = 0.0;
    ponto->y = 0.5*Lz;
    ponto->prox = proxPonto;
    ponto->ant = pontoAnt;
    ponto->fixo = 1;
    ponto->outflow = 0;
    AjustaExtremosFilamento(&(ponto->x), &(ponto->y), Lz);

    //Preparando pra proxima iteracao
    pontoAnt =  ponto;
    ponto = proxPonto;

    /// Finalizando repetindo o ponto inicial
    ponto->isNode = 1;
    ponto->fixo = 1;
    ponto->outflow = 0;
    ponto->x = 0.0;
    ponto->y = -0.5*Lz;
    ponto->ant = pontoAnt;
    ponto->prox = NULL;
    AjustaExtremosFilamento(&(ponto->x), &(ponto->y), Lz);

    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}



CURVA *AdicionaCurvaBezierUniforme(INTERFACE *Interface, CURVA *CurvaOriginal)
{
    CURVA *curva;
    PONTO *p;
    int qtdPontos;

    curva = CurvaOriginal;

    //Contando quantos pontos vou usar na curva original
    qtdPontos = 0;
    int delta_pontos = 1;
    int contador = 0;
    LOOP_CURVA(curva) {
        p = loop_stuff.ponto;

        if( contador==0 ) {
            qtdPontos++;
            contador = delta_pontos-1;
        }
        else
            contador--;
    }

    qtdPontos++; // Incrementando um pra incluir o ultimo ponto

    // Construindo os vetores x e y pra usar na funcao de Bezier
    double *x, *y;
    x = (double *)malloc( qtdPontos*sizeof(double) );
    y = (double *)malloc( qtdPontos*sizeof(double) );

    int i = 0;
    contador = 0;
    LOOP_CURVA(curva) {
        p = loop_stuff.ponto;

//        PrintDebug("HEHUEHUE %lf %lf\n", p->x, p->y);
        if( contador==0 ) {
            x[i] = p->x;
            y[i] = p->y;
            contador = delta_pontos-1;
            i++;
        }
        else
            contador--;
    }
//    PrintDebug("HEHUEHUE %lf %lf\n", p->x, p->y);

    // Inserindo o ultimo
    x[i] = p->x;
    y[i] = p->y;

//    PrintDebug("AAAAAA %lf %lf %lf %lf\n", x[0], y[0], x[i], y[i]);
//    getchar();

//    double x1 = Bezier(x, 0.5, qtdPontos-1);
//    double y1 = Bezier(y, 0.5, qtdPontos-1);
//    DeuProblema("%lf %lf\n", x1, y1);


    //Usando integracao numerica pra calcular o comprimento da curva de Bezier
    int qtdIntervalos = 500;
    double dt = 1.0/qtdIntervalos;
    double t = 0;
    double comprimento = 0.0, fa, fb, deriv_x, deriv_y, a, b;
    for( i=0; i<qtdIntervalos; i++ ) {
        a = t;
        b = t + dt;


        deriv_x = Dif_Bezier(x, a, qtdPontos-1);
        deriv_y = Dif_Bezier(y, a, qtdPontos-1);
//        deriv_x = -sin(a);
//        deriv_y = cos(a);
        fa = sqrt( deriv_x*deriv_x + deriv_y*deriv_y );

        deriv_x = Dif_Bezier(x, b, qtdPontos-1);
        deriv_y = Dif_Bezier(y, b, qtdPontos-1);
//        deriv_x = -sin(b);
//        deriv_y = cos(b);
        fb = sqrt( deriv_x*deriv_x + deriv_y*deriv_y );

//        PrintDebug("AAA %lf %lf %g\n", deriv_x, deriv_y, 0.5*(b - a)*(fa + fb));
        comprimento += 0.5*(b - a)*(fa + fb);

        t = t + dt;
    }


    // Quantidade de segmentos para a nova curva
    double espacamento_aprox = 0.005, h;
    int qtdSegmentos = (int)( comprimento/espacamento_aprox );
	h = comprimento/qtdSegmentos;

	// Criando a nova curva
	curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = NULL;
    curva->regiaoDir = NULL;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    PONTO *ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;

    curva->pontoInicial = ponto;

    double t1 = 0.0;
    PONTO *proxPonto = NULL, *pontoAnt = NULL;
    while( t1<1.0 && fabs(t1-1.0)>1e-6 ) {
//    for( i=0; i<qtdSegmentos-1; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;

        //Calculando o novo ponto
        double x1 = Bezier(x, t1, qtdPontos-1);
		double y1 = Bezier(y, t1, qtdPontos-1);

        //Setando as coordenadas da circunferencia
        ponto->x = x1;
        ponto->y = y1;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        double chute = (t1 + 0.05) > 1.0 ? 1.0 : t1 + 0.05;

        double t2 = AlgoritmoBisseccaoBezier(x, y, qtdPontos-1, t1, chute, x1, y1, h);


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        t1 = t2;
//
//        PrintDebug("T1 %.20e %lf %lf\n", t1, x1, y1);
    }

    ponto->x = Bezier(x, 1.0, qtdPontos-1);
    ponto->y = Bezier(y, 1.0, qtdPontos-1);
    ponto->prox = NULL;
    ponto->ant = pontoAnt;
    ponto->isNode = 1;

    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);

    return curva;
}

void AdicionaCurvaSessilePendingDrops(INTERFACE *Interface, double RaioCima, double RaioBaixo, double MeioX, double MeioY, double MinY, double MaxY, int BaseFixa)
{
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    CURVA *curva;
    int i;
    int qtdSegmentos = 100;

    /// === Fazendo a circunferencia de baixo
    double centroX, centroY, theta;
    double thetaMin, thetaMax, deltaTheta;

    centroX = MeioX;
    centroY = MeioY - RaioBaixo;
    thetaMin = asin( (MinY - centroY)/RaioBaixo );
    thetaMax = atan2(sin(thetaMin), - cos(thetaMin));
    thetaMax += (thetaMax < thetaMin) ? 2.0*M_PI : 0.0;
    thetaMin += 1e-7;
    thetaMax -= 1e-7;
    deltaTheta = (thetaMax - thetaMin)/(double)qtdSegmentos;

    curva = (CURVA *)malloc( sizeof(CURVA) );
    curva->regiaoEsq = NULL;
    curva->regiaoDir = NULL;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = BaseFixa;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    theta = thetaMin;
    for( i=0; i<=qtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        //Setando as coordenadas da circunferencia
        ponto->x = centroX + RaioBaixo*cos(theta);
        ponto->y = centroY + RaioBaixo*sin(theta);
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        theta += deltaTheta;
    }

    // Fechando a curva
    ponto->x = centroX + RaioBaixo*cos(thetaMin);
    ponto->y = centroY + RaioBaixo*sin(thetaMin);
    ponto->prox = NULL;
    ponto->ant = pontoAnt;
    ponto->isNode = 1;

    ponto->fixo = BaseFixa;
    pontoAnt->fixo = BaseFixa;

    AdicionaCurvaNaInterface(Interface, curva);



    /// Fazendo a curva de cima
    centroX = MeioX;
    centroY = MeioY + RaioCima;
    thetaMax = asin( (MaxY - centroY)/RaioCima );
    thetaMin = atan2(sin(thetaMax), - cos(thetaMax));
    thetaMin -= (thetaMax < thetaMin) ? 2.0*M_PI : 0.0;
    thetaMin += 1e-7;
    thetaMax -= 1e-7;
    deltaTheta = (thetaMax - thetaMin)/(double)qtdSegmentos;

    curva = (CURVA *)malloc( sizeof(CURVA) );
    curva->regiaoEsq = NULL;
    curva->regiaoDir = NULL;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    proxPonto = pontoAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = BaseFixa;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    theta = thetaMin;
    for( i=0; i<=qtdSegmentos; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        //Setando as coordenadas da circunferencia
        ponto->x = centroX + RaioCima*cos(theta);
        ponto->y = centroY + RaioCima*sin(theta);
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        theta += deltaTheta;
    }

    // Fechando a curva
    ponto->x = centroX + RaioCima*cos(thetaMin);
    ponto->y = centroY + RaioCima*sin(thetaMin);
    ponto->prox = NULL;
    ponto->ant = pontoAnt;
    ponto->isNode = 1;

    ponto->fixo = BaseFixa;
    pontoAnt->fixo = BaseFixa;

    AdicionaCurvaNaInterface(Interface, curva);


    FinalizaCurvasCiclicas(Interface);
    return;
}

CURVA *AdicionaCurvaSplinesUniforme(INTERFACE *Interface, CURVA *CurvaOriginal)
{
    CURVA *curva;
    PONTO *p;
    int qtdPontos;

    curva = CurvaOriginal;

//    IniciaTimer();
    //Contando quantos pontos vou usar na curva original
    qtdPontos = 0;
    int delta_pontos = 1;
    int contador = 0;
    for( p=curva->pontoInicial; p->prox!=NULL; p=p->prox ) {

        if( contador==0 ) {
            qtdPontos++;
            contador = delta_pontos-1;
        }
        else
            contador--;
    }

    qtdPontos++; // Incrementando um pra incluir o ultimo ponto

    // Construindo os vetores x e y pra usar na funcao de Bezier
    double *x, *y, *t_vec;
    x = (double *)malloc( qtdPontos*sizeof(double) );
    y = (double *)malloc( qtdPontos*sizeof(double) );
    t_vec = (double *)malloc( qtdPontos*sizeof(double) );

    int i = 0;
    double delta_t = 1.0/(qtdPontos-1);
    contador = 0;
    for( p=curva->pontoInicial; p->prox!=NULL; p=p->prox ) {

//        PrintDebug("HEHUEHUE %lf %lf\n", p->x, p->y);
        if( contador==0 ) {
            x[i] = p->x;
            y[i] = p->y;
            t_vec[i] = i*delta_t;
            contador = delta_pontos-1;
            i++;
        }
        else
            contador--;

    }
//    PrintDebug("HEHUEHUE %lf %lf\n", p->x, p->y);

    // Inserindo o ultimo
    x[i] = p->x;
    y[i] = p->y;
    t_vec[i] = i*delta_t;

    //Calculando os coeficientes das splines
    POL3 *s_x, *s_y;
    s_x = (POL3 *)malloc((qtdPontos)*sizeof(POL3)); //alocando memoria para os polinomios
    s_y = (POL3 *)malloc((qtdPontos)*sizeof(POL3)); //alocando memoria para os polinomios
    SplineCubica(t_vec, x, s_x, qtdPontos-1); //Calculando os polinomios
    SplineCubica(t_vec, y, s_y, qtdPontos-1); //Calculando os polinomios
//    PrintDebug("coeficientes: %lf segundos\n", FinalizaTimer());

//    IniciaTimer();
    //Usando integracao numerica pra calcular o comprimento da curva
    int qtdIntervalos = 500;
    double dt = 1.0/qtdIntervalos;
    double t = 0;
    double comprimento = 0.0, fa, fb, deriv_x, deriv_y, a, b;
    for( i=0; i<qtdIntervalos; i++ ) {
        a = t;
        b = t + dt;

        deriv_x = ValorSplineDerivada(a, t_vec, s_x, qtdPontos-1);
        deriv_y = ValorSplineDerivada(a, t_vec, s_y, qtdPontos-1);
        fa = sqrt( deriv_x*deriv_x + deriv_y*deriv_y );

        deriv_x = ValorSplineDerivada(b, t_vec, s_x, qtdPontos-1);
        deriv_y = ValorSplineDerivada(b, t_vec, s_y, qtdPontos-1);
//        deriv_x = -sin(b);
//        deriv_y = cos(b);
        fb = sqrt( deriv_x*deriv_x + deriv_y*deriv_y );

//        PrintDebug("AAA %lf %lf %g\n", deriv_x, deriv_y, 0.5*(b - a)*(fa + fb));
        comprimento += 0.5*(b - a)*(fa + fb);

        t = t + dt;
    }
//    PrintDebug("Integracao: %lf segundos\n", FinalizaTimer());

//    IniciaTimer();
    // Quantidade de segmentos para a nova curva
    double espacamento_aprox = 0.0004, h;
    int qtdSegmentos = (int)( comprimento/espacamento_aprox );
	h = comprimento/qtdSegmentos;

	// Criando a nova curva
	curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = NULL;
    curva->regiaoDir = NULL;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    PONTO *ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;
    ponto->outflow = 0;

    curva->pontoInicial = ponto;

    double t1 = 0.0;
    PONTO *proxPonto = NULL, *pontoAnt = NULL;
    while( t1<1.0 && fabs(t1-1.0)>1e-6 ) {
//    for( i=0; i<qtdSegmentos-1; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;
        proxPonto->outflow = 0;

        //Calculando o novo ponto
		double x1 = ValorSpline(t1, t_vec, s_x, qtdPontos-1);
		double y1 = ValorSpline(t1, t_vec, s_y, qtdPontos-1);

        //Setando as coordenadas da circunferencia
        ponto->x = x1;
        ponto->y = y1;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;

        double chute = (t1 + 0.05) > 1.0 ? 1.0 : t1 + 0.05;

        double t2 = AlgoritmoBisseccaoSpline(s_x, s_y, t_vec, qtdPontos-1, t1, chute, x1, y1, h);


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        t1 = t2;
//
//        PrintDebug("T1 %.20e %lf %lf\n", t1, x1, y1);
    }

    ponto->x = ValorSpline(1.0, t_vec, s_x, qtdPontos-1);
    ponto->y = ValorSpline(1.0, t_vec, s_y, qtdPontos-1);
    ponto->prox = NULL;
    ponto->ant = pontoAnt;
    ponto->isNode = 1;
//    PrintDebug("nova curva: %lf segundos\n", FinalizaTimer());

    free(x);
    free(y);
    free(t_vec);
    free(s_x);
    free(s_y);

    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);

    return curva;
}

//Funcao que resolve um sistema tridiagonal
void Thomas(int NumEquacoes, double *Diag, double *DiagInf, double *DiagSup, double *Ind, double *Solucao)
{
    double aux, *diagAux; //vamos usar uma diagonal auxiliar, pois nao queremos alterar o valor da original
    int i;
    diagAux = (double *)malloc(NumEquacoes*sizeof(double)); //Alocando espa�o

    diagAux[0] = Diag[0]; //Fazendo a copia do primeiro elemento
    for(i=1; i<NumEquacoes; i++)
    {
        diagAux[i] = Diag[i]; //Copiando
        aux = DiagInf[i]/diagAux[i-1];
        diagAux[i] = diagAux[i] - aux*DiagSup[i-1];
        Ind[i] = Ind[i] - aux*Ind[i-1];
    }
    Solucao[NumEquacoes-1] = Ind[NumEquacoes-1]/diagAux[NumEquacoes-1];
    for(i=NumEquacoes-2; i>=0; i--) {
        Solucao[i]=(Ind[i] - DiagSup[i]*Solucao[i+1])/diagAux[i];
    }

    free(diagAux); //liberando o espa�o
    return;
}

/*
    Funcao SplineCubica:
    *X apontar� para os elementos de um vetor contendo os valores de X0, x1, x2, ..., Xn
    *Y aponta para o valor dado da funcao nos pontos Xi
    *S ser� o array que vamos formar com os n polin�mios
*/
void SplineCubica(double *X, double *Y, POL3 *S, int N)
{
    double *diagInf, *diagSup, *diag, //Diagonais inferior, superior e principal respectivamente da matriz dos coeficientes do sistema
           *ind, //Termos independentes do sistema (Segundo membro da equacao)
           *g, //Derivada segunda de sk no ponto xk com k = 0...N   (Solucao do sistema)
            h, h2; //Duas variaveis usadas para armazenas valores de espacamento h
    int k;

    diagInf = (double *)malloc(N*sizeof(double));
    diagSup = (double *)malloc(N*sizeof(double));
    diag = (double *)malloc(N*sizeof(double));
    ind = (double *)malloc(N*sizeof(double));
    g = (double *)malloc((N+1)*sizeof(double));

    for(k = 1; k<N; k++) //Montando as diagonais da matriz dos coeficientes e os termos independentes do sistema
    {
        h = X[k] - X[k-1]; //h[k]
        h2 = X[k+1] - X[k]; //h[k+1]
        diagSup[k] = h2;
        diagInf[k] = h;
        diag[k] = 2*(h+h2);
        ind[k] = 6*( (Y[k+1] - Y[k]) / h2 ) - 6*( (Y[k] - Y[k-1]) / h );
    }

    Thomas(N-1, diag, diagInf, diagSup, ind, g); //Resolvendo o sistema
    g[0] = g[N] = 0; //Condicoes impostas

    for(k = 1; k<=N; k++) //Finalmente montando os polinomios
    {
        h = X[k] - X[k-1]; //h[k]
        S[k].a =  (g[k] - g[k-1])/(6.0*h); //Calculando a[k]
        S[k].b = g[k]/2.0;
        S[k].c = ( (Y[k] - Y[k-1])/h ) + ( (2*h*g[k] + g[k-1]*h)/6.0 );
        S[k].d = Y[k];
    }

    free(diagInf);
    free(diagSup);
    free(diag);
    free(ind);
    free(g);

    return;
}

void SplineCubicaNotAKnot(double *X, double *Y, POL3 *S, int N)
{
    double *diagInf, *diagSup, *diag, //Diagonais inferior, superior e principal respectivamente da matriz dos coeficientes do sistema
           *ind, //Termos independentes do sistema (Segundo membro da equacao)
           *g, //Derivada segunda de sk no ponto xk com k = 0...N   (Solucao do sistema)
            h, h2; //Duas variaveis usadas para armazenas valores de espacamento h
    int k;

    diagInf = (double *)malloc((N-1)*sizeof(double));
    diagSup = (double *)malloc((N-1)*sizeof(double));
    diag = (double *)malloc((N-1)*sizeof(double));
    ind = (double *)malloc((N-1)*sizeof(double));
    g = (double *)malloc((N+1)*sizeof(double));

    for(k = 1; k<N; k++) //Montando as diagonais da matriz dos coeficientes e os termos independentes do sistema
    {
        h = X[k] - X[k-1]; //h[k]
        h2 = X[k+1] - X[k]; //h[k+1]

        diag[k-1] = 2.0*(h+h2);
        diagSup[k-1] = h2;
        diagInf[k-1] = h;

        // Condicoes de contorno not-a-knot
        if( k==1 ) {
            diag[k-1] += h*( 1.0 + (h/h2) );
            diagSup[k-1] += h*( -h/h2 );
        }
        else if( k==N-1 ) {
            diag[k-1] += h2*( 1.0 + (h2/h) );
            diagInf[k-1] += h2*( -h2/h );
        }

        ind[k-1] = 3.0*( (Y[k+1] - Y[k])/h2 ) - 3.0*( (Y[k] - Y[k-1])/h );
    }

    Thomas(N-1, diag, diagInf, diagSup, ind, g+1); //Resolvendo o sistema

    h = X[1] - X[0]; //h[k]
    h2 = X[2] - X[1]; //h[k+1]
    g[0] = g[1] - (h/h2)*(g[2] - g[1]);

    h = X[N-1] - X[N-2]; //h[k]
    h2 = X[N] - X[N-1]; //h[k+1]
    g[N] = g[N-1] + (h2/h)*(g[N-1] - g[N-2]);


    for(k = 0; k<N; k++) //Finalmente montando os polinomios
    {
        h = X[k+1] - X[k]; //h[k]
        S[k].a = (g[k+1] - g[k])/(3.0*h); //Calculando a[k]
        S[k].b = g[k];
        S[k].c = ( (Y[k+1] - Y[k])/h ) - (h/3.0)*(2.0*g[k] + g[k+1]);
        S[k].d = Y[k];
    }

    free(diagInf);
    free(diagSup);
    free(diag);
    free(ind);
    free(g);

    return;
}

double ValorSpline(double t, double *T, POL3 *S, int N)
{
    int i;

    for(i = 0; i<N; i++) {
        if( t<=T[i+1] ) {
            t = t - T[i];
            return S[i].a*t*t*t + S[i].b*t*t + S[i].c*t + S[i].d;
        }
    }

    printf("VALOR SPLINE: PROBLEMA\n\n");
    exit(0);
    return 0.0;
}

double ValorSplineDerivada(double t, double *T, POL3 *S, int N)
{
    int i;

    for(i = 0; i<N; i++) {
        if( t<=T[i+1] ) {
            t = t - T[i];
            return 3.0*S[i].a*t*t + 2.0*S[i].b*t + S[i].c;
        }
    }

    printf("VALOR SPLINE DERIV: PROBLEMA\n\n");
    exit(0);
    return 0.0;
}

double ValorSplineDerivada2(double t, double *T, POL3 *S, int N)
{
    int i;

    for(i = 0; i<N; i++) {
        if( t<=T[i+1] ) {
            t = t - T[i];
            return 6.0*S[i].a*t + 2.0*S[i].b;
        }
    }

    printf("VALOR SPLINE DERIV2: PROBLEMA\n\n");
    exit(0);
    return 0.0;
}

CURVA *AdicionaCurvaBezierUniforme2(INTERFACE *Interface, CURVA *CurvaOriginal)
{
    CURVA *curva;
    PONTO *p;
    int qtdPontos;

    curva = CurvaOriginal;

    //Contando quantos pontos vou usar na curva original
    qtdPontos = 0;
    int delta_pontos = 5;
    int contador = 0;
    LOOP_CURVA(curva) {
        p = loop_stuff.ponto;

        if( contador==0 ) {
            qtdPontos++;
            contador = delta_pontos-1;
        }
        else
            contador--;
    }

    qtdPontos++; // Incrementando um pra incluir o ultimo ponto

    // Construindo os vetores x e y pra usar na funcao de Bezier
    double *x, *y;
    x = (double *)malloc( qtdPontos*sizeof(double) );
    y = (double *)malloc( qtdPontos*sizeof(double) );

    int i = 0;
    contador = 0;
    LOOP_CURVA(curva) {
        p = loop_stuff.ponto;

//        PrintDebug("HEHUEHUE %lf %lf\n", p->x, p->y);
        if( contador==0 ) {
            x[i] = p->x;
            y[i] = p->y;
            contador = delta_pontos-1;
            i++;
        }
        else
            contador--;
    }
//    PrintDebug("HEHUEHUE %lf %lf\n", p->x, p->y);

    // Inserindo o ultimo
    x[i] = p->x;
    y[i] = p->y;


    // Quantidade de segmentos para a nova curva
    int qtdSegmentos = 2000;
	double h = 1.0/qtdSegmentos;

	// Criando a nova curva
	curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = NULL;
    curva->regiaoDir = NULL;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    PONTO *ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;

    curva->pontoInicial = ponto;

    double t1 = 0.0;
    PONTO *proxPonto = NULL, *pontoAnt = NULL;
    for( i=0; i<=qtdSegmentos; i++ ) {
//    for( i=0; i<qtdSegmentos-1; i++ ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;

        //Calculando o novo ponto
        double x1 = Bezier(x, t1, qtdPontos-1);
		double y1 = Bezier(y, t1, qtdPontos-1);

        //Setando as coordenadas da circunferencia
        ponto->x = x1;
        ponto->y = y1;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
        t1 = t1 + h;
    }

    pontoAnt->x = Bezier(x, 0.0, qtdPontos-1);
    pontoAnt->y = Bezier(y, 0.0, qtdPontos-1);
    pontoAnt->prox = NULL;
    pontoAnt->isNode = 1;

    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);

    return curva;
}

double AlgoritmoBisseccaoBezier(double *BetaX, double *BetaY, int N, double T1, double T2, double X_Ant, double Y_Ant, double h)
{
    if( fabs(T2-T1)<1e-10 )
        return 0.5*(T2 + T1);

    double meio = 0.5*(T2 + T1);
    double bezier_x = Bezier(BetaX, T1, N);
    double bezier_y = Bezier(BetaY, T1, N);
    double valor1 = (X_Ant - bezier_x)*(X_Ant - bezier_x) + (Y_Ant - bezier_y)*(Y_Ant - bezier_y) - h*h;

    bezier_x = Bezier(BetaX, meio, N);
    bezier_y = Bezier(BetaY, meio, N);
    double valorMeio = (X_Ant - bezier_x)*(X_Ant - bezier_x) + (Y_Ant - bezier_y)*(Y_Ant - bezier_y) - h*h;

//    PrintDebug("Bissec: %lf %lf %lf %lf %lf\n", T1, meio, valor1, valorMeio, h);

    if( valor1*valorMeio<=0 )
        return AlgoritmoBisseccaoBezier(BetaX, BetaY, N, T1, meio, X_Ant, Y_Ant, h);
    else
        return AlgoritmoBisseccaoBezier(BetaX, BetaY, N, meio, T2, X_Ant, Y_Ant, h);
}

double AlgoritmoBisseccaoSpline(POL3 *S_x, POL3 *S_y, double *T_vec, int N, double T1, double T2, double X_Ant, double Y_Ant, double h)
{
    if( fabs(T2-T1)<1e-5 )
        return 0.5*(T2 + T1);

    double meio = 0.5*(T2 + T1);
    double bezier_x = ValorSpline(T1, T_vec, S_x, N);
    double bezier_y = ValorSpline(T1, T_vec, S_y, N);

    double valor1 = (X_Ant - bezier_x)*(X_Ant - bezier_x) + (Y_Ant - bezier_y)*(Y_Ant - bezier_y) - h*h;

    bezier_x = ValorSpline(meio, T_vec, S_x, N);
    bezier_y = ValorSpline(meio, T_vec, S_y, N);
    double valorMeio = (X_Ant - bezier_x)*(X_Ant - bezier_x) + (Y_Ant - bezier_y)*(Y_Ant - bezier_y) - h*h;

//    PrintDebug("Bissec: %lf %lf %lf %lf %lf\n", T1, meio, valor1, valorMeio, h);

    if( valor1*valorMeio<=0 )
        return AlgoritmoBisseccaoSpline(S_x, S_y, T_vec, N, T1, meio, X_Ant, Y_Ant, h);
    else
        return AlgoritmoBisseccaoSpline(S_x, S_y, T_vec, N, meio, T2, X_Ant, Y_Ant, h);
}

void AjustaExtremosFilamento(double *X, double *Y, double Lz)
{
//    PrintDebug("PONTO %.15lf %.15lf\n", *X, *Y);
    if( *X>=1.0  ) {
        (*X) = 1.0 - 1e-5;
//        PrintDebug("mudou PONTO %.15lf %.15lf\n", *X, *Y);
    }
    else if( *X<=0.0 )
        *X = 0.0 + 1e-8;

    if( (*Y)>=0.5*Lz )
        (*Y) = 0.5*Lz - 1e-8;
    else if( (*Y)<=-0.5*Lz )
        (*Y) = -0.5*Lz + 1e-8;
}

CURVA *AdicionaCurvaDoArquivo(INTERFACE *Interface, const char *NomeArquivo, REGIAO *RegiaoEsq, REGIAO *RegiaoDir)
{
    CURVA *curva;
    PONTO *ponto, *proxPonto = NULL, *pontoAnt = NULL;
    double x, y;
    FILE *arq;

    arq = fopen(NomeArquivo, "rt");


    curva = (CURVA *)malloc( sizeof(CURVA) );

    curva->regiaoEsq = RegiaoEsq;
    curva->regiaoDir = RegiaoDir;
    curva->proxCurva = NULL;
    curva->curvaAnt = NULL;

    //Criando o node inicial dessa curva
    ponto = (PONTO *)malloc( sizeof(PONTO) );
    ponto->isNode = 1;
    ponto->fixo = 0;

    curva->pontoInicial = ponto;

    while( 2==fscanf(arq, "%lf %lf", &x, &y) ) {
        proxPonto = (PONTO *)malloc( sizeof(PONTO) );
        proxPonto->isNode = 0;
        proxPonto->fixo = 0;

        //Setando as coordenadas da circunferencia
        ponto->x = x;
        ponto->y = y;
        ponto->prox = proxPonto;
        ponto->ant = pontoAnt;


        //Preparando pra proxima iteracao
        pontoAnt =  ponto;
        ponto = proxPonto;
    }

    fclose(arq);

    //Ultimo ponto adicionado tb eh um node
    pontoAnt->isNode = 1;
    pontoAnt->prox = NULL;

    free(proxPonto);
    AdicionaCurvaNaInterface(Interface, curva);
    FinalizaCurvasCiclicas(Interface);
    return curva;
}

double FuncaoFilamento(double z, double R0, double Lz)
{
    return 1.0 - ( (1.0 - R0)*cos(M_PI*z/Lz) );
}

double FuncaoZero_ColisaoGotas(double Theta, double C_x, double C_y, double b)
{
    double angulo1 = - Theta;
    double angulo2 = M_PI - Theta;
    double U1[2] = {cos(angulo1), sin(angulo1)};
    double U2[2] = {cos(angulo2), sin(angulo2)};
    double Ur[2];

    // Fazendo Ur = U1 - U2
    Ur[0] = U1[0] - U2[0];
    Ur[1] = U1[1] - U2[1];

    // Fazendo Ur unitario
    double norma = sqrt( Ur[0]*Ur[0] + Ur[1]*Ur[1] );
    Ur[0] = Ur[0]/norma;
    Ur[1] = Ur[1]/norma;

    double dot_prod = C_x*Ur[0] + C_y*Ur[1];
    double vetor[2] = {C_x - Ur[0]*dot_prod, C_y - Ur[1]*dot_prod};

    norma = sqrt( vetor[0]*vetor[0] + vetor[1]*vetor[1] );
    return norma - b;
}

double Bisseccao_ColisaoGotas(double T1, double T2, double C_x, double C_y, double b)
{
    if( fabs(T2-T1)<1e-10 )
        return 0.5*(T2 + T1);

    double meio = 0.5*(T2 + T1);

    double valor1 = FuncaoZero_ColisaoGotas(T1, C_x, C_y, b);
    double valorMeio = FuncaoZero_ColisaoGotas(meio, C_x, C_y, b);

    if( fabs(valor1) < 1e-10 )
		return T1;

    if( valor1*valorMeio<=0 )
        return Bisseccao_ColisaoGotas(T1, meio, C_x, C_y, b);
    else
        return Bisseccao_ColisaoGotas(meio, T2, C_x, C_y, b);
}

///FUNCOES USADAS PARA UTILIZAR A ESTRUTURA DE DADOS DA INTERFACE
void AdicionaCurvaNaInterface(INTERFACE *Interface, CURVA *Curva)
{

    //Se  a interface ainda esta vazia (nenhuma curva)
    if( Interface->curvaInicial==NULL ) {
        Interface->curvaInicial = Curva;
        Curva->proxCurva = NULL;
        Curva->curvaAnt = NULL;
    }
    else { //Vamos fazer uma insercao no INICIO da lista de curvas
        Curva->curvaAnt = NULL;
        Curva->proxCurva = Interface->curvaInicial;
        Interface->curvaInicial->curvaAnt = Curva;
        Interface->curvaInicial = Curva;
    }

    (Interface->qtdCurvas)++;

    return;
}

REGIAO *CriaNovaRegiao(INTERFACE *Interface)
{
    CURVA *curva;
    REGIAO *novaRegiao;
    int maiorNumero = -1;

    for( curva=Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        if( curva->regiaoEsq!=NULL && curva->regiaoEsq->numero>maiorNumero )
            maiorNumero = curva->regiaoEsq->numero;
        if( curva->regiaoDir!=NULL && curva->regiaoDir->numero>maiorNumero )
            maiorNumero = curva->regiaoDir->numero;
    }

    novaRegiao = (REGIAO *)malloc( sizeof(REGIAO) );
    novaRegiao->numero = maiorNumero + 1;
    novaRegiao->status = INDEFINIDA;
    novaRegiao->localRegiao = LOCAL_INDEFINIDO;
    return novaRegiao;
}

void SetaTodasRegioesComoIndefinidas(INTERFACE *Interface)
{
    CURVA *curva;

    for( curva = Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        if( curva->regiaoEsq!=NULL )
            curva->regiaoEsq->status = INDEFINIDA;
        if( curva->regiaoDir!=NULL )
            curva->regiaoDir->status = INDEFINIDA;
    }

    return;
}

void InverteDirecaoCurva(CURVA *Curva)
{
    PONTO *ponto = NULL, *proxPonto = NULL, *temp = NULL;

    ponto = Curva->pontoInicial;
    do {
        proxPonto = ponto->prox;

        //Troca o ant pelo prox
        temp = ponto->ant;
        ponto->ant = ponto->prox;
        ponto->prox = temp;

        ponto = proxPonto;
    } while(ponto!=Curva->pontoInicial);

    return;
}

void JuntaCurvasComNodesIguais(INTERFACE *Interface)
{
    CURVA *curva, *curva2, *curvaTemp, *proxIteracao1, *proxIteracao2;
    PONTO *node1Curva1, *node2Curva1, *node1Curva2, *node2Curva2;

    curva=Interface->curvaInicial;
    while( curva!=NULL ) {
        proxIteracao1 = curva->proxCurva;

        curva2 = Interface->curvaInicial;
        while( curva2!=NULL ) {
            proxIteracao2 = curva2->proxCurva;

            if( curva==curva2 ) {
                curva2 = proxIteracao2;
                continue;
            }

            node1Curva1 = curva->pontoInicial;
            node2Curva1 = UltimoPonto(node1Curva1);
            node1Curva2 = curva2->pontoInicial;
            node2Curva2 = UltimoPonto(node1Curva2);

            if( PontosIguais( node1Curva1, node1Curva2 ) ) {
                InverteDirecaoCurva(curva);

                JuntaCurvasComNodesIguais(Interface);
                return;
            }
            else if( PontosIguais( node1Curva1, node2Curva2 ) ) {
                //Mantem a curva 2 e imenda ela na curva 1
                node2Curva2->prox = node1Curva1->prox;
                node1Curva1->prox->ant = node2Curva2;
                free(node1Curva1);

                //Remove a curva 1 da lista
                for( curvaTemp=Interface->curvaInicial; curvaTemp!=NULL; curvaTemp=curvaTemp->proxCurva ) {

                    //Se curva 1 for a primeira da lista
                    if( curvaTemp==curva2 ) {


                        Interface->curvaInicial = Interface->curvaInicial->proxCurva;
                        Interface->curvaInicial->curvaAnt = NULL;

                        if( proxIteracao2==curva )
                            DeuProblema("AAHAHAH INT\n");

                        free(curva);
                        break;
                    }

                    if( curvaTemp->proxCurva==curva ) {
                        curvaTemp->proxCurva = curva->proxCurva;
                        curva->proxCurva->curvaAnt = curvaTemp;

                        if( proxIteracao2==curva )
                            DeuProblema("AAHAHAH INT\n");

                        free(curva);
                        break;
                    }
                }

                JuntaCurvasComNodesIguais(Interface);
                return;
            }
            else if( PontosIguais( node2Curva1, node1Curva2 ) ) {

                //Mantem a curva 1 e imenda ela na curva 2
                node2Curva1->prox = node1Curva2->prox;
                node1Curva2->prox->ant = node2Curva1;
                free(node1Curva2);

                //Remove a curva 2 da lista
                for( curvaTemp=Interface->curvaInicial; curvaTemp!=NULL; curvaTemp=curvaTemp->proxCurva ) {

                    //Se curva2 for a primeira da lista
                    if( curvaTemp==curva2 ) {
                        Interface->curvaInicial = Interface->curvaInicial->proxCurva;
                        Interface->curvaInicial->curvaAnt = NULL;

                        if( proxIteracao1==curva2 )
                            proxIteracao1 = curva2->proxCurva;

                        //LiberaMemoriaCurva(curva2);
                        free(curva2);
                        break;
                    }

                    if( curvaTemp->proxCurva==curva2 ) {
                        curvaTemp->proxCurva = curva2->proxCurva;
                        curva2->proxCurva->curvaAnt = curvaTemp;

                        if( proxIteracao1==curva2 )
                            proxIteracao1 = curva2->proxCurva;

                        //LiberaMemoriaCurva(curva2);
                        free(curva2);
                        break;
                    }
                }

                JuntaCurvasComNodesIguais(Interface);
                return;
            }
            if( PontosIguais( node2Curva1, node2Curva2 ) ) {
                InverteDirecaoCurva(curva);
                JuntaCurvasComNodesIguais(Interface);
                return;
            }

            curva2 = proxIteracao2;
        }

        curva = proxIteracao1;
    }

    return;
}

void FinalizaCurvasCiclicas(INTERFACE *Interface)
{
    CURVA *curva;
    PONTO *ponto;

    for( curva=Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        ponto = UltimoPonto(curva->pontoInicial);
        if( PontosIguais(curva->pontoInicial, ponto) ) {

            //Garantindo que o x, y do primeiro e ultimo node vao ser iguais. Soh pre tirar errinhos minusculos de arredondamento
            ponto->x = curva->pontoInicial->x;
            ponto->y = curva->pontoInicial->y;

            ponto->prox = NULL;
        }
    }



}

void AtualizaLocalDasRegioes(INTERFACE *Interface)
{
    CURVA *curva;
    char continuar = 0;

    for( curva=Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        if( curva->regiaoEsq->localRegiao==INTERNA && curva->regiaoDir->localRegiao==LOCAL_INDEFINIDO ) {
            curva->regiaoDir->localRegiao = EXTERNA;
            continuar = 1;
        }

        if( curva->regiaoEsq->localRegiao==EXTERNA && curva->regiaoDir->localRegiao==LOCAL_INDEFINIDO ) {
            curva->regiaoDir->localRegiao = INTERNA;
            continuar = 1;
        }

        if( curva->regiaoDir->localRegiao==INTERNA && curva->regiaoEsq->localRegiao==LOCAL_INDEFINIDO ) {
            curva->regiaoEsq->localRegiao = EXTERNA;
            continuar = 1;
        }

        if( curva->regiaoDir->localRegiao==EXTERNA && curva->regiaoEsq->localRegiao==LOCAL_INDEFINIDO ) {
            curva->regiaoEsq->localRegiao = INTERNA;
            continuar = 1;
        }
    }

    if( continuar )
        AtualizaLocalDasRegioes(Interface);

    return;
}

void MudaNodeDaCurva(INTERFACE *Interface, double X, double Y)
{
    double distancia, menorDistancia;
    PONTO *pontoEscolhido = NULL, *ponto;
    CURVA *curva;

    for( curva=Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva ) {

        /// UMA GOTA
        pontoEscolhido = NULL;
        menorDistancia = 1e+10;
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;
            distancia = (ponto->x - X)*(ponto->x - X) + (ponto->y - Y)*(ponto->y - Y);

            if( distancia < menorDistancia ) {
                pontoEscolhido = ponto;
                menorDistancia = distancia;
            }
        }


        if( pontoEscolhido==NULL || pontoEscolhido->prox==NULL /*|| pontoEscolhido->ant*/) {
            continue;
        }



        if( pontoEscolhido!=curva->pontoInicial ) {
            ponto = UltimoPonto(curva->pontoInicial)->ant;
            ponto->prox = curva->pontoInicial;
            curva->pontoInicial->ant = ponto;
            curva->pontoInicial = pontoEscolhido;
            curva->pontoInicial->isNode = 1;
            curva->pontoInicial->ant->prox = CopiaPonto(curva->pontoInicial);
            curva->pontoInicial->ant->prox->isNode = 1;
            curva->pontoInicial->ant->prox->ant = curva->pontoInicial->ant;
        }

    }

}

// void AndaComONodeDaCurva(CURVA *curva, PONTO *P, int Passos)
// {
//     PONTO *pontoEscolhido;

//     if( P==NULL ) {
//         DeuProblema("AndaComONodeDaCurva: Argumento nulo 1\n");
//     }

//     /// === Se Passos < 0, vou usar Passos=Metade da quantidade de pontos da curva
//     if( Passos<0 ) {
//         int qtdPontos = 0;
//         LOOP_CURVA(c)
//             qtdPontos++;
//         Passos = qtdPontos/2;
//     }


//     for( pontoEscolhido=P; Passos>0; Passos-- )
//         pontoEscolhido = pontoEscolhido->prox;

//     if( pontoEscolhido!=curva->nodeInicial ) {
//         PONTO *ponto = UltimoPonto(curva->nodeInicial)->ant;
//         ponto->prox = curva->nodeInicial;
//         curva->nodeInicial->ant = ponto;
//         curva->nodeInicial = pontoEscolhido;
//         curva->nodeInicial->ant->prox = CopiaPonto(curva->nodeInicial);
//         curva->nodeInicial->ant->prox->isNode = 1;
//         curva->nodeInicial->ant->prox->ant = curva->nodeInicial->ant;
//         curva->nodeInicial->ant->prox->pontoProxCurva = curva->nodeInicial->prox;
//         curva->nodeInicial->pontoAntCurva = curva->nodeInicial->ant;
//         curva->nodeInicial->ant = NULL;
//     }

//     return;
// }

void RemovePontoDaCurva(CURVA *C, PONTO *P)
{
    // Caso degenerado em que a curva tem apenas um ponto
    if( C->pontoInicial->prox==C->pontoInicial ) {
        free(P);
        C->pontoInicial = NULL;
        return;
    }

    P->ant->prox = P->prox;
    P->prox->ant = P->ant;

    // Se for o primeiro ponto da curva
    if( C->pontoInicial == P )
        C->pontoInicial = P->prox;

    free(P);
}

PONTO *UltimoPonto(PONTO *P)
{
    return P;
}

char PontosIguais(PONTO *P1, PONTO *P2)
{
    if( fabs(P1->x - P2->x)<1e-10 && fabs(P1->y - P2->y)<1e-10 )
        return 1;
	else
		return 0;
    return 0;
}

char MenorAntiHorario(PONTO *A, PONTO *B, PONTO *Centro)
{
    double vec1x, vec1y, vec2x, vec2y, ang1, ang2;

	vec1x = A->x - Centro->x;
	vec1y = A->y - Centro->y;
	vec2x = B->x - Centro->x;
	vec2y = B->y - Centro->y;

	ang1 = atan2(vec1y, vec1x);
	ang2 = atan2(vec2y, vec2x);

	if( ang1>=0 && ang2>=0 )
		return( ang1 < ang2 );
	else if( ang1<0 && ang2<0 )
		return( ang1 < ang2 );
	else if( ang1>=0 && ang2<0 )
		return 1;
	else
		return 0;
	return 0;
}

char CurvaEntrando(CURVA *Curva, PONTO *Node)
{
	//O primeiro ponto da curva eh o node, esta saindo
	if( PontosIguais(Curva->pontoInicial, Node) )
		return 0;
	else if( PontosIguais(UltimoPonto(Curva->pontoInicial), Node) )
		return 1;
	else
		return -1; //Erro
    return -1;
}

char IsNodeNH(PONTO *Node, PONTO *NodesNH)
{
    for( ; NodesNH!=NULL; NodesNH=NodesNH->prox ) {
        if( PontosIguais(Node, NodesNH) )
            return 1;
    }

    return 0;
}

PONTO *CopiaPonto(PONTO *P)
{
    PONTO *novoPonto;

    novoPonto = (PONTO *)malloc( sizeof(PONTO) );
    novoPonto->ant = NULL;
    novoPonto->prox = NULL;
    novoPonto->isNH = 0;
    novoPonto->isNode = P->isNode;
    novoPonto->x = P->x;
    novoPonto->y = P->y;
    novoPonto->fixo = P->fixo;
    novoPonto->outflow = P->outflow;

    return novoPonto;
}


///FUNCOES USADAS NO ALGORITMO DE DESEMBARACAMENTO
// void RealizaDesembaracamento(INTERFACE *Interface, PONTO *PontosNH)
// {
//     PONTO *nodesNH = NULL, *ponto;
//     char continuar, continuar2;

//     nodesNH = PossuiEmbaracamentoInterface(Interface);
//     int i = 0;
//     for(ponto=nodesNH; ponto!=NULL; ponto=ponto->prox) {
// //        PrintDebug("NODES %lf %lf \n", ponto->x, ponto->y);

//         if( PontosNH!=NULL ) {
//             PontosNH[i].x = ponto->x;
//             PontosNH[i++].y = ponto->y;
//         }
//     }

//     if( nodesNH!=NULL ) {
// //        PrintDebug("COMECOU O ALGORITMO DE DESEMBARACAMENTO\n");

//         DivideSegmentosEmbaracados(Interface, nodesNH);
//         DefineComponentesCurvas(Interface, nodesNH);
//         SetaTodasRegioesComoIndefinidas(Interface);

//         continuar = 1;
//         continuar2 = 1;
//         while( continuar ) {
//             continuar = ProcuraRegioesFisicas(Interface, nodesNH);
//             continuar2 = ProcuraRegioesNaoFisicas(Interface, nodesNH);
//             if( continuar || continuar2 )
//                 continuar = 1;
//         }

//         RemoveRegioesNaoFisicas(Interface, nodesNH);
//         JuntaCurvasComNodesIguais(Interface);
//         FinalizaCurvasCiclicas(Interface);
//         AtualizaLocalDasRegioes(Interface);

//         PrintDebug("DESEMBARACOU!\n");
//     }

//     //Dando free nos nodesNH
//     while( nodesNH!=NULL ) {
//         ponto = nodesNH->prox;
//         free(nodesNH);
//         nodesNH = ponto;
//     }



// //    CURVA *curva;
// //    int qtd;
// //
// //    qtd = 0;
// //    for( curva=Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva )
// //        qtd++;
// ////    PrintDebug("CURVAS PROX: %d\n", qtd);
// //
// //    qtd = 0;
// //    for( curva=Interface->curvaInicial; curva->proxCurva!=NULL; curva=curva->proxCurva );
// //
// //    for( ; curva!=NULL; curva=curva->curvaAnt )
// //        qtd++;
// //    PrintDebug("CURVAS ANT: %d\n", qtd);



// }

char PossuiInterseccaoSegmentos(PONTO *Seg1Ponto1, PONTO *Seg1Ponto2, PONTO *Seg2Ponto1, PONTO *Seg2Ponto2, double *X, double *Y)
{
    double x1, y1, x1L, y1L, x2, y2, x2L, y2L;
    double det;
    double lambda[2];

    //Definindo o segmento 1
	x1 = Seg1Ponto1->x; y1 = Seg1Ponto1->y;
	x1L = Seg1Ponto2->x; y1L = Seg1Ponto2->y;

	//Definindo o segmento 2
	x2 = Seg2Ponto1->x; y2 = Seg2Ponto1->y;
	x2L = Seg2Ponto2->x; y2L = Seg2Ponto2->y;



    //Construindo e resolvendo o sistema de interseccao
	double matriz[2][2] = { {x1L-x1, -(x2L-x2)}, {y1L-y1, -(y2L-y2)} };
	double vet[2] = {x2-x1, y2-y1};

	det = matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0];
	if( fabs(det)<1e-14 ) {
		*X = *Y = 0;
		return 0; //Nao encontrou, as retas eram paralelas
	}
	lambda[1] = (vet[1]*matriz[0][0] - vet[0]*matriz[1][0])/det;
	lambda[0] = (vet[0]*matriz[1][1] - vet[1]*matriz[0][1])/det;

	//Ponto de interseccao
	*X = x1 + lambda[0]*(x1L-x1);
	*Y = y1 + lambda[0]*(y1L-y1);

	if( lambda[0]<0.0 || lambda[0]>1.0 || lambda[1]<0.0 || lambda[1]>1.0 )
		return 0; //nao possui interseccao
    else if( fabs(lambda[0] - 0)<1e-10 || fabs(lambda[0] - 1)<1e-10 || fabs(lambda[1] - 0)<1e-10 || fabs(lambda[1] - 1)<1e-10 )
        return 1; //A interseccao ocorre nos extremos nos segmentos, tlvez precise tratar isso
	else {
        //printf("LAMBDA %.15lf \n", lambda[0]);
//        printf("ENCONTROU INTERSEC \n%.15lf %.15lf %.15lf %.15lf", Seg1Ponto1->x, Seg1Ponto1->y, Seg1Ponto2->x, Seg1Ponto2->y);
//        printf("\n%.15lf %.15lf %.15lf %.15lf\n", Seg2Ponto1->x, Seg2Ponto1->y, Seg2Ponto2->x, Seg2Ponto2->y);
//        printf("lambda %.15lf %.15lf\n", lambda[0], lambda[1]);
		return 1; //Encontrou interseccao
	}

    return 0;
}

// PONTO *PossuiEmbaracamentoInterface(INTERFACE *Interface)
// {
//     CURVA *curva1, *curva2;
//     PONTO *seg1P1, *seg1P2, *seg2P1, *seg2P2;
//     PONTO *nodeNH = NULL, *novoNode, *node;
//     double x, y;
//     char possui;
//     /*static aaa = 1;

//     printf("\n AAA: %d \n", aaa++);
//     getchar();*/

//     //Percorre todas as curvas
//     for( curva1=Interface->curvaInicial; curva1!=NULL; curva1=curva1->proxCurva ) {

//         //Percorre cada segmento desta curva
//         LOOP_CURVA(curva1) {
//             seg1P1 = loop_stuff.ponto;
//             seg1P2 = seg1P1->prox;

//             //Percorre todas as curvas
//             for( curva2=Interface->curvaInicial; curva2!=NULL; curva2=curva2->proxCurva ) {
//                 //Percorre cada segmento desta curva
//                 LOOP_CURVA(curva2) {
//                     seg2P1 = loop_stuff.ponto;
//                     seg2P2 = seg2P1->prox;

//                     //Ignora se os dois segmentos forem consecutivos numa mesma curva
//                     if( seg2P1==seg1P2 || seg1P1==seg2P2 )
//                         continue;
//                     else if( PontosIguais(seg2P1, seg1P2) || PontosIguais(seg1P1, seg2P2) ) //Apenas caso tente comparar o ultimo segmento de uma curva ciclica com o primeiro
//                         continue;

//                     possui = PossuiInterseccaoSegmentos(seg1P1, seg1P2, seg2P1, seg2P2, &x, &y);

//                     //Se encontrou cria um novo node NH pra cada um encontrado
//                     if( possui ) {
//                         novoNode = (PONTO *)malloc( sizeof(PONTO) );
//                         novoNode->isNode = 1;
//                         novoNode->isNH = 1;
//                         novoNode->fixo = 0;
//                         novoNode->outflow = 0;
//                         novoNode->ant = novoNode->prox = NULL;
//                         novoNode->x = x;
//                         novoNode->y = y;
//                         novoNode->curvas[0] = curva1;
//                         novoNode->curvas[1] = curva2;
//                         novoNode->seg1P1 = seg1P1;
//                         novoNode->seg2P1 = seg2P1;
//                         novoNode->curvaEntrando[0] = ENTRANDO;
//                         novoNode->curvaEntrando[1] = ENTRANDO;
//                         novoNode->curvaEntrando[2] = TIPO_INDEFINIDO;
//                         novoNode->curvaEntrando[3] = TIPO_INDEFINIDO;

//                         //Apenas verifica se esse ja foi encontrado
//                         for( node=nodeNH; node!=NULL; node=node->prox ) {
//                             if( PontosIguais(novoNode, node) )
//                                 break;
//                         }
//                         if( node==NULL ) { //Nao encontrou, insere no inicio da lista
//                             novoNode->prox = nodeNH;
//                             nodeNH = novoNode;
//                         }
//                         else
//                             free(novoNode);
//                     }
//                 }
//             }
//         }
//     }


//     return nodeNH;
// }

// void DivideSegmentosEmbaracados(INTERFACE *Interface, PONTO *NodesNH)
// {
//     PONTO *node = NodesNH, *node2, *ponto;
//     CURVA *novaCurva;
//     int j = 0;

//     //Percorre cada nodeNH encontrado
//     for( node=NodesNH; node!=NULL; node=node->prox ) {

//         //Vendo as curvas que cruzaram pra dar esse node
//         /*novaCurva = (CURVA *)node->curvas[0];
//         printf("CRUVA CRUZOU 1: [%lf %lf] [%lf %lf]\n", novaCurva->nodeInicial->x, novaCurva->nodeInicial->y, novaCurva->nodeInicial->prox->x, novaCurva->nodeInicial->prox->y);
//         novaCurva = (CURVA *)node->curvas[1];
//         printf("CRUVA CRUZOU 1: [%lf %lf] [%lf %lf]\n", novaCurva->nodeInicial->x, novaCurva->nodeInicial->y, novaCurva->nodeInicial->prox->x, novaCurva->nodeInicial->prox->y);
//         */

//         //Quebrando a primeira curva em duas curvas
//         novaCurva = (CURVA *)malloc( sizeof(CURVA) );
//         novaCurva->regiaoEsq = novaCurva->regiaoDir = NULL;
//         novaCurva->curvaAnt = NULL; //Insercao no inicio da lista de curvas
//         novaCurva->proxCurva = Interface->curvaInicial;
//         novaCurva->pontoInicial = CopiaPonto(node);
//         novaCurva->pontoInicial->prox = node->seg1P1->prox; //Imenda o pedaco da curva anterior nessa
//         novaCurva->pontoInicial->prox->ant = novaCurva->pontoInicial; //ADICIONEI AGORAAAAA

//         AdicionaCurvaNaInterface(Interface, novaCurva);
//         node->seg1P1->prox = CopiaPonto(node); //Finalizando a curva anterior nesse node
//         node->seg1P1->prox->ant = node->seg1P1;
//         node->curvas[2] = novaCurva;
//         node->curvaEntrando[2] = SAINDO;


//         //PlotInterfaceGnuplotHA(*Interface);

//         //Ajeitando os ponteiros das curvas se precisar
//         for( node2=NodesNH; node2!=NULL; node2=node2->prox ) {
//             j++;

//             //Verifica se precisa arrumar a curva0 (pelo seg1P1)
//             if( node2!=node && node2->curvas[0]==node->curvas[0] ) {//Pode precisar arrumar
//                 //Vamos ver se esse segmento ficou um pedaco na curva anterior e outro pedaco na nova
//                 ponto = UltimoPonto( ((CURVA *)node->curvas[0])->nodeInicial );
//                 ponto = ponto->ant; //penultimo ponto da curva anterior (o ultimo eh o node NH)
//                 if( ponto == node2->seg1P1 ) {
//                     node2->seg1P1 = novaCurva->pontoInicial;
//                     node2->curvas[0] = novaCurva;
//                 }

//                 //Verifica se esse segmento foi parar na curva nova
//                 for( ponto=novaCurva->nodeInicial; ponto!=NULL; ponto=ponto->prox ) {
//                     if( node2->seg1P1==ponto ) {
//                         node2->curvas[0] = novaCurva;
//                     }
//                 }
//             }

//             //Verifica se precisa arrumar a curva1 (pelo seg2P1)
//             if( node2->curvas[1]==node->curvas[0] ) {//Pode precisar arrumar
//                 //Vamo ver se esse segmento ficou um pedaco na curva anterior e outro pedaco na nova
//                 ponto = UltimoPonto( ((CURVA *)node->curvas[0])->pontoInicial );
//                 ponto = ponto->ant; //penultimo ponto da curva anterior (o ultimo eh o node NH)
//                 if( ponto == node2->seg2P1 ) {
//                     node2->seg2P1 = novaCurva->pontoInicial;
//                     node2->curvas[1] = novaCurva;
//                 }


//                 //Verifica se esse segmento foi parar na curva nova
//                 for( ponto=novaCurva->nodeInicial; ponto!=NULL; ponto=ponto->prox ) {
//                     if( node2->seg2P1==ponto )
//                         node2->curvas[1] = novaCurva;
//                 }
//             }

//         }


//         //Quebrando a segunda curva em duas curvas
//         novaCurva = (CURVA *)malloc( sizeof(CURVA) );
//         novaCurva->regiaoEsq = novaCurva->regiaoDir = NULL;
//         novaCurva->curvaAnt = NULL; //Insercao no inicio da lista de curvas
//         novaCurva->proxCurva = Interface->curvaInicial;
//         novaCurva->nodeInicial = CopiaPonto(node);
//         novaCurva->nodeInicial->prox = node->seg2P1->prox; //Imenda o pedaco da curva anterior nessa
//         novaCurva->nodeInicial->prox->ant = novaCurva->nodeInicial;//ADICIONEI AGORA;

//         AdicionaCurvaNaInterface(Interface, novaCurva);
//         node->seg2P1->prox = CopiaPonto(node); //Finalizando a curva anterior nesse node
//         node->seg2P1->prox->ant = node->seg2P1;
//         node->curvas[3] = novaCurva;
//         node->curvaEntrando[3] = SAINDO;

//         //PlotInterfaceGnuplotHA(*Interface);

//         //Ajeitando os ponteiros das curvas se precisar
//         for( node2=NodesNH; node2!=NULL; node2=node2->prox ) {
//             //Verifica se precisa arrumar a curva0 (pelo seg1P1)
//             if( node2->curvas[0]==node->curvas[1] ) {//Pode precisar arrumar
//                 //Vamo ver se esse segmento ficou um pedaco na curva anterior e outro pedaco na nova
//                 ponto = UltimoPonto( ((CURVA *)node->curvas[1])->nodeInicial );
//                 ponto = ponto->ant; //penultimo ponto da curva anterior (o ultimo eh o node NH)
//                 if( ponto == node2->seg1P1 ) {

//                     node2->seg1P1 = novaCurva->nodeInicial;
//                     node2->curvas[0] = novaCurva;
//                 }

//                 //Verifica se esse segmento foi parar na curva nova
//                 for( ponto=novaCurva->nodeInicial; ponto!=NULL; ponto=ponto->prox ) {
//                     if( node2->seg1P1==ponto )
//                         node2->curvas[0] = novaCurva;
//                 }
//             }

//             //Verifica se precisa arrumar a curva1 (pelo seg2P1)
//             if( node2!=node && node2->curvas[1]==node->curvas[1] ) {//Pode precisar arrumar
//                 //Vamo ver se esse segmento ficou um pedaco na curva anterior e outro pedaco na nova
//                 ponto = UltimoPonto( ((CURVA *)node->curvas[1])->nodeInicial );
//                 ponto = ponto->ant; //penultimo ponto da curva anterior (o ultimo eh o node NH)
//                 if( ponto == node2->seg2P1 ) {
//                     node2->seg2P1 = novaCurva->nodeInicial;
//                     node2->curvas[1] = novaCurva;
//                 }

//                 //Verifica se esse segmento foi parar na curva nova
//                 for( ponto=novaCurva->nodeInicial; ponto!=NULL; ponto=ponto->prox ) {
//                     if( node2->seg2P1==ponto )
//                         node2->curvas[1] = novaCurva;
//                 }
//             }
//         }

//     }

//     for( node=NodesNH; node!=NULL; node=node->prox )
//         OrdenaCurvasAntiHorario(node);

//     return;
// }

// void OrdenaCurvasAntiHorario(PONTO *NodeH)
// {
//     PONTO *ponto[4], *node;
//     CURVA *curva;
//     TIPO_CURVA tipoCurva;
//     int i, j;

//     //As curvas 0 e 1 sao de entrada, as 2 e 3 de saida
//     curva = (CURVA *)NodeH->curvas[0];
//     for( ponto[0]=curva->nodeInicial; ponto[0]->prox!=NULL; ponto[0]=ponto[0]->prox ); //Anda ate o final da curva
//     ponto[0] = ponto[0]->ant;

//     curva = (CURVA *)NodeH->curvas[1];
//     for( ponto[1]=curva->nodeInicial; ponto[1]->prox!=NULL; ponto[1]=ponto[1]->prox ); //Anda ate o final da curva
//     ponto[1] = ponto[1]->ant;

//     curva = (CURVA *)NodeH->curvas[2];
//     ponto[2] = curva->nodeInicial->prox;

//     curva = (CURVA *)NodeH->curvas[3];
//     ponto[3] = curva->nodeInicial->prox;


// 	//Insertion short
// 	for( i=0; i<4; i++ ) {
// 		for( j=i+1; j<4; j++ ) {

// 			if( MenorAntiHorario( ponto[j], ponto[i], NodeH ) ) { //Troca
// 				node = ponto[i];
// 				ponto[i] = ponto[j];
// 				ponto[j] = node;

// 				curva = NodeH->curvas[i];
// 				NodeH->curvas[i] = NodeH->curvas[j];
// 				NodeH->curvas[j] = curva;

// 				tipoCurva = NodeH->curvaEntrando[i];
// 				NodeH->curvaEntrando[i] = NodeH->curvaEntrando[j];
// 				NodeH->curvaEntrando[j] = tipoCurva;
//             }

// 		}
// 	}

// 	return;
// }

// void DefineComponentesCurvas(INTERFACE *Interface, PONTO *NodesNH)
// {
//     char continuar, entrando;
//     int i, iMais1, iMenos1;
//     PONTO *node;
//     CURVA *curva, *curvaAnt, *curvaProx;
//     REGIAO *regiao;

//     /*printf("\n\n REGIOES\n");
//                 for( curva = Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//                     int esq, dir;
//                     esq = curva->regiaoEsq!=NULL ? curva->regiaoEsq->numero : -1;
//                     dir = curva->regiaoDir!=NULL ? curva->regiaoDir->numero : -1;
//                     printf("CURVA [%lf %lf] [%lf %lf] Esq: %d Dir: %d\n", curva->nodeInicial->x, curva->nodeInicial->y, curva->nodeInicial->prox->x, curva->nodeInicial->prox->y, esq, dir);
//                 }*/

//     continuar = 1;
//     while( continuar ) {
//         continuar = 0;

//         //Percorre cada node NH
//         for( node=NodesNH; node!=NULL; node=node->prox ) {

//             //Percorre as quatro curvas deste node
//             for( i=0; i<4; i++ ) {
//                 curva = node->curvas[i];
//                 curvaProx = node->curvas[ (i+1)%4 ];
//                 iMais1 = (i+1)%4;
//                 if( i==0 ) {
//                     curvaAnt = node->curvas[3];
//                     iMenos1 = 3;
//                 }
//                 else {
//                     curvaAnt = node->curvas[i-1];
//                     iMenos1 = i-1;
//                 }

//                 //Primeiro vamos usar o valor de curva para preencher as componentes da curvaAnt
//                 regiao = NULL;

//                 //Verifica se curva esta saindo ou entrando
//                 entrando = node->curvaEntrando[i];//CurvaEntrando(curva, node);
//                 if( entrando==ENTRANDO ) //Entrando
// 					regiao = curva->regiaoEsq;
// 				else if( entrando==SAINDO )
// 					regiao = curva->regiaoDir;
// 				else
// 					printf("ERRO: CASO INESPERADO 1 %d\n", i);


//                 //Preenche o dir da curvaAnt
//                 entrando = node->curvaEntrando[iMenos1];//CurvaEntrando(curvaAnt, node);
//                 if( entrando==ENTRANDO && curvaAnt->regiaoDir==NULL && regiao!=NULL) {
//                     curvaAnt->regiaoDir = regiao;
//                     continuar = 1;
//                 }
//                 else if( entrando==SAINDO && curvaAnt->regiaoEsq==NULL && regiao!=NULL ) {
//                     curvaAnt->regiaoEsq = regiao;
//                     continuar = 1;
//                 }
//                 else if( entrando==TIPO_INDEFINIDO )
//                     printf("ERRO: CASO INESPERADO 2 %d\n", i);


//                 //Agora vamos usar o valor de curva para preencher as componentes de curvaProx
//                 regiao = NULL;

//                 //Verifica se curva esta saindo ou entrando
//                 entrando = node->curvaEntrando[i];//CurvaEntrando(curva, node);
//                 if( entrando==ENTRANDO ) //Entrando
// 					regiao = curva->regiaoDir;
// 				else if( entrando==SAINDO )
// 					regiao = curva->regiaoEsq;
// 				else
// 					printf("ERRO: CASO INESPERADO 3 %d\n", i);

//                 //Preenche o esq na curvaProx
//                 entrando = node->curvaEntrando[iMais1];//CurvaEntrando(curvaProx, node);
//                 if( entrando==ENTRANDO && curvaProx->regiaoEsq==NULL && regiao!=NULL ) {
//                     curvaProx->regiaoEsq = regiao;
//                     continuar = 1;
//                 }
//                 else if( entrando==SAINDO && curvaProx->regiaoDir==NULL && regiao!=NULL ) {
//                     curvaProx->regiaoDir = regiao;
//                     continuar = 1;
//                 }
//                 else if( entrando==TIPO_INDEFINIDO )
//                     printf("ERRO: CASO INESPERADO 4 %d %d\n", i, CurvaEntrando(curvaProx, node));




//                 /*printf("\n\n REGIOES\n");
//                 for( curva = Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//                     int esq, dir;
//                     esq = curva->regiaoEsq!=NULL ? curva->regiaoEsq->numero : -1;
//                     dir = curva->regiaoDir!=NULL ? curva->regiaoDir->numero : -1;
//                     printf("CURVA [%lf %lf] [%lf %lf] Esq: %d Dir: %d\n", curva->nodeInicial->x, curva->nodeInicial->y, curva->nodeInicial->prox->x, curva->nodeInicial->prox->y, esq, dir);
//                 }*/
//             }

//         }

//     }

//     //Se tiver ficado alguma curva sem regiao no Esq ou Dir, vamos criar uma nova regiao e chamar essa funcao de novo recursivamente
//     continuar = PreencheUmaNovaRegiao(Interface);
//     if( continuar )
//         DefineComponentesCurvas(Interface, NodesNH);

//     return;
// }

// char PreencheUmaNovaRegiao(INTERFACE *Interface)
// {
//     CURVA *curva;

//     for( curva=Interface->curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//         if( curva->regiaoEsq==NULL ) {
//             curva->regiaoEsq = CriaNovaRegiao(Interface);
//             return 1;
//         }
//         else if( curva->regiaoDir==NULL ) {
//             curva->regiaoDir = CriaNovaRegiao(Interface);
//             return 1;
//         }
//     }

//     return 0;
// }

// char ProcuraRegioesFisicas(INTERFACE *Interface, PONTO *NodesNH)
// {
//     PONTO *node, *node1 = NULL, *node2 = NULL;
//     CURVA *curva1, *curva2;
//     REGIAO *regiao = NULL;
//     int i, iMais1;
//     char entrando, encontrouNaoFisica;
//     int modificouAlgo = 0;


//     //Percorre cada node NH
//     for( node=NodesNH; node!=NULL; node=node->prox ) {
//         //Indica se encontrou alguma regiao ja marcada como nao-fisica ao redor desse node
//         //Se tiver alguma, todas as outras 3 sao fisicas
//         encontrouNaoFisica = 0;

//         //Percorre as 4 regioes ao redor deste node
//         for( i=0; i<4; i++ ) {

//             curva1 = node->curvas[i];
//             if( i!=3 ) {
//                 curva2 = node->curvas[i+1];
//                 iMais1 = i+1;
//             }
//             else {
//                 curva2 = node->curvas[0];
//                 iMais1 = 0;
//             }

//             //Encontrando o outro no da curva 1
//             entrando = node->curvaEntrando[i];
//             if( entrando==ENTRANDO ) {
// 				node1 = curva1->nodeInicial;
// 				regiao = curva1->regiaoDir;
//             }
// 			else if( entrando==SAINDO ) {
// 				node1 = UltimoPonto(curva1->nodeInicial);
// 				regiao = curva1->regiaoEsq;
// 			}
// 			else
// 				DeuProblema("ERRO: caso inesperado.\n");


//             //Encontrando o outro no da curva 1
//             entrando = node->curvaEntrando[iMais1];
//             if( entrando==ENTRANDO )
// 				node2 = curva2->nodeInicial;
// 			else if( entrando==SAINDO )
// 				node2 = UltimoPonto(curva2->nodeInicial);
// 			else
// 				DeuProblema("ERRO: caso inesperado.\n");


//             //Verifica a propriedade 3 do artigo, se no1 ou no2 foram nao-NH entao essa componente eh fisica
//             if( !IsNodeNH(node1, NodesNH) || !IsNodeNH(node2, NodesNH)  ) {
//                 if( regiao->status!=FISICA ) {
//                     regiao->status = FISICA;
//                     modificouAlgo = 1;
//                 }
//             }



//             //Agora vamos verificar se alguma das regioes em volta dessa curva ja foi marcada como nao-fisica
//             if( curva1->regiaoEsq->status==NAO_FISICA || curva1->regiaoDir->status==NAO_FISICA )
//                 encontrouNaoFisica = 1;
//         }

//         //Se tiver achado alguma nao-fisica, marca todas as outras como fisica
//         if( encontrouNaoFisica ) {
//             for( i=0; i<4; i++ ) {
//                 curva1 = node->curvas[i];
//                 if( curva1->regiaoEsq->status==INDEFINIDA ) {
//                     curva1->regiaoEsq->status = FISICA;
//                     modificouAlgo = 1;
//                 }
//                 if( curva1->regiaoDir->status==INDEFINIDA ) {
//                     curva1->regiaoDir->status = FISICA;
//                     modificouAlgo = 1;
//                 }
//             }
//         }
//     }

//     return modificouAlgo;
// }

// char ProcuraRegioesNaoFisicas(INTERFACE *Interface, PONTO *NodesNH)
// {
//     PONTO *node = NULL;
//     CURVA *curva = NULL;
//     REGIAO *regiao = NULL, *regiaoNaoClassificada = NULL;
//     int qtdFisica, i;
//     char entrando, modificouAlgo = 0;

//     //Percorre cada node NH
//     for( node=NodesNH; node!=NULL; node=node->prox ) {

//         qtdFisica = 0;
//         regiaoNaoClassificada = NULL;

//         //Passa por cada uma das 4 regioes desse node
//         for( i=0; i<4; i++ ) {
//             curva = node->curvas[i];

//             //Pega o ponteiro pra essa regiao
//             entrando = node->curvaEntrando[i];
//             if( entrando==ENTRANDO )
// 				regiao = curva->regiaoDir;
// 			else if( entrando==SAINDO )
// 				regiao = curva->regiaoEsq;
// 			else
// 				DeuProblema("ERRO: caso inesperado\n");

//             //Verifica se essa regiao ja foi classificada como fisica
//             if( regiao->status==FISICA )
//                 qtdFisica++;
//             else if( regiao->status==INDEFINIDA )
//                 regiaoNaoClassificada = regiao;
//         }

//         //Se tiver encontrado 3 fisicas, a que sobrou eh nao-fisica
//         if( qtdFisica==3 && regiaoNaoClassificada!=NULL ) {
//             regiaoNaoClassificada->status = NAO_FISICA;
//             modificouAlgo = 1;
//         }
//     }

//     return modificouAlgo;
// }

// void RemoveRegioesNaoFisicas(INTERFACE *Interface, PONTO *NodesNH)
// {
//     PONTO *node;
//     REGIAO *regiao, *regiaoFisica1, *regiaoFisica2;
//     CURVA *curva, *curvaProx, *proxIteracao, *curvaTemp, *curvaAnt;
//     CURVA *curvasRemover[200]; ///Tirar esse tamanho estatico depois... vai que passa de 200 neh... (improvavel)
//     int i, qtdRemover;
//     char entrando;

//     qtdRemover = 0;

//     //Percorre cada node
//     for( node=NodesNH; node!=NULL; node=node->prox ) {
//         //Percorre as 4 curvas desse node
//         for( i=0; i<4; i++ ) {
//             curva = node->curvas[i];

//             //Pega o ponteiro dessa regiao
//             regiao = NULL;
//             regiaoFisica1 = NULL;
//             entrando = node->curvaEntrando[i];
//             if( entrando==ENTRANDO ) {
// 				regiao = curva->regiaoDir;
//                 regiaoFisica1 = curva->regiaoEsq;
//             }
// 			else if( entrando==SAINDO ) {
// 				regiao = curva->regiaoEsq;
//                 regiaoFisica1 = curva->regiaoDir;
// 			}
// 			else
// 				printf("ERRO: caso inesperado\n");

// 			//Se for uma regiao nao fisica, vai remover esta curva e a proxima (i+1)
//             if( regiao->status==NAO_FISICA ) {
//                 regiao->status = FISICA;

//                 //Pegando o ponteiro pra regiao fisica que cerca a curva (i+1)
//                 regiaoFisica2 = NULL;
//                 curvaProx = node->curvas[ (i+1)%4 ];
//                 entrando = node->curvaEntrando[ (i+1)%4 ];
//                 if( entrando==ENTRANDO )
//                     regiaoFisica2 = curvaProx->regiaoDir;
//                 else if( entrando==SAINDO )
//                     regiaoFisica2 = curvaProx->regiaoEsq;
//                 else
//                     printf("ERRO: caso inesperado\n");

//                 //Removendo as curvas
//                 //percorre todas as curvas da interface e troca a regiaoFisica1 e regiaoFisica2 pela regiao
//                 for( curvaTemp=Interface->curvaInicial; curvaTemp!=NULL; curvaTemp=curvaTemp->proxCurva ) {
//                     if( curvaTemp->regiaoEsq==regiaoFisica1 || curvaTemp->regiaoEsq==regiaoFisica2 )
//                         curvaTemp->regiaoEsq = regiao;
//                     if( curvaTemp->regiaoDir==regiaoFisica1 || curvaTemp->regiaoDir==regiaoFisica2 )
//                         curvaTemp->regiaoDir = regiao;
//                 }

//                 free(regiaoFisica1);
//                 if( regiaoFisica2!=regiaoFisica1 )
//                     free(regiaoFisica2);


//                 //Agora remove as curvas da lista de curvas da interface
//                 curvaAnt = NULL;
//                 curvaTemp = Interface->curvaInicial;
//                 while( curvaTemp!=NULL ) {
//                     proxIteracao = curvaTemp->proxCurva;

//                     if( curvaTemp==curva || curvaTemp==curvaProx ) {
//                         //Se for a primeira curva da interface
//                         if( curvaTemp==Interface->curvaInicial ) {
//                             Interface->curvaInicial = curvaTemp->proxCurva;
//                             Interface->curvaInicial->curvaAnt = NULL;
//                             curvasRemover[qtdRemover++] = curvaTemp;
//                         }
//                         else {
//                             curvaAnt->proxCurva = curvaTemp->proxCurva;
//                             curvaTemp->proxCurva->curvaAnt = curvaAnt;
//                             curvasRemover[qtdRemover++] = curvaTemp;
//                         }
//                     }
//                     else
//                         curvaAnt = curvaTemp;

//                     curvaTemp = proxIteracao;
//                 }
//             }
//         }
//     }


//     for( i=0; i<qtdRemover; i++ )
//         LiberaMemoriaCurva(curvasRemover[i]);

//     return;
// }

void LiberaMemoriaCurva(CURVA *C)
{
    PONTO *ponto, *proxPonto;

    ponto = C->pontoInicial;
    while( ponto!=NULL ) {
        proxPonto = ponto->prox;
        free(ponto);
        ponto = proxPonto;
    }

    free(C);
    return;
}




