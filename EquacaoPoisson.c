#include "EquacaoPoisson.h"

extern int passoTemp;

extern MALHA malha_temp;

extern int velocidadeNewtoniana;

extern double **debugsCell;

PetscErrorCode InicializaMatrizPoisson(Mat A, int NumLinhas, MALHA M, int iStart, int iEnd)
{
    int linha, indiceVetor = 0, coluna;
    int i, j;
    int d_nnz[iEnd - iStart], o_nnz[iEnd - iStart], qtd[iEnd - iStart];
    PetscErrorCode ierr;

    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesP[linha].i;
        j = M.indicesP[linha].j;

        o_nnz[indiceVetor] = 0;
        //Diagonal principal sempre vai ter
        qtd[indiceVetor] = 1;
        d_nnz[indiceVetor] = 1;


        //Caso surface
        if( M.pontosP[i][j].tipo==SURFACE ) {
            AlocaLinhaCondicaoContornoPoisson(&M, i, j, indiceVetor, iStart, iEnd, d_nnz, o_nnz, qtd);
            indiceVetor++;
            continue;
        }
//        if( M.pontosP[i][j].tipo==SURFACE ) {
//            qtd[indiceVetor] = 1;
//            d_nnz[indiceVetor] = 1;
//            o_nnz[indiceVetor] = 0;
//            indiceVetor++;
//            continue;
//        }


        /// Verificando se eh condicao de periodicidade
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE ||
           M.pontosU[i+1][j].tipoBoundary==PERIODICIDADE ) {
            qtd[indiceVetor]++;
            o_nnz[indiceVetor]++;
        }


        //Verificando o da esquerda
        if( i!=0 && (M.pontosP[i-1][j].tipo==FULL || M.pontosP[i-1][j].tipo==SURFACE) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesP[i-1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }

        //Verificando o de baixo
        if( j!=0 && (M.pontosP[i][j-1].tipo==FULL || M.pontosP[i][j-1].tipo==SURFACE) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesP[i][j-1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }

        //Verificando o da direita
        if( i!=M.Nx-1 && (M.pontosP[i+1][j].tipo==FULL || M.pontosP[i+1][j].tipo==SURFACE) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesP[i+1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }

        //Verificando o de cima
        if( j!=M.Ny-1 && (M.pontosP[i][j+1].tipo==FULL || M.pontosP[i][j+1].tipo==SURFACE) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesP[i][j+1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }

        indiceVetor++;
    }

    ierr = MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A, 0, qtd); CHKERRQ(ierr);
    ierr = MatSeqSBAIJSetPreallocation(A, 1, 0, qtd); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode MontaMatrizPoisson(Mat A, MALHA M, double **ViscosityNewtGen)
{
    int i, j;
    int qtd;
    int linha=0, colunas[5], centro;
    double valores[5], acrescentaCentro;
    double coefEsq, coefDir, coefBaixo, coefCima, coefCentro;
    double hx1, hx2, hy1, hy2;
    double cilindrico;
    int iStart, iEnd;

    MatGetOwnershipRange(A, &iStart, &iEnd);

    cilindrico = (M.tipoCoord==AXI_CILINDRICO) ? 1.0 : 0.0;

    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesP[linha].i;
        j = M.indicesP[linha].j;

        qtd = 0;
        valores[0]=valores[1]=valores[2]=valores[3]=valores[4]=0.0;
        acrescentaCentro = 0.0;

        //Caso de condicao de contorno, vou chamar a outra funcao que coloca essa linha da matriz
        if( M.pontosP[i][j].tipo==SURFACE ) {

            // Formulacao tradicional. Com equacao normal de sup livre
            LinhaCondicaoContornoPoisson(A, &M, ViscosityNewtGen, i, j, linha);

            // Formulacao laplaciano superficial explicita - debora
            //colunas[0] = linha;
            //valores[0] = 1.0;
            //qtd = 1;
            //MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES);

            // Formulacao laplaciano superficial implicita - debora
            //LaplacianoSupImplicitoMatriz(A, &M, i, j, linha);

            continue;
        }

        hx1 = (i==0) ? M.dx[i] : 0.5*(M.dx[i-1] + M.dx[i]);
        hx2 = (i==M.Nx-1) ? M.dx[i] : 0.5*(M.dx[i] + M.dx[i+1]);
        hy1 = (j==0) ? M.dy[j] : 0.5*(M.dy[j-1] + M.dy[j]);
        hy2 = (j==M.Ny-1) ? M.dy[j] : 0.5*(M.dy[j] + M.dy[j+1]);

        coefEsq = ( (-cilindrico*hx2/M.r[i]) + 2.0 )/(hx1*(hx1+hx2));
        coefDir = ( (cilindrico*hx1/M.r[i]) + 2.0 )/(hx2*(hx1+hx2));
        coefBaixo = 2.0/(hy1*(hy1+hy2));
        coefCima = 2.0/(hy2*(hy1+hy2));
        coefCentro = ( ((cilindrico*(hx2-hx1)/M.r[i]) - 2.0)/(hx1*hx2) ) - (2.0/(hy1*hy2));

        /// === Primeiramente verificando se ta na direita e eh condicao de PERIODICIDADE
        /// === Coloca um termo do lado esquerdo do dominio, caso seja
        if( M.pontosU[i+1][j].tipoBoundary == PERIODICIDADE ) {
            valores[qtd] = coefDir;
            colunas[qtd++] = M.matIndicesP[0][j];
        }


        //Colocando os coeficientes do Psi(i-1, j)
        if( i!=0 && (M.pontosP[i-1][j].tipo==FULL || M.pontosP[i-1][j].tipo==SURFACE) )
        {
            valores[qtd] = coefEsq;
            colunas[qtd++] = M.matIndicesP[i-1][j];
        }
        else {
            if( M.pontosU[i][j].tipoBoundary == NEUMANN )
                acrescentaCentro += -coefEsq; //Dirichlet na velocidade vira neumann na Psi
            else if( M.pontosU[i][j].tipoBoundary == NOSLIP )
                acrescentaCentro += coefEsq; //Dirichlet na velocidade vira neumann na Psi
            else if( M.pontosU[i][j].tipoBoundary == SLIP )
                acrescentaCentro += coefEsq; //Dirichlet na velocidade vira neumann na Psi
            else if( M.pontosU[i][j].tipoBoundary == SIMETRIA )
                acrescentaCentro += coefEsq;
            else if( M.pontosU[i][j].tipoBoundary == INFLOW )
                acrescentaCentro += coefEsq; //Dirichlet na velocidade vira neumann na Psi
            else if( M.pontosU[i][j].tipoBoundary == KNOWN_PRESSURE ) {
                // Nao faz nada. Assumindo que Psi = 0 neste teste de canal ()
            }
            else if( M.pontosU[i][j].tipoBoundary == PERIODICIDADE ) {
                // Nao faz nada, pois ja fiz la embaixo
            }
            else if( M.pontosU[i][j].tipoBoundary == EIXO_AXISSIMETRIA )
                acrescentaCentro += coefEsq;
            else {
                DesenhaMalhaVTK(M, 0);
                DeuProblema("\n\n MontaMatrizPoisson: Problema Esquerda %d %d\n\n", i, j);
            }
        }

        // Colocando o coeficiente do Psi(i, j-1)
        if( j!=0 && (M.pontosP[i][j-1].tipo==FULL || M.pontosP[i][j-1].tipo==SURFACE) ) {
            valores[qtd] = coefBaixo;
            colunas[qtd++] = M.matIndicesP[i][j-1];
        }
        else {
            if( M.pontosV[i][j].tipoBoundary == NOSLIP )
                acrescentaCentro += coefBaixo; //Dirichlet na velocidade vira neumann na Psi
            else if( M.pontosV[i][j].tipoBoundary == SLIP )
                acrescentaCentro += coefBaixo; //Dirichlet na velocidade vira neumann na Psi
            else if( M.pontosV[i][j].tipoBoundary == INFLOW )
                acrescentaCentro += coefBaixo; //Dirichlet na velocidade vira neumann na Psi
            else if( M.pontosV[i][j].tipoBoundary == NEUMANN )
                acrescentaCentro += -coefBaixo;
            else if( M.pontosV[i][j].tipoBoundary == SIMETRIA )
                acrescentaCentro += coefBaixo;
            else {
                DeuProblema("\n\n MontaMatrizPoisson: Problema Baixo\n\n");
            }
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
        }

        //Colocando o coeficiente do Psi(i, j)
        centro = qtd;
        valores[qtd] += coefCentro;
        colunas[qtd++] = M.matIndicesP[i][j];

        //Colocando o coeficiente do Psi(i, j+1)
        if( j!=M.Ny-1 && (M.pontosP[i][j+1].tipo==FULL || M.pontosP[i][j+1].tipo==SURFACE) ) {
            valores[qtd] = coefCima;
            colunas[qtd++] = M.matIndicesP[i][j+1];
        }
        else {
            if( M.pontosV[i][j+1].tipoBoundary == NOSLIP )
                acrescentaCentro += coefCima; //dirichlet na velocidade vira neumann na Psi
            else if( M.pontosV[i][j+1].tipoBoundary == SLIP )
                acrescentaCentro += coefCima; //dirichlet na velocidade vira neumann na Psi
            else if( M.pontosV[i][j+1].tipoBoundary == INFLOW )
                acrescentaCentro += coefCima; //dirichlet na velocidade vira neumann na Psi
            else if( M.pontosV[i][j+1].tipoBoundary == NEUMANN )
                acrescentaCentro += -coefCima;
            else if( M.pontosV[i][j+1].tipoBoundary == SIMETRIA )
                acrescentaCentro += coefCima;
            else DeuProblema("\n MontaMatrizPoisson: Problema Cima poisson %d %d\n", i, j);
        }

        //Colocando o coeficiente do Psi(i+1, j)
        if( i!=M.Nx-1 && (M.pontosP[i+1][j].tipo==FULL || M.pontosP[i+1][j].tipo==SURFACE) ) {
            valores[qtd] = coefDir;
            colunas[qtd++] = M.matIndicesP[i+1][j];
        }
        else {
            if( M.pontosU[i+1][j].tipoBoundary == NEUMANN )
                acrescentaCentro += -coefDir; //Neumann na velocidade, vira dirichlet na Psi
            else if( M.pontosU[i+1][j].tipoBoundary == NOSLIP )
                acrescentaCentro += coefDir; //Dirichlet na velocidade, vira neumann na Psi
            else if( M.pontosU[i+1][j].tipoBoundary == SLIP )
                acrescentaCentro += coefDir; //Dirichlet na velocidade, vira neumann na Psi
            else if( M.pontosU[i+1][j].tipoBoundary == SIMETRIA )
                acrescentaCentro += coefDir;
            else if( M.pontosU[i+1][j].tipoBoundary == INFLOW )
                acrescentaCentro += coefDir; //Dirichlet na velocidade, vira neumann na Psi
            else if( M.pontosU[i+1][j].tipoBoundary == PERIODICIDADE ) {
                // Nao faz nada, pois ja fiz la em cima
            }
            else {
                DesenhaMalhaVTK(M, 0);
                DeuProblema("\n MontaMatrizPoisson: Problema direita poisson %d %d\n", i, j);
            }
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
        }

        /// === Finaliza verificando se ta na esquerda e eh condicao de PERIODICIDADE
        /// === Coloca um termo do lado direito do dominio, caso seja
        if( M.pontosU[i][j].tipoBoundary == PERIODICIDADE ) {
            valores[qtd] = coefEsq;
            colunas[qtd++] = M.matIndicesP[M.Nx-1][j];
        }

        valores[centro] += acrescentaCentro;
        MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES);
    }

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    return 0;
}

PetscErrorCode VetorIndependentePoisson(Vec VetVetor, double **U, double **V, double **PVelho, double **PNovo, double **Txx, double **Txy, double **Tyy, double **Psi, double **ViscosityNewtGen, MALHA M, LISTA_CELULAS_SURFACE *ListaSurf)
{
    int linha=0;
    int i, j;
    double valor;
    int iStart, iEnd;
    double rImaisMeio, rImenosMeio;
    PetscErrorCode ierr;

    ierr = VecGetOwnershipRange(VetVetor, &iStart, &iEnd); CHKERRQ(ierr);


    /// Se for uma simulacao com tensao superficial, primeiramente vou calcular a curvatura em cada celula surface
    if( M.tensaoSuperficial ) {
        DeterminaVetorNormalDeCadaCelula(&M, ListaSurf, 1.5*M.minDxDy);
        DeterminaCurvaturaDeCadaCelulaVersao2(&M, ListaSurf, 100);

//        DeterminaCurvaturaDeCadaCelulaVersaoBezier(&M, ListaSurf, 100);
    }

    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesP[linha].i;
        j = M.indicesP[linha].j;

        if( M.tipoCoord==AXI_CILINDRICO ) {
            rImaisMeio = M.r[i] + 0.5*M.dx[i];
            rImenosMeio = M.r[i] - 0.5*M.dx[i];
        }
        else
            rImaisMeio = rImenosMeio = 1.0;


        if( M.pontosP[i][j].tipo==EMPTY )
            DeuProblema("\n PROBLEMA \n");
        else if( M.pontosP[i][j].tipo==SURFACE ) {

            // Formulacao tradicional - tensao normal
            valor = LadoDireitoContornoPoisson(M, i, j, U, V, PVelho, Txx, Txy, Tyy, ViscosityNewtGen);
            ierr = VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES); CHKERRQ(ierr);

            // Formulacao laplaciano superf explicito - debora
            //valor = Psi[i][j];
            //ierr = VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES); CHKERRQ(ierr);

            // Formulacao laplaciano superf implicito - debora
            //LaplacianoSupImplicitoRHS(VetVetor, &M, i, j, linha, U, V, PVelho);
        }
        else {
            valor = (1.0/M.r[i])*( (rImaisMeio*U[i+1][j] - rImenosMeio*U[i][j])/M.dx[i] ) + ( (V[i][j+1]-V[i][j])/M.dy[j] );
            ierr = VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(VetVetor); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(VetVetor); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode LinhaCondicaoContornoPoisson(Mat A, MALHA *M, double **ViscosityNewtGen, int i, int j, int Linha)
{
    static int colunas[6], indiceAcrescenta1, indiceAcrescenta2, indiceAcrescenta3, indiceAcrescenta4, indiceAcrescentaCentro;
    static double valores[6], acrescenta1, acrescenta2, acrescenta3, acrescenta4, acrescentaCentro, nx, ny;
    static double h1, h2, coef1, coef2, valor;
    static double dx1, dx2;
    int qtd = 0;

    //double rImenosMeio, rImaisMeio, rI, 
    double rI, cilindrico;
    double coef3;
    cilindrico = (M->tipoCoord==AXI_CILINDRICO) ? 1.0 : 0.0;

    if( M->celulas[i][j]!=SURFACE ) {
        DeuProblema("1: PROBLEMA NO CONTORNO POISSON: %d %d\n", i, j);
    }

    //Colocando zero soh pra testar
//    if( LADO_HORIZONTAL_BAIXO(i, j) ) {
//        PrintDebug("BAM\n");
//        valores[0] = 1.0;
//        colunas[0] = M->matIndicesP[i][j];
//        qtd = 1;
//        MatSetValues(A, 1, &Linha, qtd, colunas, valores, INSERT_VALUES);
//        return 0;
//    }

    /// We adapt the "beta" that multiplies the diffusion term depending on the model
    double beta = ( M->tipo_modelo==MODELO_EVPT ) ? M->eta_inf : M->beta;

    //Valores da coordenada r usados no caso axisimetrico
    if( M->tipoCoord==AXI_CILINDRICO ) {
    	rI = M->r[i];
        //rImenosMeio = M.r[i] - 0.5*M.dx[i];
        //rImaisMeio = M.r[i] + 0.5*M.dx[i];
    }
    else
        rI = 1.0;

    /// Fazendo por casos da superficie livre
    if( LADO_HORIZONTAL_CIMA(i, j) ) {
        qtd = 0;
        nx = 0.0;
        ny = 1.0;
        coef1 = (j==0) ? (2.0*beta*M->dt*nx*ny)/(M->Re*0.5*(M->dy[j] + M->dy[j])) : (2.0*beta*M->dt*nx*ny)/(M->Re*0.5*(M->dy[j-1] + M->dy[j]));
        coef2 = (2.0*beta*M->dt*(nx*nx - ny*ny))/(M->Re);
	coef3 = cilindrico*(2.0*beta*M->dt*((ny*ny)/rI))/(M->Re);
        h1 = (i==0) ? M->dx[i] : 0.5*(M->dx[i-1] + M->dx[i]);
        h2 = (i==M->Nx-1) ? M->dx[i] : 0.5*(M->dx[i] + M->dx[i+1]);

        indiceAcrescenta1 = indiceAcrescenta2 = indiceAcrescenta3 = indiceAcrescenta4 = -1;
        acrescenta1 = acrescenta2 = acrescenta3 = acrescenta4 = 0.0;

        //Phi(i-1, j-1)
        valor = - 2.0*coef1*h2/(h1*(h1+h2));
        if( (i!=0) && (j!=0) && M->celulas[i-1][j-1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j-1];
        }
//        else if( j!=0 && M->pontosU[i][j-1].tipo==BOUNDARY ){
//            if( M->pontosU[i][j-1].tipoBoundary==NOSLIP )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i][j-1].tipoBoundary==INFLOW )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i][j-1].tipoBoundary==NEUMANN )
//                acrescenta1 -= valor; //neumann na velocidade vira dirichlet na Psi
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                ImprimeInterfaceVTK(*M, 10000);
//                DeuProblema("2a: Problema Contorno Poisson %d %d\n", i, j);
//            }
//        }
//        else if( i!=0 && M->pontosV[i-1][j].tipo==BOUNDARY ) { //Joga pra cima
//            if( M->pontosV[i-1][j].tipoBoundary==NOSLIP )
//                acrescenta3 += valor;
//            else if( M->pontosV[i-1][j].tipoBoundary==INFLOW )
//                acrescenta3 += valor;
//            else if( M->pontosV[i-1][j].tipoBoundary==NEUMANN )
//                acrescenta3 -= valor;
//            else DeuProblema("2b: Problema Contorno Poisson %d %d\n", i, j);
//        }

        //Phi(i-1, j)
        valor = ( 2.0*coef1*h2/(h1*(h1+h2)) ) - ( coef2*2.0/(h1*(h1+h2)) ) - ( coef3*h2/(h1*(h1+h2)) );
        if( (i!=0) && M->celulas[i-1][j]!=EMPTY ) {
            indiceAcrescenta3 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j];
        }
//        else {
//            if( M->pontosU[i][j].tipoBoundary==NOSLIP )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i][j].tipoBoundary==INFLOW )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i][j].tipoBoundary==NEUMANN )
//                acrescenta2 -= valor; //neumann na velocidade vira dirichlet na Psi
//            else DeuProblema("3: Problema Contorno Poisson %d %d\n", i, j);
//        }

        //Phi(i, j-1)
        valor = ( 2.0*coef1*(h2-h1)/(h1*h2) );
        if( (j!=0) && M->celulas[i][j-1]!=EMPTY ) {
            indiceAcrescenta1 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j-1];
        }
//        else if( i!=0 && M->pontosV[i][j].tipo==BOUNDARY ) { //Joga pra cima
//            if( M->pontosV[i][j].tipoBoundary==NOSLIP )
//                acrescenta2 += valor;
//            else if( M->pontosV[i][j].tipoBoundary==INFLOW )
//                acrescenta2 += valor;
//            else if( M->pontosV[i][j].tipoBoundary==NEUMANN )
//                acrescenta2 -= valor;
//            else DeuProblema("2b: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else DeuProblema("4: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j)
        valor = ( -2.0*coef1*(h2-h1)/(h1*h2) ) + ( coef2*2.0/(h1*h2) ) + ( coef3*(h2-h1)/(h1*h2) ) - 1.0;
        if( M->celulas[i][j]!=EMPTY ) {
            indiceAcrescenta2 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j];
        }
//        else DeuProblema("5: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i+1, j-1)
        valor = 2.0*coef1*h1/(h2*(h1+h2));
        if( (j!=0) && (i!=M->Nx-1) && M->celulas[i+1][j-1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j-1];
        }
//        else if( CELULA_FLUIDO(i, j-1) ) { //Joga pra esquerda
//            if( M->pontosU[i+1][j-1].tipoBoundary==NEUMANN )
//                acrescenta1 += -valor; //Vira dirichlet na psi
//            else if( M->pontosU[i+1][j-1].tipoBoundary==NOSLIP )
//                acrescenta1 += valor; //Vira neumann na psi
//            else if( M->pontosU[i+1][j-1].tipoBoundary==INFLOW )
//                acrescenta1 += valor; //Vira neumann na psi
//            else if( M->celulas[i][j-1]==SURFACE ) { } //Acontece muito de vez em quando. vou deixar ZERO
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                DeuProblema("\n 1. Problema CondicaoDeContornoPoisson: tipo de contorno %d %d\n\n", i, j);
//            }
//        }
//        else if( CELULA_FLUIDO(i+1, j) ) { //Joga pra cima
//            if( M->pontosV[i+1][j].tipoBoundary==NEUMANN )
//                acrescenta4 += -valor; //Vira dirichlet na psi
//            else if( M->pontosV[i+1][j].tipoBoundary==NOSLIP )
//                acrescenta4 += valor; //Vira neumann na psi
//            else if( M->pontosV[i+1][j].tipoBoundary==INFLOW )
//                acrescenta4 += valor; //Vira neumann na psi
//            else DeuProblema("\n\n Problema Contorno poisson: (i+1, j-1) horizontal cima %d %d\n", i, j);
//        }
//        else DeuProblema("\n\n Problema Contorno poisson: (i+1, j-1) horizontal cima %d %d\n", i, j);

        //Phi(i+1, j)
        valor = ( - 2.0*coef1*h1/(h2*(h1+h2)) ) - ( coef2*2.0/(h2*(h1+h2)) ) + ( coef3*h1/(h2*(h1+h2)) );
        if( (i!=M->Nx-1) && M->celulas[i+1][j]!=EMPTY ) {
            indiceAcrescenta4 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j];
        }
//        else {
//            if( M->pontosU[i+1][j].tipoBoundary==NEUMANN )
//                acrescenta2 += -valor; //Vira dirichlet na psi
//            else if( M->pontosU[i+1][j].tipoBoundary==INFLOW )
//                acrescenta2 += valor; //Vira dirichlet na psi
//            else if( M->pontosU[i+1][j].tipoBoundary==NOSLIP )
//                acrescenta2 += valor; //Vira dirichlet na psi
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                DeuProblema("\n 2. Problema CondicaoDeContornoPoisson: tipo de contorno %d %d\n\n", i, j);
//            }
//        }

        valores[indiceAcrescenta1] += acrescenta1;
        valores[indiceAcrescenta2] += acrescenta2;
        valores[indiceAcrescenta3] += acrescenta3;
        valores[indiceAcrescenta4] += acrescenta4;
    }
    else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
        qtd = 0;
        acrescenta1 = acrescenta2 = acrescenta3 = acrescentaCentro = 0.0;
        indiceAcrescenta1 = indiceAcrescenta2 = indiceAcrescenta3 = indiceAcrescentaCentro = -1;
        nx = 0.0;
        ny = -1.0;
        coef1 = (j==M->Ny-1) ? (2.0*beta*M->dt*nx*ny)/(M->Re*M->dy[j]) : (2.0*beta*M->dt*nx*ny)/(M->Re*0.5*(M->dy[j] + M->dy[j+1]));
        coef2 = (2.0*beta*M->dt*(nx*nx - ny*ny))/(M->Re);
	coef3 = cilindrico*(2.0*beta*M->dt*((ny*ny)/rI))/(M->Re);
        h1 = (i==0) ? M->dx[i] : 0.5*(M->dx[i-1] + M->dx[i]);
        h2 = (i==M->Nx-1) ? M->dx[i] : 0.5*(M->dx[i] + M->dx[i+1]);


        //Phi(i-1, j)
        valor = ( - 2.0*coef1*h2/(h1*(h1+h2)) ) - ( coef2*2.0/(h1*(h1+h2)) ) - ( coef3*h2/(h1*(h1+h2)) );
        if( (i!=0) && M->celulas[i-1][j]!=EMPTY ) {
            indiceAcrescenta1 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j];
        }
//        else {
//            if( M->pontosU[i][j].tipoBoundary==NOSLIP )
//                acrescentaCentro += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i][j].tipoBoundary==INFLOW )
//                acrescentaCentro += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i][j].tipoBoundary==NEUMANN )
//                acrescentaCentro -= valor; //neumann na velocidade vira dirichlet na Psi
//            else DeuProblema("7: Problema Contorno Poisson %d %d\n", i, j);
//        }


        //Phi(i-1, j+1)
        valor = ( 2.0*coef1*h2/(h1*(h1+h2)) );
        if( (j!=M->Ny-1) && (i!=0) && M->celulas[i-1][j+1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j+1];
        }
//        else if( (j!=M->Ny-1) && CELULA_FLUIDO(i, j+1) ) { //Joga pra direita
//            if( M->pontosU[i][j+1].tipoBoundary==NOSLIP )
//                acrescenta3 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i][j+1].tipoBoundary==INFLOW )
//                acrescenta3 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i][j+1].tipoBoundary==NEUMANN )
//                acrescenta3 -= valor; //neumann na velocidade vira dirichlet na Psi
//            else if( M->celulas[i][j+1]==SURFACE ) { } //Acontece muito de vez emq uando. vou deixar ZERO
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                DeuProblema("1. i-1 j+1: Problema Contorno Poisson %d %d\n", i, j);
//            }
//        }
//        else if ( CELULA_FLUIDO(i-1, j) ) {
//            if( M->pontosV[i-1][j+1].tipoBoundary == NOSLIP )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i-1][j+1].tipoBoundary == INFLOW )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i-1][j+1].tipoBoundary == NEUMANN )
//                acrescenta1 += -valor;
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                DeuProblema("2. i-1 j+1 problema contorno poisson: %d %d\n\n", i, j);
//            }
//        }

        //Phi(i, j)
        valor = ( 2.0*coef1*(h2-h1)/(h1*h2) ) + ( coef2*2.0/(h1*h2) ) + ( coef3*(h2-h1)/(h1*h2) ) - 1.0;
        if( M->celulas[i][j]!=EMPTY ) {
            indiceAcrescentaCentro = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j];
        }
//        else DeuProblema("9: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j+1)
        valor = ( -2.0*coef1*(h2-h1)/(h1*h2) );
        if( (j!=M->Ny-1) && M->celulas[i][j+1]!=EMPTY ) {
            indiceAcrescenta3 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j+1];
        }
//        else {
//            if( M->pontosV[i][j+1].tipoBoundary == INFLOW )
//                acrescentaCentro += valor; //dirichlet na velocidade vira neumann na Psi
//            else {
//                //ImprimeInterfaceVTK(*M, 10000);
//                //DesenhaMalhaVTK(*M, 0);
//                DeuProblema("10: Problema Contorno Poisson %d %d\n", i, j);
//            }
//        }
         //Phi(i+1, j)
        valor = ( 2.0*coef1*h1/(h2*(h1+h2)) ) - ( coef2*2.0/(h2*(h1+h2)) ) + ( coef3*h1/(h2*(h1+h2)) );
        if( (i!=M->Nx-1) && M->celulas[i+1][j]!=EMPTY ) {
            indiceAcrescenta2 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j];
        }
//        else {
//            if( M->pontosU[i+1][j].tipoBoundary==NEUMANN )
//                acrescentaCentro += -valor; //Vira dirichlet na psi
//            else if( M->pontosU[i+1][j].tipoBoundary==NOSLIP )
//                acrescentaCentro += valor; //Vira dirichlet na psi
//            else if( M->pontosU[i+1][j].tipoBoundary==INFLOW )
//                acrescentaCentro += valor; //Vira dirichlet na psi
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                DeuProblema("\n 3. Problema CondicaoDeContornoPoisson: tipo de contorno %d %d\n\n", i, j);
//            }
//        }

        //Phi(i+1, j+1)
        valor = - 2.0*coef1*h1/(h2*(h1+h2));
        if( (j!=M->Ny-1) && (i!=M->Nx-1) && M->celulas[i+1][j+1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j+1];
        }
//        else if( (i!=M->Nx-1) && CELULA_FLUIDO(i+1, j) && M->pontosV[i+1][j+1].tipo==BOUNDARY ) { //Joga pra baixo
//            if( M->pontosV[i+1][j+1].tipoBoundary == NOSLIP )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i+1][j+1].tipoBoundary == INFLOW )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i+1][j+1].tipoBoundary == NEUMANN )
//                acrescenta2 += -valor;
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                ImprimeInterfaceVTK(*M, 100000);
//                DeuProblema("12: Problema Contorno Poisson %d %d\n", i, j);
//            }
//        }
//        else if( CELULA_FLUIDO(i, j+1) && M->pontosU[i+1][j+1].tipo==BOUNDARY ) { //Joga pra esquerda
//            if( M->pontosU[i+1][j+1].tipoBoundary == NOSLIP )
//                acrescenta3 += valor;
//            else if( M->pontosU[i+1][j+1].tipoBoundary == INFLOW )
//                acrescenta3 += valor;
//            else if( M->pontosU[i+1][j+1].tipoBoundary == NEUMANN )
//                acrescenta3 += -valor;
//            else DeuProblema("\n\n13: Problema no contorno poisson %d %d\n\n", i, j);
//        }
        //else DeuProblema("\n\n Problema contorno poisson: (i+1, j+1) horizontal baixo %d %d\n", i, j);

        valores[indiceAcrescenta1] += acrescenta1;
        valores[indiceAcrescenta2] += acrescenta2;
        valores[indiceAcrescenta3] += acrescenta3;
        valores[indiceAcrescentaCentro] += acrescentaCentro;
    }
    else if( LADO_VERTICAL_DIREITA(i, j) ) {
        qtd = 0;
        nx = 1.0;
        ny = 0.0;
        coef1 = (i==0) ? (2.0*beta*M->dt*nx*ny)/(M->Re*(M->dx[i])) : (2.0*beta*M->dt*nx*ny)/(M->Re*0.5*(M->dx[i-1] + M->dx[i]));
        coef2 = (2.0*beta*M->dt*(ny*ny - nx*nx))/(M->Re);
	coef3 = cilindrico*(2.0*beta*M->dt*((nx*nx)/rI))/(M->Re);
        h1 = (j==0) ? M->dy[j] : 0.5*(M->dy[j-1] + M->dy[j]);
        h2 = (j==M->Ny-1) ? M->dy[j] : 0.5*(M->dy[j] + M->dy[j+1]);

        indiceAcrescenta1 = indiceAcrescenta2 = indiceAcrescenta3 = indiceAcrescenta4 = -1;
        acrescenta1 = acrescenta2 = acrescenta3 = acrescenta4 = 0.0;

        //Phi(i-1, j-1)
        valor = - 2.0*coef1*h2/(h1*(h1+h2));
        if( (i!=0) && (j!=0) && CELULA_FLUIDO(i-1, j-1) ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j-1];
        }
//        else if( CELULA_FLUIDO(i, j-1) ) { //Joga pra direita
//            if( M->pontosU[i][j-1].tipoBoundary==INFLOW )
//                acrescenta1 += valor;
//            else if( M->pontosU[i][j-1].tipoBoundary==NOSLIP )
//                acrescenta1 += valor;
//            else if( M->pontosU[i][j-1].tipoBoundary==NEUMANN )
//                acrescenta1 += - valor;
//            else if( M->celulas[i][j-1]==SURFACE ) { } //Acontece de veeez em quando. vou apenas colocar zero aqui
//            else DeuProblema("1dir: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else if( (i!=0) && CELULA_FLUIDO(i-1, j) ) { //Joga pra cima
//            if( M->pontosV[i-1][j].tipoBoundary==INFLOW )
//                acrescenta4 += valor;
//            else if( M->pontosV[i-1][j].tipoBoundary==NOSLIP )
//                acrescenta4 += valor;
//            else if( M->pontosV[i-1][j].tipoBoundary==NEUMANN )
//                acrescenta4 += - valor;
//            else DeuProblema("1dir: Problema Contorno Poisson %d %d\n", i, j);
//        }
        //else DeuProblema("1dir: Problema Contorno Poisson %d %d\n", i, j);



        //Phi(i-1, j)
        valor = 2.0*coef1*(h2-h1)/(h1*h2) - ( coef3/h1 );
        if( (i!=0) && CELULA_FLUIDO(i-1, j) ) {
            indiceAcrescenta4 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j];
        }
//        else if( CELULA_FLUIDO(i, j) ) { //Joga pra direita
//            if( M->pontosU[i][j].tipoBoundary==INFLOW )
//                acrescenta2 += valor;
//            else if( M->pontosU[i][j].tipoBoundary==NOSLIP )
//                acrescenta2 += valor;
//            else if( M->pontosU[i][j].tipoBoundary==NEUMANN )
//                acrescenta2 += - valor;
//            else DeuProblema("2.1dir: Problema Contorno Poisson %d %d\n", i, j);
//        }



        //Phi(i-1, j+1)
        valor = 2.0*coef1*h1/(h2*(h1+h2));
        if( (i!=0) && (j!=M->Ny-1) && CELULA_FLUIDO(i-1, j+1) ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j+1];
        }
//        else if( (j!=M->Ny-1) && CELULA_FLUIDO(i, j+1) ) { //Joga pra direita
//            if( M->pontosU[i][j+1].tipoBoundary==INFLOW )
//                acrescenta3 += valor;
//            else if( M->pontosU[i][j+1].tipoBoundary==NOSLIP )
//                acrescenta3 += valor;
//            else if( M->pontosU[i][j+1].tipoBoundary==NEUMANN )
//                acrescenta3 += - valor;
//            else if( M->celulas[i][j+1]==SURFACE ) { } //AContece muito de vez em quando, vou deixar zero aqui
//            else DeuProblema("3.1dir: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else if( (i!=0) && CELULA_FLUIDO(i-1, j) ) { //Joga pra baixo
//            if( M->pontosV[i-1][j+1].tipoBoundary==INFLOW )
//                acrescenta4 += valor;
//            else if( M->pontosV[i-1][j+1].tipoBoundary==NOSLIP )
//                acrescenta4 += valor;
//            else if( M->pontosV[i-1][j+1].tipoBoundary==NEUMANN )
//                acrescenta4 += - valor;
//            else DeuProblema("1dir: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else {
//            DesenhaMalhaVTK(*M, 0);
//            ImprimeInterfaceVTK(*M, 10000);
//            DeuProblema("2.2dir: Problema Contorno Poisson %d %d\n", i, j);
//        }



         //Phi(i, j-1)
        valor = ( 2.0*coef1*h2/(h1*(h1+h2)) ) - ( coef2*2.0/(h1*(h1+h2)) );
        if( (j!=0) && CELULA_FLUIDO(i, j-1) ) {
            indiceAcrescenta1 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j-1];
        }
//        else if( CELULA_FLUIDO(i, j) ){
//            if( M->pontosV[i][j].tipoBoundary==INFLOW )
//                acrescenta2 += valor;
//            else if( M->pontosV[i][j].tipoBoundary==NOSLIP )
//                acrescenta2 += valor;
//            else if( M->pontosV[i][j].tipoBoundary==NEUMANN )
//                acrescenta2 += - valor;
//            else DeuProblema("3.2dir: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else DeuProblema("4dir: Problema Contorno Poisson %d %d\n", i, j);



        //Phi(i, j)
        valor = ( - 2.0*coef1*(h2-h1)/(h1*h2) ) + ( coef2*2.0/(h1*h2) ) + ( coef3/h1 ) - 1.0;
        if( CELULA_FLUIDO(i, j) ) {
            indiceAcrescenta2 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j];
        }
//        else DeuProblema("5dir: Problema Contorno Poisson %d %d\n", i, j);



        //Phi(i, j+1)
        valor = ( - 2.0*coef1*h1/(h2*(h1+h2)) ) - ( coef2*2.0/(h2*(h1+h2)) );
        if( (j!=M->Ny-1) && CELULA_FLUIDO(i, j+1) ) {
            indiceAcrescenta3 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j+1];
        }
//        else if( CELULA_FLUIDO(i, j) ){
//            if( M->pontosV[i][j+1].tipoBoundary==INFLOW )
//                acrescenta2 += valor;
//            else if( M->pontosV[i][j+1].tipoBoundary==NOSLIP )
//                acrescenta2 += valor;
//            else if( M->pontosV[i][j+1].tipoBoundary==NEUMANN )
//                acrescenta2 += - valor;
//            else DeuProblema("3.3dir: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else DeuProblema("5dir: Problema Contorno Poisson %d %d\n", i, j);

        valores[indiceAcrescenta1] += acrescenta1;
        valores[indiceAcrescenta2] += acrescenta2;
        valores[indiceAcrescenta3] += acrescenta3;
        valores[indiceAcrescenta4] += acrescenta4;
    }
    else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
        qtd = 0;
        nx = -1.0;
        ny = 0.0;
        coef1 = (i==M->Nx-1) ? (2.0*beta*M->dt*nx*ny)/(M->Re*M->dx[i]) : (2.0*beta*M->dt*nx*ny)/(M->Re*0.5*(M->dx[i] + M->dx[i+1]));
        coef2 = (2.0*beta*M->dt*(ny*ny - nx*nx))/(M->Re);
	coef3 = cilindrico*(2.0*beta*M->dt*((nx*nx)/rI))/(M->Re);
        h1 = (j==0) ? M->dy[j] : 0.5*(M->dy[j-1] + M->dy[j]);
        h2 = (j==M->Ny-1) ? M->dy[j] : 0.5*(M->dy[j] + M->dy[j+1]);

        indiceAcrescenta1 = indiceAcrescenta2 = indiceAcrescenta3 = indiceAcrescenta4 = -1;
        acrescenta1 = acrescenta2 = acrescenta3 = acrescenta4 = 0.0;

        //Phi(i, j-1)
        valor = ( - 2.0*coef1*h2/(h1*(h1+h2)) ) - ( coef2*2.0/(h1*(h1+h2)) );
        if( (j!=0) && M->celulas[i][j-1]!=EMPTY ) {
            indiceAcrescenta1 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j-1];
        }
//        else if( CELULA_FLUIDO(i, j) ) { //Joga pra cima
//            if( M->pontosV[i][j].tipoBoundary==INFLOW )
//                acrescenta4 += valor;
//            else if( M->pontosV[i][j].tipoBoundary==NOSLIP )
//                acrescenta4 += valor;
//            else if( M->pontosV[i][j].tipoBoundary==NEUMANN )
//                acrescenta4 += - valor;
//            else DeuProblema("1esq: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else DeuProblema("1esq: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j)
        valor = ( 2.0*coef1*(h2-h1)/(h1*h2) ) + ( coef2*2.0/(h1*h2) ) - ( coef3/h2 ) - 1.0;
        if( M->celulas[i][j]!=EMPTY ) {
            indiceAcrescenta4 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j];
        }
//        else DeuProblema("2esq: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j+1)
        valor = ( 2.0*coef1*h1/(h2*(h1+h2)) ) - ( coef2*2.0/(h2*(h1+h2)) );
        if( (j!=M->Ny-1) && CELULA_FLUIDO(i, j+1) ) {
            indiceAcrescenta2 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j+1];
        }
//        else if( CELULA_FLUIDO(i, j) ) { //Joga pra baixo
//            if( M->pontosV[i][j+1].tipoBoundary==INFLOW )
//                acrescenta4 += valor;
//            else if( M->pontosV[i][j+1].tipoBoundary==NOSLIP )
//                acrescenta4 += valor;
//            else if( M->pontosV[i][j+1].tipoBoundary==NEUMANN )
//                acrescenta4 += - valor;
//            else DeuProblema("3esq: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else DeuProblema("3esq: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i+1, j-1)
        valor = 2.0*coef1*h2/(h1*(h1+h2));
        if( (i!=M->Nx-1) && (j!=0) && M->celulas[i+1][j-1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j-1];
        }
//        else if( CELULA_FLUIDO(i, j-1) ) { //Joga pra esquerda
//            if( M->pontosU[i+1][j-1].tipoBoundary==INFLOW )
//                acrescenta1 += valor;
//            else if( M->pontosU[i+1][j-1].tipoBoundary==NOSLIP )
//                acrescenta1 += valor;
//            else if( M->pontosU[i+1][j-1].tipoBoundary==NEUMANN )
//                acrescenta1 += -valor;
//            else if( M->celulas[i][j-1] ) { } //Acontece mto de vez em qdo. vou deixar zero
//            else DeuProblema("4esq: Problema Contorno Poisson %d %d\n", i, j);
//        }
        //else DeuProblema("4esq: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i+1, j)
        valor = - 2.0*coef1*(h2-h1)/(h1*h2) + ( coef3/h2 );
        if( (i!=M->Nx-1) && M->celulas[i+1][j]!=EMPTY ) {
            indiceAcrescenta3 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j];
        }
//        else if( CELULA_FLUIDO(i, j) ) { //Joga pra esquerda
//            if( M->pontosU[i+1][j].tipoBoundary==INFLOW )
//                acrescenta4 += valor;
//            else if( M->pontosU[i+1][j].tipoBoundary==NOSLIP )
//                acrescenta4 += valor;
//            else if( M->pontosU[i+1][j].tipoBoundary==NEUMANN )
//                acrescenta4 += -valor;
//            else DeuProblema("4.5esq: Problema Contorno Poisson %d %d\n", i, j);
//        }

        //Phi(i+1, j+1)
        valor = - 2.0*coef1*h1/(h2*(h1+h2));
        if( (i!=M->Nx-1) && (j!=M->Ny-1) && CELULA_FLUIDO(i+1, j+1) ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j+1];
        }
//        else if( (j!=M->Ny-1) && CELULA_FLUIDO(i, j+1) ) { //Joga pra esquerda
//            if( M->pontosU[i+1][j+1].tipoBoundary==INFLOW )
//                acrescenta2 += valor;
//            else if( M->pontosU[i+1][j+1].tipoBoundary==NOSLIP )
//                acrescenta2 += valor;
//            else if( M->pontosU[i+1][j+1].tipoBoundary==NEUMANN )
//                acrescenta2 += -valor;
//            else if( M->celulas[i][j+1]==SURFACE ) { } //Acontece de veeeez em quando. vou simplesmente botar zero aqui
//            else DeuProblema("5esq: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else if( (i!=M->Nx-1) && CELULA_FLUIDO(i+1, j) ) { //Joga pra baixo
//            if( M->pontosV[i+1][j+1].tipoBoundary==INFLOW )
//                acrescenta3 += valor;
//            else if( M->pontosV[i+1][j+1].tipoBoundary==NOSLIP )
//                acrescenta3 += valor;
//            else if( M->pontosV[i+1][j+1].tipoBoundary==NEUMANN )
//                acrescenta3 += - valor;
//            else DeuProblema("5esq: Problema Contorno Poisson %d %d\n", i, j);
//        }
        //dont bother in this case. phi=0   else DeuProblema("5esq: Problema Contorno Poisson %d %d\n", i, j);

        valores[indiceAcrescenta1] += acrescenta1;
        valores[indiceAcrescenta2] += acrescenta2;
        valores[indiceAcrescenta3] += acrescenta3;
        valores[indiceAcrescenta4] += acrescenta4;
    }
    else if( QUINA_CIMA_DIREITA(i, j) ) {
        qtd = 0;
        nx = sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;
        h1 = (j==0) ? 2.0*M->dy[j] : (M->dy[j-1] + M->dy[j]);
        coef1 = (i==0) ? (2.0*beta*M->dt*nx*ny)/(M->Re*0.50*(M->dx[i])*h1) : (2.0*beta*M->dt*nx*ny)/(M->Re*0.25*(M->dx[i-1]+M->dx[i])*h1);
        coef2 = (2.0*beta*M->dt*(nx*nx - ny*ny))/(M->Re);
	coef3 = cilindrico*(2.0*beta*M->dt*((ny*ny)/rI))/(M->Re);
        h1 = (i==0) ? M->dx[i] : 0.5*(M->dx[i-1] + M->dx[i]);//Observe que ele redefine esse valor!
        h2 = (i==M->Nx-1) ? M->dx[i] : 0.5*(M->dx[i] + M->dx[i+1]);

        indiceAcrescenta1 = indiceAcrescenta2 = indiceAcrescenta3 = -1;
        acrescenta1 = acrescenta2 = acrescenta3 = 0.0;

        //Phi(i-1, j-1)
        valor = - 2.0*coef1;
        if( i!=0 && j!=0 && M->celulas[i-1][j-1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j-1];
        }
//        else if( i==0 && j!=0 ) {
//            if( M->pontosU[i][j-1].tipoBoundary==INFLOW )
//                acrescenta1 += valor;
//            else if( M->pontosU[i][j-1].tipoBoundary==NOSLIP )
//                acrescenta1 += valor;
//            else if( M->pontosU[i][j-1].tipoBoundary==NEUMANN )
//                acrescenta1 += - valor;
//            else DeuProblema("1acimadir: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else {
//            DesenhaMalhaVTK(*M, 0);
//            DeuProblema("1bcimadir: Problema Contorno Poisson %d %d\n", i, j);
//        }

        //Phi(i-1, j)
        valor = ( 2.0*coef1 ) - ( coef2*2.0/(h1*(h1+h2)) ) - ( coef3/h1 );
        if( i!=0 && M->celulas[i-1][j]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j];
        }
//        else if( i==0 ) {
//            if( M->pontosU[i][j].tipoBoundary==INFLOW )
//                acrescenta2 += valor;
//            else if( M->pontosU[i][j].tipoBoundary==NOSLIP )
//                acrescenta2 += valor;
//            else if( M->pontosU[i][j].tipoBoundary==NEUMANN )
//                acrescenta2 += - valor;
//            else DeuProblema("2cimadir: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else DeuProblema("2cimadir: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j-1)
        valor = 2.0*coef1;
        if( j!=0 && M->celulas[i][j-1]!=EMPTY ) {
            indiceAcrescenta1 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j-1];
        }
        //else DeuProblema("3cimadir: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j)
        valor = ( - 2.0*coef1 ) + ( coef2*2.0/(h1*h2) ) + ( coef3/h1 ) - 1.0;
        if( M->celulas[i][j]!=EMPTY ) {
            indiceAcrescenta2 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j];
        }
//        else DeuProblema("4cimadir: Problema Contorno Poisson %d %d\n", i, j);

        valores[indiceAcrescenta1] += acrescenta1;
        valores[indiceAcrescenta2] += acrescenta2;
    }
    else if( QUINA_CIMA_ESQUERDA(i, j) ) {
        // if( j==0 ) {
        //     DesenhaMalhaVTK(*M, 0);
        //     ImprimeInterfaceVTK(*M, 100000);
        //     DeuProblema("HUFAHUFA\n");
        // }
            

        qtd = 0;
        nx = -sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;
        dx1 = M->dx[i];
        dx2 = (i==M->Nx-1) ? M->dx[i] : M->dx[i+1];

        h1 = (j==0) ? 2.0*M->dy[j] : (M->dy[j-1] + M->dy[j]);
        coef1 = (2.0*beta*M->dt*nx*ny)/(M->Re*0.25*(dx1+dx2)*h1);
        coef2 = (2.0*beta*M->dt*(nx*nx - ny*ny))/(M->Re);
	coef3 = cilindrico*(2.0*beta*M->dt*((ny*ny)/rI))/(M->Re);
        h1 = (i==0) ? M->dx[i] : 0.5*(M->dx[i-1] + M->dx[i]);
        h2 = (i==M->Nx-1) ? M->dx[i] : 0.5*(M->dx[i] + M->dx[i+1]);

        acrescenta1 = acrescenta2 = 0.0;
        indiceAcrescenta1 = indiceAcrescenta2 = -1;

        //Phi(i, j-1)
        valor = - 2.0*coef1;
        if( (j!=0) && M->celulas[i][j-1]!=EMPTY ) {
            indiceAcrescenta1 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j-1];
        }
        //else DeuProblema("1esqCima: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j)
        valor = ( 2.0*coef1 ) + ( coef2*2.0/(h1*h2) ) - ( coef3/h2 ) - 1.0;
        if( M->celulas[i][j]!=EMPTY ) {
            indiceAcrescenta2 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j];
        }
//        else DeuProblema("2esqCima Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i+1, j-1)
        valor = 2.0*coef1;
        if( (j!=0) && (i!=M->Nx-1) && M->celulas[i+1][j-1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j-1];
        }
//        else if( CELULA_FLUIDO(i, j-1) && M->pontosU[i+1][j-1].tipo==BOUNDARY ) {
//            if( M->pontosU[i+1][j-1].tipoBoundary == NOSLIP )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i+1][j-1].tipoBoundary == INFLOW )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i+1][j-1].tipoBoundary == NEUMANN )
//                acrescenta1 += - valor;
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                ImprimeInterfaceVTK(*M, 100000);
//                DeuProblema("5.1esqcima: Problema Contorno Poisson %d %d\n", i, j);
//            }
//        }
        //else DeuProblema("5esqCima: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i+1, j)
        valor = ( - 2.0*coef1 ) - ( coef2*2.0/(h2*(h1+h2)) ) + ( coef3/h2 );
        if( (i!=M->Nx-1) && M->celulas[i+1][j]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j];
        }
//        else if( CELULA_FLUIDO(i, j) ) {
//            if( M->pontosU[i+1][j].tipoBoundary == NOSLIP )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i+1][j].tipoBoundary == INFLOW )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i+1][j].tipoBoundary == NEUMANN )
//                acrescenta2 += - valor;
//            else DeuProblema("4.1esqcima: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else DeuProblema("4esqCima: Problema Contorno Poisson %d %d\n", i, j);

        valores[indiceAcrescenta1] += acrescenta1;
        valores[indiceAcrescenta2] += acrescenta2;
    }
    else if( QUINA_BAIXO_DIREITA(i, j) ) {
        qtd = 0;
        acrescenta1 = acrescenta2 = acrescenta3 = 0.0;
        indiceAcrescenta1 = indiceAcrescenta2 = indiceAcrescenta3 = -1;
        nx = sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;

        double denominador = M->Re;
        denominador = (j==M->Ny-1) ? denominador*M->dy[j] : denominador*0.50*(M->dy[j]+M->dy[j+1]);
        denominador = (i==0) ? denominador*M->dx[i] : denominador*0.50*(M->dx[i-1]+M->dx[i]);
        coef1 = (2.0*beta*M->dt*nx*ny)/denominador;
        coef2 = (2.0*beta*M->dt*(nx*nx - ny*ny))/(M->Re);
	coef3 = cilindrico*(2.0*beta*M->dt*((ny*ny)/rI))/(M->Re);
        h1 = (i==0) ? M->dx[i] : 0.5*(M->dx[i-1] + M->dx[i]);
        h2 = (i==M->Nx-1) ? M->dx[i] : 0.5*(M->dx[i] + M->dx[i+1]);

        //Phi(i-1, j)
        valor = ( - 2.0*coef1 ) - ( coef2*2.0/(h1*(h1+h2)) ) - ( coef3/h1 );
        if( (i!=0) && M->celulas[i-1][j]!=EMPTY ) {
            indiceAcrescenta1 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j];
        }
        //else DeuProblema("1baixodir: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i-1, j+1)
        valor = 2.0*coef1;
        if( (i!=0) && (j!=M->Ny-1) && M->celulas[i-1][j+1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i-1][j+1];
        }
//        else if( (i==0) && (j!=M->Ny-1) ) {
//            if( M->pontosU[i][j+1].tipoBoundary==INFLOW )
//                acrescenta3 += valor;
//            else if( M->pontosU[i][j+1].tipoBoundary==NOSLIP )
//                acrescenta3 += valor;
//            else if( M->pontosU[i][j+1].tipoBoundary==NEUMANN )
//                acrescenta3 -= valor;
//            else DeuProblema("2.1baixodir: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else { //Joga pra baixo
//            if( M->pontosV[i-1][j+1].tipoBoundary == NOSLIP )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i-1][j+1].tipoBoundary == INFLOW )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i-1][j+1].tipoBoundary == NEUMANN )
//                acrescenta1 += - valor;
//            else if( M->celulas[i-1][j]==SURFACE ) { } //Acontece muito de vez em quando. vou deixar ZERO
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                DeuProblema("2.2baixodir: Problema Contorno Poisson %d %d\n", i, j);
//            }
//        }

        //Phi(i, j)
        valor = ( 2.0*coef1 ) + ( coef2*2.0/(h1*h2) ) + ( coef3/h1 ) - 1.0;
        if( M->celulas[i][j]!=EMPTY ) {
            indiceAcrescenta2 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j];
        }
//        else DeuProblema("3baixodir: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j+1)
        valor = - 2.0*coef1;
        if( (j!=M->Ny-1) && M->celulas[i][j+1]!=EMPTY ) {
            indiceAcrescenta3 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j+1];
        }
//        else {
//            if( M->pontosV[i][j+1].tipoBoundary == NOSLIP )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i][j+1].tipoBoundary == INFLOW )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i][j+1].tipoBoundary == NEUMANN )
//                acrescenta2 += - valor;
//            else DeuProblema("4baixodir: Problema Contorno Poisson %d %d\n", i, j);
//        }

        valores[indiceAcrescenta1] += acrescenta1;
        valores[indiceAcrescenta2] += acrescenta2;
        valores[indiceAcrescenta3] += acrescenta3;
    }
    else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
        qtd = 0;
        acrescenta1 = acrescenta2 = acrescenta3 = 0.0;
        indiceAcrescenta1 = indiceAcrescenta2 = indiceAcrescenta3 = -1;
        nx = -sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;
        dx1 = M->dx[i];
        dx2 = (i==M->Nx-1) ? M->dx[i] : M->dx[i+1];
        coef1 = (j==M->Ny-1) ? (2.0*beta*M->dt*nx*ny)/(M->Re*0.5*(dx1+dx2)*(M->dy[j])) : (2.0*beta*M->dt*nx*ny)/(M->Re*0.25*(dx1+dx2)*(M->dy[j]+M->dy[j+1]));
        coef2 = (2.0*beta*M->dt*(nx*nx - ny*ny))/(M->Re);
	coef3 = cilindrico*(2.0*beta*M->dt*((ny*ny)/rI))/(M->Re);
        h1 = (i==0) ? M->dx[i] : 0.5*(M->dx[i-1] + M->dx[i]);
        h2 = (i==M->Nx-1) ? M->dx[i] : 0.5*(M->dx[i] + M->dx[i+1]);


        //Phi(i, j)
        valor = ( - 2.0*coef1 ) + ( coef2*2.0/(h1*h2) ) - ( coef3/h2 ) - 1.0;
        if( M->celulas[i][j]!=EMPTY ) {
            indiceAcrescenta1 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j];
        }
//        else DeuProblema("1baixoesq: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i, j+1)
        valor = 2.0*coef1;
        if( (j!=M->Ny-1) && M->celulas[i][j+1]!=EMPTY ) {
            indiceAcrescenta3 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i][j+1];
        }
//        else {
//            if( M->pontosV[i][j+1].tipoBoundary == NOSLIP )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i][j+1].tipoBoundary == INFLOW )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i][j+1].tipoBoundary == NEUMANN )
//                acrescenta1 += - valor;
//            else DeuProblema("2baixoesq: Problema Contorno Poisson %d %d\n", i, j);
//        }

        //Phi(i+1, j)
        valor = ( 2.0*coef1 ) - ( coef2*2.0/(h2*(h1+h2)) ) + ( coef3/h2 );
        if( (i!=M->Nx-1) && M->celulas[i+1][j]!=EMPTY ) {
            indiceAcrescenta2 = qtd;
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j];
        }
//        else if( CELULA_FLUIDO(i, j) ) {
//            if( M->pontosU[i+1][j].tipoBoundary == NOSLIP )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i+1][j].tipoBoundary == INFLOW )
//                acrescenta1 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i+1][j].tipoBoundary == NEUMANN )
//                acrescenta1 += - valor;
//            else DeuProblema("3.1baixoesq: Problema Contorno Poisson %d %d\n", i, j);
//        }
//        else DeuProblema("3baixoesq: Problema Contorno Poisson %d %d\n", i, j);

        //Phi(i+1, j+1)
        valor = - 2.0*coef1;
        if( (j!=M->Ny-1) && (i!=M->Nx-1) && M->celulas[i+1][j+1]!=EMPTY ) {
            valores[qtd] = valor;
            colunas[qtd++] = M->matIndicesP[i+1][j+1];
        }
//        else if( (i!=M->Nx-1) && CELULA_FLUIDO(i+1, j) ) { //Joga pra baixo
//            if( M->pontosV[i+1][j+1].tipoBoundary == NOSLIP )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i+1][j+1].tipoBoundary == INFLOW )
//                acrescenta2 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosV[i+1][j+1].tipoBoundary == NEUMANN )
//                acrescenta2 += -valor;
//            else {
//                DesenhaMalhaVTK(*M, 0);
//                ImprimeInterfaceVTK(*M, 100000);
//                DeuProblema("4baixoesq: Problema Contorno Poisson %d %d\n", i, j);
//            }
//        }
//        else if( (j!=M->Ny-1) && CELULA_FLUIDO(i, j+1) ) {
//            if( M->pontosU[i+1][j+1].tipoBoundary == NOSLIP )
//                acrescenta3 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i+1][j+1].tipoBoundary == INFLOW )
//                acrescenta3 += valor; //dirichlet na velocidade vira neumann na Psi
//            else if( M->pontosU[i+1][j+1].tipoBoundary == NEUMANN )
//                acrescenta3 += - valor;
//            else DeuProblema("4.1baixoesq: Problema Contorno Poisson %d %d\n", i, j);
//        }

        valores[indiceAcrescenta1] += acrescenta1;
        valores[indiceAcrescenta2] += acrescenta2;
        valores[indiceAcrescenta3] += acrescenta3;
    }
    else if( DEGENERADO_CIMA_BAIXO(i, j) ) {
        qtd = 0;

        //Apenas colocando psi(i, j)=0 neste caso
        valores[qtd] = 1.0;
        colunas[qtd++] = M->matIndicesP[i][j];
    }
    else if( DEGENERADO_ESQ_DIR(i, j) ) {
        qtd = 0;

        //Apenas colocando psi(i, j)=0 neste caso
        valores[qtd] = 1.0;
        colunas[qtd++] = M->matIndicesP[i][j];
    }
    else DeuProblema("\n Problema Linha Contorno Poisson: CASOS SUPERFICIE %d %d\n", i, j);

    MatSetValues(A, 1, &Linha, qtd, colunas, valores, INSERT_VALUES);
    return 0;
}

double LadoDireitoContornoPoisson(MALHA M, int i, int j, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, double **ViscosityNewtGen)
{
    static double valor = 0;
    static double uDirBaixo = 0.0, uDirCima = 0.0, uEsqBaixo = 0.0, uEsqCima = 0.0;
    static double vDirBaixo = 0.0, vDirCima = 0.0, vEsqBaixo = 0.0, vEsqCima = 0.0;
    static double h1, h2, velAnt, velCentro, velProx;
    static double coef1, coef2, nx, ny;
    static char direcaoDerivada = 'a';

    //double rImenosMeio, rImaisMeio, rI, 
    double rI, cilindrico;
    double coef3;
    cilindrico = (M.tipoCoord==AXI_CILINDRICO) ? 1.0 : 0.0;

    ///Colocando Psi=0 apenas para testes e debugs as vezes. NAO DEIXAR ISSO!
//    if( LADO_VERTICAL_DIREITA_func(i, j, &M) ) {
//        PrintDebug("BAM22\n");
//        valor = 0.0;
//        return valor;
//    }

    //Valores da coordenada r usados no caso axisimetrico
    if( M.tipoCoord==AXI_CILINDRICO ) {
    	rI = M.r[i];
        //rImenosMeio = M.r[i] - 0.5*M.dx[i];
        //rImaisMeio = M.r[i] + 0.5*M.dx[i];
    }
    else
        rI = 1.0;


    //Calculo da normal
    if( LADO_HORIZONTAL_CIMA_func(i, j, &M) ) {
        nx = 0.0;
        ny = 1.0;
        direcaoDerivada = 'u';
    }
    else if( LADO_HORIZONTAL_BAIXO_func(i, j, &M) ) {
        nx = 0.0;
        ny = -1.0;
        direcaoDerivada = 'u';
    }
    else if( LADO_VERTICAL_DIREITA_func(i, j, &M) ) {
        nx = 1.0;
        ny = 0.0;
        direcaoDerivada = 'v';
    }
    else if( LADO_VERTICAL_ESQUERDA_func(i, j, &M) ) {
        nx = -1.0;
        ny = 0.0;
        direcaoDerivada = 'v';
    }
    else if( QUINA_CIMA_DIREITA_func(i, j, &M) ) {
        nx = sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;
        direcaoDerivada = 'u';
    }
    else if( QUINA_CIMA_ESQUERDA_func(i, j, &M) ) {
        nx = -sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;
        direcaoDerivada = 'u';
    }
    else if( QUINA_BAIXO_DIREITA_func(i, j, &M) ) {
        nx = sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;
        direcaoDerivada = 'u';
    }
    else if( QUINA_BAIXO_ESQUERDA_func(i, j, &M) ) {
        nx = -sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;
        direcaoDerivada = 'u';
    }
    else if( DEGENERADO_CIMA_BAIXO_func(i, j, &M) ) {
        //Vou apenas colocar psi(i, j)=0 nesse caso degenerado
        return 0.0;
    }
    else if( DEGENERADO_ESQ_DIR_func(i, j, &M) ) {
        //Vou apenas colocar psi(i, j)=0 nesse caso degenerado
        return 0.0;
    }
    else DeuProblema("\n Problema lado direito contorno poisson %d %d\n", i, j);

    /// === Pegando os valores de U e V que vao ser usados na formula
    if( M.pontosV[i][j].tipoBoundary==NOSLIP || M.pontosV[i][j].tipoBoundary==INFLOW )
        uDirBaixo = - U[i+1][j]; //Condicao de dirichlet
    else if( M.pontosV[i][j].tipoBoundary==NEUMANN )
        uDirBaixo = U[i+1][j]; //Condicao de neumann
    else if( M.pontosV[i][j].tipoBoundary==SIMETRIA )
        uDirBaixo = U[i+1][j]; //Condicao de neumann
    else if( M.pontosU[i+1][j-1].tipo==SURFACE || M.pontosU[i+1][j-1].tipo==FULL || M.pontosU[i+1][j-1].tipoBoundary==NOSLIP || M.pontosU[i+1][j-1].tipoBoundary==NEUMANN || M.pontosU[i+1][j-1].tipoBoundary==INFLOW || M.pontosU[i+1][j-1].tipoBoundary==SIMETRIA || M.pontosU[i+1][j-1].tipoBoundary==PERIODICIDADE )
        uDirBaixo = U[i+1][j-1];
    else DeuProblema("\n 1 - RHS: Contorno Poisson %d %d\n", i, j);

    if( M.pontosV[i][j].tipoBoundary==NOSLIP || M.pontosV[i][j].tipoBoundary==INFLOW )
        uEsqBaixo =  - U[i][j];
    else if( M.pontosV[i][j].tipoBoundary==NEUMANN )
        uEsqBaixo =  U[i][j];
    else if( M.pontosV[i][j].tipoBoundary==SIMETRIA )
        uEsqBaixo =  U[i][j];
    else if( M.pontosU[i][j-1].tipo==SURFACE || M.pontosU[i][j-1].tipo==FULL || M.pontosU[i][j-1].tipoBoundary==NOSLIP || M.pontosU[i][j-1].tipoBoundary==INFLOW || M.pontosU[i][j-1].tipoBoundary==SIMETRIA || M.pontosU[i][j-1].tipoBoundary==PERIODICIDADE || M.pontosU[i][j-1].tipoBoundary==EIXO_AXISSIMETRIA )
        uEsqBaixo = U[i][j-1];
    else DeuProblema("\n 2 - RHS: Contorno Poisson %d %d\n", i, j);

    if( M.pontosV[i][j+1].tipoBoundary==NOSLIP || M.pontosV[i][j+1].tipoBoundary==INFLOW )
        uDirCima = - U[i+1][j];
    else if( M.pontosV[i][j+1].tipoBoundary==NEUMANN )
        uDirCima = U[i+1][j];
    else if( M.pontosV[i][j+1].tipoBoundary==SIMETRIA )
        uDirCima = U[i+1][j];
    else if( M.pontosU[i+1][j+1].tipo==SURFACE || M.pontosU[i+1][j+1].tipo==FULL || M.pontosU[i+1][j+1].tipoBoundary==NOSLIP || M.pontosU[i+1][j+1].tipoBoundary==NEUMANN || M.pontosU[i+1][j+1].tipoBoundary==INFLOW || M.pontosU[i+1][j+1].tipoBoundary==SIMETRIA || M.pontosU[i+1][j+1].tipoBoundary==PERIODICIDADE )
        uDirCima = U[i+1][j+1];
    // Caso esquisito e raro que acontece soh quando ta numa quina de parede (tipo na simulacao gota-buraco)
    else if( M.pontosU[i+1][j].tipo==BOUNDARY && M.pontosU[i+1][j+1].tipo==EMPTY )
        uDirCima = U[i+1][j]; //Extrapolha pra baixo
    else {
//        uEsqCima = 0.0; //nao sei se isso eh bom...
//        ImprimeInterfaceVTK(M, 10000000);
//        DesenhaMalhaVTK(M, 0);
        DeuProblema("\n 3 - RHS: Contorno Poisson %d %d %d\n", i, j, M.pontosU[i+1][j+1].tipo);
    }

    if( M.pontosV[i][j+1].tipoBoundary==NOSLIP || M.pontosV[i][j+1].tipoBoundary==INFLOW )
        uEsqCima = - U[i][j];
    else if( M.pontosV[i][j+1].tipoBoundary==NEUMANN )
        uEsqCima = U[i][j];
    else if( M.pontosV[i][j+1].tipoBoundary==SIMETRIA )
        uEsqCima = U[i][j];
    else if( M.pontosU[i][j+1].tipo==SURFACE || M.pontosU[i][j+1].tipo==FULL || M.pontosU[i][j+1].tipoBoundary==NOSLIP || M.pontosU[i][j+1].tipoBoundary==INFLOW || M.pontosU[i][j+1].tipoBoundary==SIMETRIA || M.pontosU[i][j+1].tipoBoundary==PERIODICIDADE || M.pontosU[i][j+1].tipoBoundary==EIXO_AXISSIMETRIA )
        uEsqCima = U[i][j+1];
    else {
//        uEsqCima = 0.0; //(nao sei se isso eh uma boa ideia...)
        DeuProblema("\n 4 - RHS: Contorno Poisson %d %d\n", i, j);
    }

    if( M.pontosU[i+1][j].tipoBoundary==NOSLIP || M.pontosU[i+1][j].tipoBoundary==INFLOW )
        vDirCima = - V[i][j+1];
    else if( M.pontosU[i+1][j].tipoBoundary==NEUMANN || M.pontosU[i+1][j].tipoBoundary==SIMETRIA )
        vDirCima = V[i][j+1];
    else if( M.pontosU[i+1][j].tipoBoundary==PERIODICIDADE )
        vDirCima = V[0][j+1];
    else if( M.pontosV[i+1][j+1].tipo==SURFACE || M.pontosV[i+1][j+1].tipo==FULL || M.pontosV[i+1][j+1].tipoBoundary==NOSLIP || M.pontosV[i+1][j+1].tipoBoundary==INFLOW || M.pontosV[i+1][j+1].tipoBoundary==NEUMANN || M.pontosV[i+1][j+1].tipoBoundary==SIMETRIA )
        vDirCima = V[i+1][j+1];
    else DeuProblema("\n 5 - RHS: Contorno Poisson %d %d %d\n", i, j, M.pontosV[i+1][j+1].tipo);

    if( M.pontosU[i+1][j].tipoBoundary==NOSLIP || M.pontosU[i+1][j].tipoBoundary==INFLOW )
        vDirBaixo = - V[i][j];
    else if( M.pontosU[i+1][j].tipoBoundary==NEUMANN || M.pontosU[i+1][j].tipoBoundary==SIMETRIA )
        vDirBaixo = V[i][j];
    else if( M.pontosU[i+1][j].tipoBoundary==PERIODICIDADE )
        vDirBaixo = V[0][j];
    else if( M.pontosV[i+1][j].tipo==SURFACE || M.pontosV[i+1][j].tipo==FULL || M.pontosV[i+1][j].tipoBoundary==NOSLIP || M.pontosV[i+1][j].tipoBoundary==INFLOW || M.pontosV[i+1][j].tipoBoundary==NEUMANN || M.pontosV[i+1][j].tipoBoundary==SIMETRIA )
        vDirBaixo = V[i+1][j];
    // Caso esquisito e raro que acontece soh quando ta numa quina de parede (tipo na simulacao gota-buraco)
    else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i+1][j].tipo==EMPTY )
        vDirBaixo = V[i][j]; //Extrapolei pra esquerda (valor da parede...)
    else {
        DesenhaMalhaVTK(M, 0);
        DeuProblema("\n 6 - RHS: Contorno Poisson %d %d\n", i, j);
    }


    if( ESQUERDA_V(i, j+1) && (M.pontosU[i][j].tipoBoundary==NOSLIP || M.pontosU[i][j].tipoBoundary==INFLOW) )
        vEsqCima = - V[i][j+1];
    else if( M.pontosU[i][j].tipoBoundary==SIMETRIA )
        vEsqCima = V[i][j+1];
    else if( M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA )
        vEsqCima = V[i][j+1];
    else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
        vEsqCima = V[M.Nx-1][j+1];
    else if( M.pontosV[i-1][j+1].tipo==SURFACE || M.pontosV[i-1][j+1].tipo==FULL || M.pontosV[i-1][j+1].tipoBoundary==NOSLIP || M.pontosV[i-1][j+1].tipoBoundary==INFLOW || M.pontosV[i-1][j+1].tipoBoundary==NEUMANN || M.pontosV[i-1][j+1].tipoBoundary==SIMETRIA )
        vEsqCima = V[i-1][j+1];
    else DeuProblema("\n 7 - RHS: Contorno Poisson %d %d\n", i, j);

    if( ESQUERDA_V(i, j) && (M.pontosU[i][j].tipoBoundary==NOSLIP || M.pontosU[i][j].tipoBoundary==INFLOW) )
        vEsqBaixo = - V[i][j];
    else if( M.pontosU[i][j].tipoBoundary==SIMETRIA )
        vEsqBaixo = V[i][j];
    else if( M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA )
        vEsqBaixo = V[i][j];
    else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
        vEsqBaixo = V[M.Nx-1][j];
    else if( M.pontosV[i-1][j].tipo==SURFACE || M.pontosV[i-1][j].tipo==FULL || M.pontosV[i-1][j].tipoBoundary==NOSLIP || M.pontosV[i-1][j].tipoBoundary==INFLOW || M.pontosV[i-1][j].tipoBoundary==NEUMANN || M.pontosV[i-1][j].tipoBoundary==SIMETRIA )
        vEsqBaixo = V[i-1][j];
    // Caso esquisito e raro que acontece soh quando ta numa quina de parede (tipo na simulacao gota-buraco)
    else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i-1][j].tipo==EMPTY )
        vEsqBaixo = V[i][j]; //Extrapolei pra direita (valor da parede...)
    else DeuProblema("\n 8 - RHS: Contorno Poisson %d %d %d\n\n", i, j, M.pontosV[i-1][j].tipo);

    /// We adapt the "beta" that multiplies the diffusion term depending on the model
    double beta = ( M.tipo_modelo==MODELO_EVPT ) ? M.eta_inf : M.beta;

    if( direcaoDerivada=='u' ) {
        coef1 = ( 2.0*beta*M.dt*(nx*nx - ny*ny) )/(M.Re);
        valor = -(coef1*(U[i+1][j] - U[i][j]))/(M.dx[i]);

	//Adicao para o caso axi-cilindrico
	coef3 = cilindrico*(2.0*beta*M.dt*((ny*ny)/rI))/(M.Re);
	valor += coef3*0.5*(U[i+1][j] + U[i][j]);
    }
    else if( direcaoDerivada=='v' ) {
        coef1 = ( 2.0*beta*M.dt*(ny*ny - nx*nx) )/(M.Re);
        valor = -(coef1*(V[i][j+1] - V[i][j]))/(M.dy[j]);

	//Adicao para o caso axi-cilindrico
	coef3 = cilindrico*(2.0*beta*M.dt*((nx*nx)/rI))/(M.Re);
	valor += coef3*0.5*(U[i+1][j] + U[i][j]);
    }
    else DeuProblema("\n 9 - RHS: Contorno Poisson %d %d\n", i, j);

    coef2 = (2.0*beta*M.dt*nx*ny)/(M.Re);


    //Derivada du/dy
    h1 = ( j==0 ) ? M.dy[j] : 0.5*(M.dy[j-1] + M.dy[j]);
    h2 = ( j==M.Ny-1 ) ? M.dy[j] : 0.5*(M.dy[j] + M.dy[j+1]);
    velAnt = 0.5*( uEsqBaixo + uDirBaixo );
    velProx = 0.5*( uEsqCima + uDirCima );
    velCentro = 0.5*( U[i][j] + U[i+1][j] );
    valor += -coef2*( ( -velAnt*h2/(h1*(h1+h2)) ) + ( velCentro*(h2-h1)/(h1*h2) ) + ( velProx*h1/(h2*(h1+h2)) ) );

    //Derivada dv/dx
    h1 = ( i==0 ) ? M.dx[i] : 0.5*(M.dx[i-1] + M.dx[i]);
    h2 = ( i==M.Nx-1 ) ? M.dx[i] : 0.5*(M.dx[i] + M.dx[i+1]);
    velAnt = 0.5*( vEsqCima + vEsqBaixo );
    velProx = 0.5*( vDirCima + vDirBaixo );
    velCentro = 0.5*( V[i][j] + V[i][j+1] );
    valor += -coef2*( ( -velAnt*h2/(h1*(h1+h2)) ) + ( velCentro*(h2-h1)/(h1*h2) ) + ( velProx*h1/(h2*(h1+h2)) ) );

    /// === Used only as a test thing for the channel poeiseuille flow when i want to impose a pressure at the boundary
    /// Usually this variable should be kept as ZERO
    static const int pressure_condition = 0;

    // PrintDebug("RETIRAR ISTO AQUI!!!! LadoDireitoPoisson\n\n");
    if( pressure_condition )
        valor = 0.0;

    //Pressao
    valor += M.dt*P[i][j];

    //Parte nao-newtoniana
    if( !velocidadeNewtoniana && !pressure_condition ){
	    //printf("velocidadeNewtoniana=%d\n", velocidadeNewtoniana);getchar();
        valor += (- M.dt/M.Re)*( nx*nx*Txx[i][j] + 2.0*nx*ny*Txy[i][j] + ny*ny*Tyy[i][j] );
    }

    //Parte da tsensao superficial
    double curvatura = 0.0;

    if( M.tensaoSuperficial ) {

//        if( 1 ) {
//            double curvaturaFF = CalculaCurvatura(&(M.interface), 0.5*(M.x[i] + M.x[i+1]), 0.5*(M.y[j] + M.y[j+1]), 1.5*M.minDxDy, 1.5*M.minDxDy, nx, ny);
////           // curvatura = CalculaCurvaturaPeloAnguloDeContato(&(M.interface), 0.5*(M.x[i] + M.x[i+1]), 0.5*(M.y[j] + M.y[j+1]), 1.5*M.minDxDy);
////
//            curvatura = CalculaCurvaturaBezier(&(M.interface), 0.5*(M.x[i] + M.x[i+1]), 0.5*(M.y[j] + M.y[j+1]));
////
//            PrintDebug("%d %d %lf %lf\n", i, j, curvatura, curvaturaFF);
//            getchar();
//        }
//        else
//        double curvaturaFreef;

//        int tsup_extremos = BuscaParametro(&M, "tsup_extremos");
        int tsup_extremos = 0;
        int aplicar_tensao = 0;
        if( tsup_extremos )
            aplicar_tensao = 1;
        else if( /*i!=0 && i!=1 &&*/ i!=M.Nx-2 && i!=M.Nx-1 &&
           j!=0 && j!=1 && j!=M.Ny-1 && j!=M.Ny-2 )
           aplicar_tensao = 1;


        if( aplicar_tensao ) {

            //Estrategia freeflow para curvatura
//            curvatura = CalculaCurvatura(&(M.interface), 0.5*(M.x[i] + M.x[i+1]), 0.5*(M.y[j] + M.y[j+1]), 1.5*M.minDxDy, 1.5*M.minDxDy, nx, ny);
//            curvatura = CalculaCurvaturaML(&(malha_temp.interface), 0.5*(M.x[i] + M.x[i+1]), 0.5*(M.y[j] + M.y[j+1]));

//            if( i==7 && j==1 ) {
//                DeuProblema("AAAA %lf %lf... %lf %lf\n", nx, ny, M.vetNx[i][j], M.vetNy[i][j]);
//            }

            curvatura = M.curvaturaCelulas[i][j];
            if( (nx*M.vetNx[i][j]+ny*M.vetNy[i][j]) < 0.0 ) {
                curvatura = -curvatura;


                /// Vou corrigir o sentido de vetNx pra usar na curvatura do caso axissimetrico (logo ali embaixo)
                M.vetNx[i][j] *= -1.0;
                M.vetNy[i][j] *= -1.0;
            }

            /// Adicionando o termo da curvatura pro caso axissimetrico
            if( (M.tipoCoord==AXI_CILINDRICO) ) {
                double normal_r = M.vetNx[i][j];
                double coord_r = M.r[i];

                // curvatura = 0.0;

//                if( j==0 )
//                    curvatura = fabs(curvatura);
                // if( i==0 && j<500 ) {
                //     PrintDebug("CURV %d... %d %d: %lf %lf\n", M.passoTemporal, i, j, curvatura, normal_r);
                // }
//                curvatura += ((coord_r<1e-8) ? 0.0 : fabs(normal_r)/coord_r);
                // if( i==0 )
                    // PrintDebug("%lf %lf: %lf %lf ... %lf %lf ... %lf\n", 0.5*(M.x[i] + M.x[i+1]), 0.5*(M.y[j] + M.y[j+1]), curvatura, curvatura + normal_r/coord_r, M.vetNx[i][j], M.vetNy[i][j], coord_r);
                curvatura += ((coord_r<1e-8) ? 0.0 : normal_r/coord_r);
                M.curvaturaCelulas[i][j] = curvatura;
                

//                PrintDebug("CUYRVATURA %d %d: %lf.... [%lf %lf]\n", i, j, curvatura, M.vetNx[i][j], M.vetNy[i][j]);



//                if( (normal_r < 0.0) && (j!=0) ) {
//                    DesenhaMalhaVTK(M, 0);
//                    ImprimeInterfaceVTK(M, 100000);
//                    DeuProblema("Celula [i, j] = [%d %d]:\n \t Normal = [%.15lf %.15lf] \n \t Coordenada r: %.15lf \n \t [kappa_1, kappa_2] = [%.15lf %.15lf]\n \t kappa = %.15lf\n", i, j, normal_r, M.vetNy[i][j], coord_r, curvatura - normal_r/coord_r, normal_r/coord_r, curvatura);
//                }
//                if( (i==0) && (M.passoTemporal%50==0) )
//                    PrintDebug("Celula [i, j] = [%d %d]:\n \t Normal = [%.15lf %.15lf] \n \t Coordenada r: %.15lf \n \t [kappa_1, kappa_2] = [%.15lf %.15lf]\n \t kappa = %.15lf\n", i, j, normal_r, M.vetNy[i][j], coord_r, curvatura - normal_r/coord_r, normal_r/coord_r, curvatura);
//
//                if( (j==M.Ny/2) )
//                    PrintDebug("Celula [i, j] = [%d %d]:\n \t Normal = [%.15lf %.15lf] \n \t Coordenada r: %.15lf \n \t [kappa_1, kappa_2] = [%.15lf %.15lf]\n \t kappa = %.15lf\n", i, j, normal_r, M.vetNy[i][j], coord_r, curvatura - normal_r/coord_r, normal_r/coord_r, curvatura);
            }

            /// colocando angulo de contato... teste...
            int usar_angulo_contato = 0;
            if( usar_angulo_contato && j==0 ) {
                int iS = -1, jS = 1;

                // Encontrando a celula SURFACE acima dessa
                if( M.celulas[i-1][j+1]==SURFACE )
                    iS = i - 1;
                else if( M.celulas[i][j+1]==SURFACE )
                    iS = i;
                else if( M.celulas[i+1][j+1]==SURFACE )
                    iS = i + 1;
                else {
                    // ImprimeInterfaceVTK(M, 1000000);
                    // DesenhaMalhaVTK(M, 0);
                    DeuProblema("Angulo de contato: nao encontrou a celula acima...\n");
                }

//                double pontoS_x = M.pontoCelula[iS][jS]->x;
                double pontoS_y = M.pontoCelula[iS][jS]->y;
//                double vetS_x = M.vetNx[iS][jS];
                double vetS_y = M.vetNy[iS][jS];
                double betaContato = pontoS_y - M.yMin;

                double angulo = BuscaParametro(&M, "contact_angle");
                angulo *= M_PI/180.0;
//                double vet1_x = sin(M_PI - angulo);
                double vet1_y = -cos(M_PI - angulo);

                /// Sobrescrevendo a curvatura
                // PrintDebug("BETA: %e\n", betaContato);
                curvatura = (vetS_y - vet1_y)/betaContato;


//                DeuProblema("CELULA %d, %d: %lf %lf ... %lf %lf ... %lf %lf... %lf\n", iS, jS, pontoS_x, pontoS_y, vetS_x, vetS_y, vet1_x, vet1_y, curvatura);
            }

        }
    }
    

    // if( i<=3 ) {
    //     PrintDebug("i: %d     curv=%lf   normal=[%lf %lf]\n", i, curvatura, M.vetNx[i][j], M.vetNy[i][j]);
    // }
    
    valor += -(M.dt)*curvatura/(M.weber);
    
    /// CANAL: colocando a pressao imposta no canal suplivre
    // static double p_in = 30.0, p_out = 0.0;
    // double x = 0.5*(M.x[i] + M.x[i+1]);
    // double pressure_suplivre = p_in - (x + 5.0)*(p_in - p_out)/10.0;
    // valor += - (M.dt)*pressure_suplivre;

    return valor;
}

void DeterminaVetorNormalDeCadaCelula(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, double RaioParticulas)
{
    double xc, yc, distancia;
    CELULA_SURFACE *celula;
    PONTO *p;
    CURVA *c;

    /// Escolhendo o raio de particulas em volta de cada centro de celula
    RaioParticulas = RaioParticulas*RaioParticulas;

    /// Inicializando a matriz de cada celula surface
    double valorInicial = 0.0;
    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        celula->menorDistancia = 1e+10;
        celula->matrizCurvatura = (double **)AlocaMatriz(3, 4, sizeof(double), &valorInicial);
    }

    /// Percorrendo cada ponto da interface e verificando qual esta mais proximo do centro de cada celula SURFACE
    for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            /// No caso axissimetrico, ignora pontos sobre o eixo de simetria
            if( (M->tipoCoord==AXI_CILINDRICO) && (p->x<1e-6) )
                continue;

            /// Ignora pontos fixos no contorno (por exemplo na simulacao CABER)
            if( p->fixo )
                continue;

            // Cada celula surface...
            for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
                xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
                yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
                distancia = ((xc - p->x)*(xc - p->x)) + ((yc - p->y)*(yc - p->y));

                if( distancia > celula->menorDistancia )
                    continue;

                celula->menorDistancia = distancia;
                celula->pontoEscolhido = p;
                M->pontoCelula[celula->i][celula->j] = p;
            }

        }
    }

    int QtdParticulas = 100;


    /// Construindo o sistema de cada celula
    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        int i;
        PONTO *pTras, *pFrente;
        int pAxiSimetriaFrente = 0, pAxiSimetriaTras = 0; // Usado para "espelhar" a interface no caso axissimetrico

        xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
        yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
        pTras = celula->pontoEscolhido;
        pFrente = celula->pontoEscolhido;
        double distanciaTotal = 0.0;
        for( i=0; i<QtdParticulas; i++ ) {
            PONTO *pTras2;
            PONTO *pFrente2;

            

            if( !pAxiSimetriaFrente )
                pFrente2 = pFrente->prox;
            else
                pFrente2 = pFrente->ant;

            if( !pAxiSimetriaTras )
                pTras2 = pTras->ant;
            else
                pTras2 = pTras->prox;

            /// Chegou no eixo de simetria. Vou comecar a espelhar pontos "fantasmas"
            if( (M->tipoCoord==AXI_CILINDRICO) && (pFrente2->x<=1e-6) ) {
                pAxiSimetriaFrente = 1;
                // continue;
            }
            if( (M->tipoCoord==AXI_CILINDRICO) && (pTras2->x<=1e-6) ) {
                pAxiSimetriaTras = 1;
                // continue;
            }

            // if( celula->i==0 ) {
            //     PrintDebug("PONTO %d %d: %lf %lf\n", celula->i, celula->j, pFrente2->x, pFrente2->y);
            // }

            distanciaTotal += sqrt( ((pTras->x - pTras2->x)*(pTras->x - pTras2->x)) + ((pTras->y - pTras2->y)*(pTras->y - pTras2->y)) );
            distanciaTotal += sqrt( ((pFrente->x - pFrente2->x)*(pFrente->x - pFrente2->x)) + ((pFrente->y - pFrente2->y)*(pFrente->y - pFrente2->y)) );

            if( distanciaTotal>=0.025 )
                break;
//            if( distanciaTotal>=0.1 )
//                break;
//            if( distanciaTotal>=3.0*M->minDxDy )
//                break;

            pTras = pTras2;
            pFrente = pFrente2;


            /// Freeflow - Ponto pra tras
            p = pTras;
            if( !(p->fixo) ) {
                
                double x = p->x;
                if( pAxiSimetriaTras )
                    x = -x;

                //lhs
                celula->matrizCurvatura[0][0] += x*x;
                celula->matrizCurvatura[0][1] += x;
                celula->matrizCurvatura[1][0] += x;
                celula->matrizCurvatura[1][1] += 1.0;
                //rhs
                celula->matrizCurvatura[0][2] += x*(p->y);
                celula->matrizCurvatura[1][2] += p->y;
            }



            /// Freeflow - Ponto pra frente
            p = pFrente;
            if( !(p->fixo) ) {

                double x = p->x;
                if( pAxiSimetriaFrente ) 
                    x = -x;

                //lhs
                celula->matrizCurvatura[0][0] += x*x;
                celula->matrizCurvatura[0][1] += x;
                celula->matrizCurvatura[1][0] += x;
                celula->matrizCurvatura[1][1] += 1.0;
                //rhs
                celula->matrizCurvatura[0][2] += x*(p->y);
                celula->matrizCurvatura[1][2] += p->y;
            }

        }

        // Adicionando tb o ponto central
        p = celula->pontoEscolhido;
        //lhs
        celula->matrizCurvatura[0][0] += (p->x)*(p->x);
        celula->matrizCurvatura[0][1] += p->x;
        celula->matrizCurvatura[1][0] += p->x;
        celula->matrizCurvatura[1][1] += 1.0;
        //rhs
        celula->matrizCurvatura[0][2] += (p->x)*(p->y);
        celula->matrizCurvatura[1][2] += p->y;
        
    }




    /// Percorrendo cada ponto da interface e verificando as celulas que usarao ele na sua curvatura e vetor normal
//    for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
//        for( p=c->nodeInicial; p!=NULL; p=p->prox ) {
//
//            /// No caso axissimetrico, ignora pontos sobre o eixo de simetria
//            if( (M->tipoCoord==AXI_CILINDRICO) && (p->x<1e-6) )
//                continue;
//
//            /// Ignora pontos fixos no contorno (por exemplo na simulacao CABER)
//            if( p->fixo )
//                continue;
//
//            // Cada celula surface...
//            for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
//                xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
//                yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
//                distancia = ((xc - p->x)*(xc - p->x)) + ((yc - p->y)*(yc - p->y));
//
//
//
//                if( distancia > RaioParticulas )
//                    continue;
//
//                //lhs
//                celula->matrizCurvatura[0][0] += (p->x)*(p->x);
//                celula->matrizCurvatura[0][1] += p->x;
//                celula->matrizCurvatura[1][0] += p->x;
//                celula->matrizCurvatura[1][1] += 1.0;
//                //rhs
//                celula->matrizCurvatura[0][2] += (p->x)*(p->y);
//                celula->matrizCurvatura[1][2] += p->y;
//            }
//        }
//    }

    /// Resolvendo o sistema de cada celula
    double **A, det, a;
    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        A = celula->matrizCurvatura;

        det = A[0][0]*A[1][1] - A[0][1]*A[1][0]; //Determinante
        if( fabs(det)<1e-9 ) {
            M->vetNx[celula->i][celula->j] = 1.0;
            M->vetNy[celula->i][celula->j] = 0.0;
            //DeuProblema("PROB: VETOR NORMAL MINIMOS QUADRADOS %lf %lf\n", Xc, Yc);
            continue;
        }

        a = (A[0][2]*A[1][1] - A[1][2]*A[0][1])/det;
        //b = (rhs[1]*A[0][0] - rhs[0]*A[1][0])/det;


        double norma = sqrt(a*a + 1.0);
        M->vetNx[celula->i][celula->j] = -a/norma;
        M->vetNy[celula->i][celula->j] = 1.0/norma;



        //Freeflow
    //    *ResultX = a/norma;
    //    *ResultY = -1.0/norma;
    }

    return;
}

void DeterminaCurvaturaDeCadaCelula(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, double RaioParticulas)
{
    double xc, yc, distancia, xi, eta;
    CELULA_SURFACE *celula;
    PONTO *p;
    CURVA *c;

    /// Escolhendo o raio de particulas em volta de cada centro de celula
    RaioParticulas = RaioParticulas*RaioParticulas;

    /// Inicializando a matriz de cada celula surface
    /// Estou supondo que foi alocado memoria previamente para esta matriz
    /// Geralmente estou alocando na funcao DeterminaVetorNormalDeCadaCelula
    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        celula->matrizCurvatura[0][0] = 0.0;
        celula->matrizCurvatura[0][1] = 0.0;
        celula->matrizCurvatura[0][2] = 0.0;
        celula->matrizCurvatura[1][0] = 0.0;
        celula->matrizCurvatura[1][1] = 0.0;
        celula->matrizCurvatura[1][2] = 0.0;
        celula->matrizCurvatura[2][0] = 0.0;
        celula->matrizCurvatura[2][1] = 0.0;
        celula->matrizCurvatura[2][2] = 0.0;

        //Vetor do lado direito
        celula->matrizCurvatura[0][3] = 0.0;
        celula->matrizCurvatura[1][3] = 0.0;
        celula->matrizCurvatura[2][3] = 0.0;
    }

    /// Percorrendo cada ponto da interface e verificando as celulas que usarao ele na sua curvatura
    for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            // Cada celula surface...
            for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
                xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
                yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
                distancia = ((xc - p->x)*(xc - p->x)) + ((yc - p->y)*(yc - p->y));

                if( distancia > RaioParticulas )
                    continue;

                //Freeflow
                xi = (p->x-xc)*M->vetNy[celula->i][celula->j] - (p->y-yc)*M->vetNx[celula->i][celula->j];
                eta = (p->x-xc)*M->vetNx[celula->i][celula->j] + (p->y-yc)*M->vetNy[celula->i][celula->j];

                //Elementos da matriz
                celula->matrizCurvatura[0][0] += xi*xi*xi*xi;
                celula->matrizCurvatura[0][1] += xi*xi*xi;
                celula->matrizCurvatura[0][2] += xi*xi;
                celula->matrizCurvatura[1][0] += xi*xi*xi;
                celula->matrizCurvatura[1][1] += xi*xi;
                celula->matrizCurvatura[1][2] += xi;
                celula->matrizCurvatura[2][0] += xi*xi;
                celula->matrizCurvatura[2][1] += xi;
                celula->matrizCurvatura[2][2] += 1.0;

                //Vetor do lado direito
                celula->matrizCurvatura[0][3] += xi*xi*eta;
                celula->matrizCurvatura[1][3] += xi*eta;
                celula->matrizCurvatura[2][3] += eta;
            }

        }
    }

    /// Resolvendo o sistema 3x3 de cada celula pra calcular a curvatura
    double **A, det, curvatura, d;
    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        A = celula->matrizCurvatura;

        det = A[0][0]*A[1][1]*A[2][2] + A[2][0]*A[0][1]*A[1][2] + A[1][0]*A[2][1]*A[0][2]
            - A[2][0]*A[1][1]*A[0][2] - A[0][0]*A[2][1]*A[1][2] - A[2][2]*A[1][0]*A[0][1]; //Determinante

        //PrintDebug("%e ", det);
        if( fabs(det)<1e-9 ) {
            curvatura = 0.0;
        }
        else {
            //Usando regra de cramer
            d = A[0][3]*A[1][1]*A[2][2] + A[2][3]*A[0][1]*A[1][2] + A[1][3]*A[2][1]*A[0][2]
                    - A[2][3]*A[1][1]*A[0][2] - A[0][3]*A[2][1]*A[1][2] - A[2][2]*A[1][3]*A[0][1]; //Determinante
            d = d/det;

            curvatura = -2.0*d;

            /// O sinal da curvatura vai ser ajustado na funcao LadoDireitoContornoPoisson
            /// Fiz isso pois eu ja estou calculando o vetor normal ruim la, entao nao queria calcular de novo aqui
//            if( (nx*NxCelulas+ny*NyCelulas) < 0.0 )
//                curvatura = -curvatura;

            M->curvaturaCelulas[celula->i][celula->j] = curvatura;
        }

        DesalocaMatriz((void**)celula->matrizCurvatura, 3, 4);
    }


}

void DeterminaCurvaturaDeCadaCelulaVersao2(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, int QtdParticulas)
{
    double xc, yc, distancia, xi, eta;
    CELULA_SURFACE *celula;
    PONTO *p;
    CURVA *c;

    // Apagar depois
    //FILE *arq = fopen("curvaturas.txt", "a");

    /// Inicializando...
    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        celula->menorDistancia = 1e+10;

        celula->matrizCurvatura[0][0] = 0.0;
        celula->matrizCurvatura[0][1] = 0.0;
        celula->matrizCurvatura[0][2] = 0.0;
        celula->matrizCurvatura[1][0] = 0.0;
        celula->matrizCurvatura[1][1] = 0.0;
        celula->matrizCurvatura[1][2] = 0.0;
        celula->matrizCurvatura[2][0] = 0.0;
        celula->matrizCurvatura[2][1] = 0.0;
        celula->matrizCurvatura[2][2] = 0.0;

        //Vetor do lado direito
        celula->matrizCurvatura[0][3] = 0.0;
        celula->matrizCurvatura[1][3] = 0.0;
        celula->matrizCurvatura[2][3] = 0.0;
    }

    /// Percorrendo cada ponto da interface e verificando qual esta mais proximo do centro de cada celula SURFACE
    for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            /// No caso axissimetrico, ignora pontos sobre o eixo de simetria
            if( (M->tipoCoord==AXI_CILINDRICO) && (p->x<1e-6) )
                continue;

            /// Ignora pontos fixos no contorno (por exemplo na simulacao CABER)
//            if( p->fixo )
//                continue;

            // Cada celula surface...
            for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
                xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
                yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
                distancia = ((xc - p->x)*(xc - p->x)) + ((yc - p->y)*(yc - p->y));

                if( distancia > celula->menorDistancia )
                    continue;

                celula->menorDistancia = distancia;
                celula->pontoEscolhido = p;
                M->pontoCelula[celula->i][celula->j] = p;
            }

        }
    }

    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        int i;
        PONTO *pTras, *pFrente;
        int pAxiSimetriaFrente = 0, pAxiSimetriaTras = 0;

        xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
        yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
        pTras = celula->pontoEscolhido;
        pFrente = celula->pontoEscolhido;
        double distanciaTotal = 0.0;
        for( i=0; i<QtdParticulas; i++ ) {
            PONTO *pTras2; //= (pTras->ant) ? pTras->ant : pTras->pontoAntCurva;
            PONTO *pFrente2; //= (pFrente->prox) ? pFrente->prox : pFrente->pontoProxCurva;

            if( !pAxiSimetriaFrente )
                pFrente2 = pFrente->prox;
            else
                pFrente2 = pFrente->ant;

            if( !pAxiSimetriaTras )
                pTras2 = pTras->ant;
            else
                pTras2 = pTras->prox;

            /// Chegou no eixo de simetria. Vou comecar a espelhar pontos "fantasmas"
            if( (M->tipoCoord==AXI_CILINDRICO) && (pFrente2->x<=1e-6) ) {
                pAxiSimetriaFrente = 1;
                // continue;
            }
            if( (M->tipoCoord==AXI_CILINDRICO) && (pTras2->x<=1e-6) ) {
                pAxiSimetriaTras = 1;
                // continue;
            }

            // if( celula->i==0 )
            //     PrintDebug("ponto frente %d %d: %lf %lf\n", celula->i, celula->j, pFrente2->x, pFrente2->y);

            distanciaTotal += sqrt( ((pTras->x - pTras2->x)*(pTras->x - pTras2->x)) + ((pTras->y - pTras2->y)*(pTras->y - pTras2->y)) );
            distanciaTotal += sqrt( ((pFrente->x - pFrente2->x)*(pFrente->x - pFrente2->x)) + ((pFrente->y - pFrente2->y)*(pFrente->y - pFrente2->y)) );

//            if( distanciaTotal>=0.025 )
//                break;
//            if( distanciaTotal>=0.1 )
//                break;
            if( distanciaTotal>=2.0*M->minDxDy )
                break;

            pTras = pTras2;
            pFrente = pFrente2;


            /// Freeflow - Ponto pra tras
            p = pTras;
            /// No caso axissimetrico, ignora pontos sobre o eixo de simetria
            if( /*!(p->fixo) &&*/ ((M->tipoCoord!=AXI_CILINDRICO) || (p->x>1e-6)) ) {
                double x = p->x;
                if( pAxiSimetriaTras )
                    x = -x;

                xi = (x-xc)*M->vetNy[celula->i][celula->j] - (p->y-yc)*M->vetNx[celula->i][celula->j];
                eta = (x-xc)*M->vetNx[celula->i][celula->j] + (p->y-yc)*M->vetNy[celula->i][celula->j];


                //Elementos da matriz
                celula->matrizCurvatura[0][0] += xi*xi*xi*xi;
                celula->matrizCurvatura[0][1] += xi*xi*xi;
                celula->matrizCurvatura[0][2] += xi*xi;
                celula->matrizCurvatura[1][0] += xi*xi*xi;
                celula->matrizCurvatura[1][1] += xi*xi;
                celula->matrizCurvatura[1][2] += xi;
                celula->matrizCurvatura[2][0] += xi*xi;
                celula->matrizCurvatura[2][1] += xi;
                celula->matrizCurvatura[2][2] += 1.0;
                //Vetor do lado direito
                celula->matrizCurvatura[0][3] += xi*xi*eta;
                celula->matrizCurvatura[1][3] += xi*eta;
                celula->matrizCurvatura[2][3] += eta;
            }



            /// Freeflow - Ponto pra frente
            p = pFrente;
            if( /*!(p->fixo) &&*/ ((M->tipoCoord!=AXI_CILINDRICO) || (p->x>1e-6)) ) {
                double x = p->x;
                if( pAxiSimetriaFrente )
                    x = -x;

                xi = (x-xc)*M->vetNy[celula->i][celula->j] - (p->y-yc)*M->vetNx[celula->i][celula->j];
                eta = (x-xc)*M->vetNx[celula->i][celula->j] + (p->y-yc)*M->vetNy[celula->i][celula->j];


                //Elementos da matriz
                celula->matrizCurvatura[0][0] += xi*xi*xi*xi;
                celula->matrizCurvatura[0][1] += xi*xi*xi;
                celula->matrizCurvatura[0][2] += xi*xi;
                celula->matrizCurvatura[1][0] += xi*xi*xi;
                celula->matrizCurvatura[1][1] += xi*xi;
                celula->matrizCurvatura[1][2] += xi;
                celula->matrizCurvatura[2][0] += xi*xi;
                celula->matrizCurvatura[2][1] += xi;
                celula->matrizCurvatura[2][2] += 1.0;
                //Vetor do lado direito
                celula->matrizCurvatura[0][3] += xi*xi*eta;
                celula->matrizCurvatura[1][3] += xi*eta;
                celula->matrizCurvatura[2][3] += eta;
            }

        }

//        PrintDebug("DISTANCIA TOTAL: %lf\n", distanciaTotal);

//        DeuProblema("CELULA %d %d\n\n", celula->i, celula->j);

        // Adicionando tb o ponto central
        p = celula->pontoEscolhido;
        xi = (p->x-xc)*M->vetNy[celula->i][celula->j] - (p->y-yc)*M->vetNx[celula->i][celula->j];
        eta = (p->x-xc)*M->vetNx[celula->i][celula->j] + (p->y-yc)*M->vetNy[celula->i][celula->j];

        //Elementos da matriz
        celula->matrizCurvatura[0][0] += xi*xi*xi*xi;
        celula->matrizCurvatura[0][1] += xi*xi*xi;
        celula->matrizCurvatura[0][2] += xi*xi;
        celula->matrizCurvatura[1][0] += xi*xi*xi;
        celula->matrizCurvatura[1][1] += xi*xi;
        celula->matrizCurvatura[1][2] += xi;
        celula->matrizCurvatura[2][0] += xi*xi;
        celula->matrizCurvatura[2][1] += xi;
        celula->matrizCurvatura[2][2] += 1.0;
        //Vetor do lado direito
        celula->matrizCurvatura[0][3] += xi*xi*eta;
        celula->matrizCurvatura[1][3] += xi*eta;
        celula->matrizCurvatura[2][3] += eta;
    }

    /// Resolvendo o sistema 3x3 de cada celula pra calcular a curvatura
    double **A, det, curvatura, d;
    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        A = celula->matrizCurvatura;

        det = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1]
            - A[0][2]*A[1][1]*A[2][0] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2]; //Determinante

        if( fabs(det)<1e-9 ) {
            xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
            yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
            curvatura = 0.0;
        }
        else {
            //Usando regra de cramer
            d = A[0][3]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][3] + A[0][2]*A[1][3]*A[2][1]
               - A[0][2]*A[1][1]*A[2][3] - A[0][3]*A[1][2]*A[2][1] - A[0][1]*A[1][3]*A[2][2]; //Determinante
            
            d = d/det;

            curvatura = -2.0*d;

            /// O sinal da curvatura vai ser ajustado na funcao LadoDireitoContornoPoisson
            /// Fiz isso pois eu ja estou calculando o vetor normal ruim la, entao nao queria calcular de novo aqui
//            if( (nx*NxCelulas+ny*NyCelulas) < 0.0 )
//                curvatura = -curvatura;
        }

        M->curvaturaCelulas[celula->i][celula->j] = curvatura;
        /// Se for axissimetrico tem que adicionar mais um termo na curvatura
        /// (isso vai ser feito na funcao LadoDireitoContornoPoisson)

        // Apagar depois
        //fprintf(arq, "%f ", curvatura);
    

        DesalocaMatriz((void**)celula->matrizCurvatura, 3, 4);
    }
    // Apagar depois
    //fprintf(arq, "\n\n");
    //fclose(arq);


    /// Corrigindo o sinal da curvatura, caso necessario
    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
        int i = celula->i;
        int j = celula->j;
        double nx = 0.0, ny = 0.0;

        //Calculo da normal pela classificacao de celulas
        if( LADO_HORIZONTAL_CIMA_func(i, j, M) ) {
            nx = 0.0;
            ny = 1.0;
        }
        else if( LADO_HORIZONTAL_BAIXO_func(i, j, M) ) {
            nx = 0.0;
            ny = -1.0;
        }
        else if( LADO_VERTICAL_DIREITA_func(i, j, M) ) {
            nx = 1.0;
            ny = 0.0;
        }
        else if( LADO_VERTICAL_ESQUERDA_func(i, j, M) ) {
            nx = -1.0;
            ny = 0.0;
        }
        else if( QUINA_CIMA_DIREITA_func(i, j, M) ) {
            nx = sqrt(2.0)/2.0;
            ny = sqrt(2.0)/2.0;
        }
        else if( QUINA_CIMA_ESQUERDA_func(i, j, M) ) {
            nx = -sqrt(2.0)/2.0;
            ny = sqrt(2.0)/2.0;
        }
        else if( QUINA_BAIXO_DIREITA_func(i, j, M) ) {
            nx = sqrt(2.0)/2.0;
            ny = -sqrt(2.0)/2.0;
        }
        else if( QUINA_BAIXO_ESQUERDA_func(i, j, M) ) {
            nx = -sqrt(2.0)/2.0;
            ny = -sqrt(2.0)/2.0;
        }
        else if( DEGENERADO_CIMA_BAIXO_func(i, j, M) ) {
            //Vou apenas colocar psi(i, j)=0 nesse caso degenerado
//            return 0.0;
        }
        else if( DEGENERADO_ESQ_DIR_func(i, j, M) ) {
            //Vou apenas colocar psi(i, j)=0 nesse caso degenerado
//            return 0.0;
        }
        else DeuProblema("\n CalculaCurvatura: problema classificacao celulas %d %d\n", i, j);

        if( (nx*M->vetNx[i][j]+ny*M->vetNy[i][j]) < 0.0 ) {
            M->curvaturaCelulas[celula->i][celula->j] *= - 1.0;

            /// Vou corrigir o sentido de vetNx pra usar na curvatura do caso axissimetrico (logo ali embaixo)
            M->vetNx[i][j] *= -1.0;
            M->vetNy[i][j] *= -1.0;
        }
    }

}

PONTO *pontoTeste = NULL;

void DeterminaCurvaturaDeCadaCelulaVersaoSplines(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, int QtdParticulas)
{
    double xc, yc, distancia;
    CELULA_SURFACE *celula;
    PONTO *p;
    CURVA *c;

    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox )
        celula->menorDistancia = 1e+10;


    /// Percorrendo cada ponto da interface e verificando qual esta mais proximo do centro de cada celula SURFACE
    for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            // Cada celula surface...
            for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
                xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
                yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
                distancia = ((xc - p->x)*(xc - p->x)) + ((yc - p->y)*(yc - p->y));

                if( distancia > celula->menorDistancia )
                    continue;

                celula->menorDistancia = distancia;
                celula->pontoEscolhido = p;
            }

        }
    }

    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
//        if( M->passoTemporal==300 && (celula->i==51 && celula->j==30) )
//            pontoTeste = celula->pontoEscolhido;

        M->curvaturaCelulas[celula->i][celula->j] = CalculaCurvaturaDaParticulaSpline(M, QtdParticulas, celula->pontoEscolhido);
        M->curvaturaCelulas[celula->i][celula->j] = fabs(M->curvaturaCelulas[celula->i][celula->j]);
    }


    return;
}

void DeterminaCurvaturaDeCadaCelulaVersaoBezier(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, int QtdParticulas)
{
    double xc, yc, distancia;
    CELULA_SURFACE *celula;
    PONTO *p;
    CURVA *c;

    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox )
        celula->menorDistancia = 1e+10;


    /// Percorrendo cada ponto da interface e verificando qual esta mais proximo do centro de cada celula SURFACE
    for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            // Cada celula surface...
            for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
                xc = 0.5*(M->x[celula->i] + M->x[celula->i + 1]);
                yc = 0.5*(M->y[celula->j] + M->y[celula->j + 1]);
                distancia = ((xc - p->x)*(xc - p->x)) + ((yc - p->y)*(yc - p->y));

                if( distancia > celula->menorDistancia )
                    continue;

                celula->menorDistancia = distancia;
                celula->pontoEscolhido = p;
            }

        }
    }

    for( celula=ListaSurf->prim; celula!=NULL; celula=celula->prox ) {
//        if( M->passoTemporal==300 && (celula->i==51 && celula->j==30) )
//            pontoTeste = celula->pontoEscolhido;

        M->curvaturaCelulas[celula->i][celula->j] = CalculaCurvaturaBezier(M, QtdParticulas, celula->pontoEscolhido);
        M->curvaturaCelulas[celula->i][celula->j] = fabs(M->curvaturaCelulas[celula->i][celula->j]);
    }


    return;
}

double CalculaCurvaturaDaParticulaSpline(MALHA *M, int QtdParticulas, PONTO *PontoEscolhido)
{
    int i;

    // Andando pra tras nas particulas ate chegar na particula inicial
    PONTO *p = PontoEscolhido;
    for( i=0; i<QtdParticulas; i++ )
        p = p->ant;

    PONTO *pInicial = p;

    // Alocando memoria para os vetores
    double *x, *y, *t_vec;
    int totalPontos = 2*QtdParticulas + 1;
    x = (double *)malloc( totalPontos*sizeof(double) );
    y = (double *)malloc( totalPontos*sizeof(double) );
    t_vec = (double *)malloc( totalPontos*sizeof(double) );


    // Fazendo a parametrizacao com o comprimento de cada segmento
    double tEscolhido = -1.0;
    t_vec[0] = 0.0;
    for( i=1; i<totalPontos; i++ ) {
        PONTO *pProx = p->prox;

        t_vec[i] = t_vec[i-1] + sqrt( (p->x - pProx->x)*(p->x - pProx->x) + (p->y - pProx->y)*(p->y - pProx->y) );

        if( pProx==PontoEscolhido )
            tEscolhido = t_vec[i];

        p = pProx;
    }

    if( tEscolhido<0.0 )
        DeuProblema("tEscolhido negativo\n\n");

    // Percorrendo os pontos que serao utilizdos na interpolacao splines e jogando em vetores
    p = pInicial;
    for( i=0; i<totalPontos; i++ ) {
        x[i] = p->x;
        y[i] = p->y;

        p = p->prox;
    }

    //Calculando os coeficientes das splines
    POL3 *s_x, *s_y;
    s_x = (POL3 *)malloc((totalPontos)*sizeof(POL3)); //alocando memoria para os polinomios
    s_y = (POL3 *)malloc((totalPontos)*sizeof(POL3)); //alocando memoria para os polinomios

    SplineCubicaNotAKnot(t_vec, x, s_x, totalPontos-1); //Calculando os polinomios
    SplineCubicaNotAKnot(t_vec, y, s_y, totalPontos-1); //Calculando os polinomios

//    SplineCubicaUriAscher(t_vec, x, s_x, totalPontos-1); //Calculando os polinomios
//    SplineCubicaUriAscher(t_vec, y, s_y, totalPontos-1); //Calculando os polinomios


    if( pontoTeste ) {
        printf("T ESCOLHIDO: %lf\n", tEscolhido);
        printf("T TOTAL: %lf\n", t_vec[totalPontos-1]);
    }

    // Calculando as derivadas (ordem 1 e ordem 2) das splines
    double t = tEscolhido;
    double x_d1, x_d2, y_d1, y_d2;
    x_d1 = ValorSplineDerivada(t, t_vec, s_x, totalPontos-1);
    x_d2 = ValorSplineDerivada2(t, t_vec, s_x, totalPontos-1);
    y_d1 = ValorSplineDerivada(t, t_vec, s_y, totalPontos-1);
    y_d2 = ValorSplineDerivada2(t, t_vec, s_y, totalPontos-1);

//    p = pInicial;
//    for( i=0; i<totalPontos; i++ ) {
//
//        t = t_vec[i];
//        x_d1 = ValorSplineDerivada(t, t_vec, s_x, totalPontos-1);
//        x_d2 = ValorSplineDerivada2(t, t_vec, s_x, totalPontos-1);
//        y_d1 = ValorSplineDerivada(t, t_vec, s_y, totalPontos-1);
//        y_d2 = ValorSplineDerivada2(t, t_vec, s_y, totalPontos-1);
//
//        double f_derivada = cos(p->x);
//        double f_derivada2 = - sin(p->x);
//        double curv_exata = - f_derivada2/( pow( 1.0 + f_derivada*f_derivada, 1.5) );
//
//        printf("curv [%.15lf %.15lf]: %lf %lf\n", p->x, p->y, (x_d1*y_d2 - y_d1*x_d2)/pow( (x_d1*x_d1 + y_d1*y_d1), 1.5), curv_exata);
//
////        printf("[%.15lf %.15lf] [%.15lf %.15lf]\n", p->x, p->y, ValorSpline(t_vec[i], t_vec, s_x, totalPontos-1), ValorSpline(t_vec[i], t_vec, s_y, totalPontos-1));
//
////        printf("DERV [%lf %lf] [%lf %lf]\n", p->x, p->y, ValorSplineDerivada(t, t_vec, s_x, totalPontos-1), ValorSplineDerivada(t, t_vec, s_y, totalPontos-1));
//        p = (p->prox) ? p->prox : p->pontoProxCurva;
//    }

    if( pontoTeste ) {
        PrintDebug("ENCONTROU %lf\n", (x_d1*y_d2 - y_d1*x_d2)/pow( (x_d1*x_d1 + y_d1*y_d1), 1.5));

        p = pInicial;
        for( i=0; i<totalPontos; i++ ) {
            printf("PONTO %lf %lf\n", p->x, p->y);

            p = p->prox;
        }

        FILE *arq;
        char nomeArq[500];
        sprintf(nomeArq, "ArquivosPlot/VTK/%s/spline.vtk", M->nomeArquivo);
        arq = fopen(nomeArq, "wt");

        int qtdPontosPlot = 100;
        double delta_t = t_vec[totalPontos-1]/(qtdPontosPlot - 1);

        fprintf(arq, "# vtk DataFile Version 2.0\nStructured Grid Example\nASCII\n");
        fprintf(arq, "DATASET POLYDATA\n");
        fprintf(arq, "POINTS %d float\n", qtdPontosPlot);

        t = 0.0;
        for(i=0; i<qtdPontosPlot; i++) {
            double x_plot = ValorSpline(t, t_vec, s_x, totalPontos-1);
            double y_plot = ValorSpline(t, t_vec, s_y, totalPontos-1);
            fprintf(arq, "%e %e %e\n", x_plot, y_plot, 0.0);

            t += delta_t;
        }

        fprintf(arq, "POLYGONS %d %d\n", 1, qtdPontosPlot+1);
        fprintf(arq, "%d", qtdPontosPlot);
        for(i=0; i<qtdPontosPlot; i++)
            fprintf(arq, " %d", i);
        fprintf(arq, "\n");
        fclose(arq);


//        DeuProblema("parou\n");
        PrintDebug("imprimiu arquivo spline\n");
    }

    // Finalmente, a curvatura
    return (x_d1*y_d2 - y_d1*x_d2)/pow( (x_d1*x_d1 + y_d1*y_d1), 1.5);
}

double CalculaCurvaturaBezier(MALHA *M, int QtdParticulas, PONTO *PontoEscolhido)
{
    int i;

    // Andando pra tras nas particulas ate chegar na particula inicial
    PONTO *p = PontoEscolhido;
    for( i=0; i<QtdParticulas; i++ )
        p = p->ant;

    PONTO *pInicial = p;

    // Alocando memoria para os vetores
    double *x, *y;
    int totalPontos = 2*QtdParticulas + 1;
    x = (double *)malloc( totalPontos*sizeof(double) );
    y = (double *)malloc( totalPontos*sizeof(double) );


    // Percorrendo os pontos que serao utilizdos na interpolacao splines e jogando em vetores
    p = pInicial;
    for( i=0; i<totalPontos; i++ ) {
        x[i] = p->x;
        y[i] = p->y;

        p = p->prox;
    }

    double t = 0.5;
    double x_d1, x_d2, y_d1, y_d2;
    x_d1 = Dif_Bezier(x, t, totalPontos-1);
    x_d2 = Dif2_Bezier(x, t, totalPontos-1);
    y_d1 = Dif_Bezier(y, t, totalPontos-1);
    y_d2 = Dif2_Bezier(y, t, totalPontos-1);

//    if( 1/*(Yc<0.02)*/ ) {
//        double tTeste = 0.5;
//        PrintDebug("    PIVO  BEZIER: %lf %lf\n", pivo->x, pivo->y);
//        //PrintDebug("PONTO BEZIER: %lf %lf\n", Bezier(x, 0.0, qtdPontos-1), Bezier(y, 0.0, qtdPontos-1));
//        PrintDebug("    PONTO BEZIER: %lf %lf\n", Bezier(x, tTeste, qtdPontos-1), Bezier(y, tTeste, qtdPontos-1));
//        PrintDebug("    P_FINAL BEZIER: %lf %lf\n", Bezier(x, 1.0, qtdPontos-1), Bezier(y, 1.0, qtdPontos-1));
//       // getchar();
//    }


    if( pontoTeste ) {
        PrintDebug("ENCONTROU %lf\n", (x_d1*y_d2 - y_d1*x_d2)/pow( (x_d1*x_d1 + y_d1*y_d1), 1.5));

        p = pInicial;
        for( i=0; i<totalPontos; i++ ) {
            printf("PONTO %lf %lf\n", p->x, p->y);

            p = p->prox;
        }

        FILE *arq;
        char nomeArq[500];
        sprintf(nomeArq, "ArquivosPlot/VTK/%s/bezier.vtk", M->nomeArquivo);
        arq = fopen(nomeArq, "wt");

        int qtdPontosPlot = 100;
        double delta_t = 1.0/(qtdPontosPlot - 1);

        fprintf(arq, "# vtk DataFile Version 2.0\nStructured Grid Example\nASCII\n");
        fprintf(arq, "DATASET POLYDATA\n");
        fprintf(arq, "POINTS %d float\n", qtdPontosPlot);

        t = 0.0;
        for(i=0; i<qtdPontosPlot; i++) {
            double x_plot = Bezier(x, t, totalPontos-1);
            double y_plot = Bezier(y, t, totalPontos-1);
            fprintf(arq, "%e %e %e\n", x_plot, y_plot, 0.0);

            t += delta_t;
        }

        fprintf(arq, "POLYGONS %d %d\n", 1, qtdPontosPlot+1);
        fprintf(arq, "%d", qtdPontosPlot);
        for(i=0; i<qtdPontosPlot; i++)
            fprintf(arq, " %d", i);
        fprintf(arq, "\n");
        fclose(arq);


        DeuProblema("parou arquvo bezier\n");
        PrintDebug("imprimiu arquivo bezier\n");
    }

    free(x);
    free(y);

    //curvatura = fabs(x_d1*y_d2 - y_d1*x_d2)/pow( (x_d1*x_d1 + y_d1*y_d1), 3.0/2.0);
    double curvatura = (x_d1*y_d2 - y_d1*x_d2)/pow( (x_d1*x_d1 + y_d1*y_d1), 1.5);

    return curvatura;
}

void AlocaLinhaCondicaoContornoPoisson(MALHA *M, int i, int j, int indiceVetor, int iStart, int iEnd, int *d_nnz, int *o_nnz, int *qtd)
{
    static int i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0, i6 = 0;
    static int j1= 0, j2 = 0, j3 = 0, j4 = 0, j5 = 0, j6 = 0;
    int coluna;

    if( LADO_HORIZONTAL_CIMA(i, j) ) {
        i1 = i-1; j1 = j-1;
        i2 = i-1; j2 = j;
        i3 = i; j3 = j-1;
        i4 = i; j4 = j;
        i5 = i+1; j5 = j-1;
        i6 = i+1; j6 = j;
    }
    else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
        i1 = i-1; j1 = j;
        i2 = i-1; j2 = j+1;
        i3 = i; j3 = j;
        i4 = i; j4 = j+1;
        i5 = i+1; j5 = j;
        i6 = i+1; j6 = j+1;
    }
    else if( LADO_VERTICAL_DIREITA(i, j) ) {
        i1 = i-1; j1 = j-1;
        i2 = i-1; j2 = j;
        i3 = i-1; j3 = j+1;
        i4 = i; j4 = j-1;
        i5 = i; j5 = j;
        i6 = i; j6 = j+1;
    }
    else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
        i1 = i; j1 = j-1;
        i2 = i; j2 = j;
        i3 = i; j3 = j+1;
        i4 = i+1; j4 = j-1;
        i5 = i+1; j5 = j;
        i6 = i+1; j6 = j+1;
    }
    else if( QUINA_CIMA_DIREITA(i, j) ) {
        i1 = -1; j1 = -1;
        i2 = i-1; j2 = j-1;
        i3 = i-1; j3 = j;
        i4 = i; j4 = j-1;
        i5 = i; j5 = j;
        i6 = -1; j6 = -1;
    }
    else if( QUINA_CIMA_ESQUERDA(i, j) ) {
        i1 = -1; j1 = -1;
        i2 = i; j2 = j-1;
        i3 = i; j3 = j;
        i4 = i+1; j4 = j-1;
        i5 = i+1; j5 = j;
        i6 = -1; j6 = -1;
    }
    else if( QUINA_BAIXO_DIREITA(i, j) ) {
        i1 = -1; j1 = -1;
        i2 = i-1; j2 = j;
        i3 = i-1; j3 = j+1;
        i4 = i; j4 = j;
        i5 = i; j5 = j+1;
        i6 = -1; j6 = -1;
    }
    else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
        i1 = -1; j1 = -1;
        i2 = i; j2 = j+1;
        i3 = i; j3 = j;
        i4 = i+1; j4 = j;
        i5 = i+1; j5 = j+1;
        i6 = -1; j6 = -1;
    }
    else if( DEGENERADO_CIMA_BAIXO(i, j) ) {
        //Vou apenas colocar uma equacao psi(i, j) = 0 neste caso
        i1 = i; j1 = j;
        i2 = -1; j2 = -1;
        i3 = -1; j3 = -1;
        i4 = -1; j4 = -1;
        i5 = -1; j5 = -1;
        i6 = -1; j6 = -1;
    }
    else if( DEGENERADO_ESQ_DIR(i, j) ) {
        //Vou apenas colocar uma equacao psi(i, j) = 0 neste caso
        i1 = i; j1 = j;
        i2 = -1; j2 = -1;
        i3 = -1; j3 = -1;
        i4 = -1; j4 = -1;
        i5 = -1; j5 = -1;
        i6 = -1; j6 = -1;
    }
    else {
        DesenhaMalhaVTK(*M, 0);
        DeuProblema("\n Problema aloca matriz contorno poisson %d %d\n", i, j);
    }


    if( i1>=0 && i1<M->Nx && j1>=0 && j1<M->Ny && (M->pontosP[i1][j1].tipo==FULL || M->pontosP[i1][j1].tipo==SURFACE) ) {
        qtd[indiceVetor]++;
        coluna = M->matIndicesP[i1][j1];
        if( coluna>=iStart && coluna<iEnd )
            d_nnz[indiceVetor]++;
        else
            o_nnz[indiceVetor]++;
    }

    if( i2>=0 && i2<M->Nx && j2>=0 && j2<M->Ny && (M->pontosP[i2][j2].tipo==FULL || M->pontosP[i2][j2].tipo==SURFACE) ) {
        qtd[indiceVetor]++;
        coluna = M->matIndicesP[i2][j2];
        if( coluna>=iStart && coluna<iEnd )
            d_nnz[indiceVetor]++;
        else
            o_nnz[indiceVetor]++;
    }

    if( i3>=0 && i3<M->Nx && j3>=0 && j3<M->Ny && (M->pontosP[i3][j3].tipo==FULL || M->pontosP[i3][j3].tipo==SURFACE) ) {
        qtd[indiceVetor]++;
        coluna = M->matIndicesP[i3][j3];
        if( coluna>=iStart && coluna<iEnd )
            d_nnz[indiceVetor]++;
        else
            o_nnz[indiceVetor]++;
    }

    if( i4>=0 && i4<M->Nx && j4>=0 && j4<M->Ny && (M->pontosP[i4][j4].tipo==FULL || M->pontosP[i4][j4].tipo==SURFACE) ) {
        qtd[indiceVetor]++;
        coluna = M->matIndicesP[i4][j4];
        if( coluna>=iStart && coluna<iEnd )
            d_nnz[indiceVetor]++;
        else
            o_nnz[indiceVetor]++;
    }

    if( i5>=0 && i5<M->Nx && j5>=0 && j5<M->Ny && (M->pontosP[i5][j5].tipo==FULL || M->pontosP[i5][j5].tipo==SURFACE) ) {
        qtd[indiceVetor]++;
        coluna = M->matIndicesP[i5][j5];
        if( coluna>=iStart && coluna<iEnd )
            d_nnz[indiceVetor]++;
        else
            o_nnz[indiceVetor]++;
    }

    if( i6>=0 && i6<M->Nx && j6>=0 && j6<M->Ny && (M->pontosP[i6][j6].tipo==FULL || M->pontosP[i6][j6].tipo==SURFACE) ) {
        qtd[indiceVetor]++;
        coluna = M->matIndicesP[i6][j6];
        if( coluna>=iStart && coluna<iEnd )
            d_nnz[indiceVetor]++;
        else
            o_nnz[indiceVetor]++;
    }

    return;
}

void LaplacianoSupExplicito(MALHA *M, LISTA_CELULAS_SURFACE *CelulasSurf, double **U, double **V, double **P, double **Psi)
{
    CELULA_SURFACE *celula;
    int qtdIncognitas;
    Mat A;
    Vec sol, rhs;
    int i, j;

    double valoresMat[3];
    int colunas[3], linha;
    double valorRHS;
    double curvatura;


    qtdIncognitas = 0;
    /// Contando a qtd de incognitas do sistema
    for( celula=CelulasSurf->prim; celula!=NULL; celula=celula->prox )
        qtdIncognitas++;

    /// Alocando a matriz
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, qtdIncognitas, qtdIncognitas);
    MatSetFromOptions(A);
    MatSeqAIJSetPreallocation(A, 3, NULL);

    /// Alocando os vetores
    VecCreate(PETSC_COMM_WORLD, &sol);
	VecSetSizes(sol, PETSC_DECIDE, qtdIncognitas);
    VecSetFromOptions(sol);
    VecDuplicate(sol, &(rhs));

    /// Colocando os elementos da matriz e do rhs
    linha = 0;
    for( celula=CelulasSurf->prim; celula!=NULL; celula=celula->prox ) {
        i = celula->i;
        j = celula->j;

        if( linha==0 ) {
            colunas[0] = linha;
            colunas[1] = linha + 1;
            colunas[2] = qtdIncognitas-1;

            valoresMat[0] = (1.0/M->dt) + (4.0/(M->Re*M->dx[0]*M->dx[0]));
            valoresMat[1] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
            valoresMat[2] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        }
        else if( linha==qtdIncognitas-1 ) {
            colunas[0] = 0;
            colunas[1] = linha - 1;
            colunas[2] = linha;

            valoresMat[2] = (1.0/M->dt) + (4.0/(M->Re*M->dx[0]*M->dx[0]));
            valoresMat[0] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
            valoresMat[1] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        }
        else {
            colunas[0] = linha-1;
            colunas[1] = linha;
            colunas[2] = linha+1;

            // Fazendo as diagonais
            valoresMat[0] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
            valoresMat[1] = (1.0/M->dt) + (4.0/(M->Re*M->dx[0]*M->dx[0]));
            valoresMat[2] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        }





        double nx = 0.0, ny = 0.0;
        double u, v, u_n;
        double del_um = 0, um1, um2;

        if( LADO_VERTICAL_DIREITA(i, j) ) {
            nx = 1.0;
            ny = 0.0;

            u = 0.5*( U[i][j+1] + U[i+1][j+1] );
            v = 0.5*( V[i][j+1] + V[i][j+2] );
            um1 = - u*ny + v*nx;

            u = 0.5*( U[i][j-1] + U[i+1][j-1] );
            v = 0.5*( V[i][j] + V[i][j-1] );
            um2 = - u*ny + v*nx;

            del_um = (um1 - um2)/(2.0*M->dy[0]);
        }
        else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
            nx = -1.0;
            ny = 0.0;

            u = 0.5*( U[i][j+1] + U[i+1][j+1] );
            v = 0.5*( V[i][j+1] + V[i][j+2] );
            um1 = - u*ny + v*nx;

            u = 0.5*( U[i][j-1] + U[i+1][j-1] );
            v = 0.5*( V[i][j] + V[i][j-1] );
            um2 = - u*ny + v*nx;

            del_um = - (um1 - um2)/(2.0*M->dy[0]);
        }
        else if( LADO_HORIZONTAL_CIMA(i, j) ) {
            nx = 0.0;
            ny = 1.0;

            u = 0.5*( U[i+1][j] + U[i+2][j] );
            v = 0.5*( V[i+1][j] + V[i+1][j+1] );
            um1 = - u*ny + v*nx;

            u = 0.5*( U[i-1][j] + U[i][j] );
            v = 0.5*( V[i-1][j] + V[i-1][j+1] );
            um2 = - u*ny + v*nx;

            del_um = - (um1 - um2)/(2.0*M->dy[0]);
        }
        else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
            nx = 0.0;
            ny = -1.0;

            u = 0.5*( U[i+1][j] + U[i+2][j] );
            v = 0.5*( V[i+1][j] + V[i+1][j+1] );
            um1 = - u*ny + v*nx;

            u = 0.5*( U[i-1][j] + U[i][j] );
            v = 0.5*( V[i-1][j] + V[i-1][j+1] );
            um2 = - u*ny + v*nx;

            del_um = (um1 - um2)/(2.0*M->dy[0]);
        }
        else if( QUINA_CIMA_DIREITA(i, j) ) {
            nx = sqrt(2.0)/2.0;
            ny = sqrt(2.0)/2.0;


            u = 0.5*( U[i][j] + U[i-1][j] );
            v = 0.5*( V[i-1][j] + V[i-1][j+1] );
            um1 = - u*ny + v*nx;

            u = 0.5*( U[i][j-1] + U[i+1][j-1] );
            v = 0.5*( V[i][j] + V[i][j-1] );
            um2 = - u*ny + v*nx;

            del_um = (um1 - um2)/(2.0*M->dy[0]);
        }
        else if( QUINA_CIMA_ESQUERDA(i, j) ) {
            nx = -sqrt(2.0)/2.0;
            ny = sqrt(2.0)/2.0;

            u = 0.5*( U[i][j-1] + U[i+1][j-1] );
            v = 0.5*( V[i][j] + V[i][j-1] );
            um1 = - u*ny + v*nx;

            u = 0.5*( U[i+1][j] + U[i+2][j] );
            v = 0.5*( V[i+1][j] + V[i+1][j+1] );
            um2 = - u*ny + v*nx;

            del_um = (um1 - um2)/(2.0*M->dy[0]);
        }
        else if( QUINA_BAIXO_DIREITA(i, j) ) {
            nx = sqrt(2.0)/2.0;
            ny = -sqrt(2.0)/2.0;

            u = 0.5*( U[i][j+1] + U[i+1][j+1] );
            v = 0.5*( V[i][j+1] + V[i][j+2] );
            um1 = - u*ny + v*nx;

            u = 0.5*( U[i-1][j] + U[i][j] );
            v = 0.5*( V[i-1][j] + V[i-1][j+1] );
            um2 = - u*ny + v*nx;

            del_um = (um1 - um2)/(2.0*M->dy[0]);
        }
        else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
            nx = -sqrt(2.0)/2.0;
            ny = -sqrt(2.0)/2.0;

            u = 0.5*( U[i+1][j] + U[i+2][j] );
            v = 0.5*( V[i+1][j] + V[i+1][j+1] );
            um1 = - u*ny + v*nx;

            u = 0.5*( U[i][j+1] + U[i+1][j+1] );
            v = 0.5*( V[i][j+1] + V[i][j+2] );
            um2 = - u*ny + v*nx;

            del_um = (um1 - um2)/(2.0*M->dy[0]);
        }
        else DeuProblema("CASO INESPERADO LAP SURF.\n\n");

        u = 0.5*( U[i][j] + U[i+1][j] );
        v = 0.5*( V[i][j] + V[i][j+1] );
        u_n = u*nx + v*ny;

        // Fazendo a diagonal inferior
        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

        valorRHS = (2.0/M->Re)*( curvatura*u_n - del_um ) + (curvatura/M->weber) - P[i][j];

        MatSetValues(A, 1, &linha, 3, colunas, valoresMat, INSERT_VALUES);
        VecSetValues(rhs, 1, &linha, &valorRHS, INSERT_VALUES);
        linha++;
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);


    KSP kspSurf;

    InicializaSolverKSP(&kspSurf, A);
    KSPSolve(kspSurf, rhs, sol);

    linha = 0;
    for( celula=CelulasSurf->prim; celula!=NULL; celula=celula->prox ) {
        VecGetValues(sol, 1, &linha, &(Psi[celula->i][celula->j]));
        linha++;
    }

    VecDestroy(&rhs);
    VecDestroy(&sol);
    MatDestroy(&A);
    KSPDestroy(&kspSurf);


    //ImprimeMatrizPraMatlab(A, "debugs/matriz.m");
    //ImprimeVetorPraMatlab(rhs, "debugs/ladoDireito.m");
    //ImprimeVetorPraMatlab(sol, "debugs/sol.m");

    return;
}

void LaplacianoSupImplicitoMatriz(Mat A, MALHA *M, int i, int j, int Linha)
{
    double curvatura;
    int colunas[4];
    double valores[4];
    double nx, ny;


    if( LADO_VERTICAL_DIREITA(i, j) ) {
        nx = 1.0;
        ny = 0.0;

        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

        valores[0] = - ( 2.0*curvatura/(M->Re*M->dx[0]) );
        valores[1] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        valores[2] = (1.0/M->dt) + ( 4.0/(M->Re*M->dy[0]*M->dy[0]) ) + ( 2.0*curvatura/(M->Re*M->dx[0]) );
        valores[3] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));

        colunas[0] = M->matIndicesP[i-1][j];
        colunas[1] = M->matIndicesP[i][j-1];
        colunas[2] = M->matIndicesP[i][j];
        colunas[3] = M->matIndicesP[i][j+1];
    }
    else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
        nx = -1.0;
        ny = 0.0;

        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

        valores[3] = - ( 2.0*curvatura/(M->Re*M->dx[0]) );
        valores[2] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        valores[1] = (1.0/M->dt) + ( 4.0/(M->Re*M->dy[0]*M->dy[0]) ) + ( 2.0*curvatura/(M->Re*M->dx[0]) );
        valores[0] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));

        colunas[0] = M->matIndicesP[i][j-1];
        colunas[1] = M->matIndicesP[i][j];
        colunas[2] = M->matIndicesP[i][j+1];
        colunas[3] = M->matIndicesP[i+1][j];
    }
    else if( LADO_HORIZONTAL_CIMA(i, j) ) {
        nx = 0.0;
        ny = 1.0;

        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

        valores[1] = - ( 2.0*curvatura/(M->Re*M->dx[0]) );
        valores[0] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        valores[2] = (1.0/M->dt) + ( 4.0/(M->Re*M->dy[0]*M->dy[0]) ) + ( 2.0*curvatura/(M->Re*M->dx[0]) );
        valores[3] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));

        colunas[0] = M->matIndicesP[i-1][j];
        colunas[1] = M->matIndicesP[i][j-1];
        colunas[2] = M->matIndicesP[i][j];
        colunas[3] = M->matIndicesP[i+1][j];
    }
    else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
        nx = 0.0;
        ny = -1.0;

        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

        valores[2] = - ( 2.0*curvatura/(M->Re*M->dx[0]) );
        valores[0] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        valores[1] = (1.0/M->dt) + ( 4.0/(M->Re*M->dy[0]*M->dy[0]) ) + ( 2.0*curvatura/(M->Re*M->dx[0]) );
        valores[3] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));

        colunas[0] = M->matIndicesP[i-1][j];
        colunas[1] = M->matIndicesP[i][j];
        colunas[2] = M->matIndicesP[i][j+1];
        colunas[3] = M->matIndicesP[i+1][j];
    }
    else if( QUINA_CIMA_DIREITA(i, j) ) {
        nx = sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;

        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

        double delta = sqrt( (M->dx[0]*M->dx[0]) + (M->dy[0]*M->dy[0]) );

        valores[0] = - ( 2.0*curvatura/(M->Re*delta) );
        valores[1] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        valores[3] = (1.0/M->dt) + ( 4.0/(M->Re*M->dy[0]*M->dy[0]) ) + ( 2.0*curvatura/(M->Re*delta) );
        valores[2] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));

        colunas[0] = M->matIndicesP[i-1][j-1];
        colunas[1] = M->matIndicesP[i-1][j];
        colunas[2] = M->matIndicesP[i][j-1];
        colunas[3] = M->matIndicesP[i][j];
    }
    else if( QUINA_CIMA_ESQUERDA(i, j) ) {
        nx = -sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;

        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

        double delta = sqrt( (M->dx[0]*M->dx[0]) + (M->dy[0]*M->dy[0]) );

        valores[2] = - ( 2.0*curvatura/(M->Re*delta) );
        valores[0] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        valores[1] = (1.0/M->dt) + ( 4.0/(M->Re*M->dy[0]*M->dy[0]) ) + ( 2.0*curvatura/(M->Re*delta) );
        valores[3] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));

        colunas[0] = M->matIndicesP[i][j-1];
        colunas[1] = M->matIndicesP[i][j];
        colunas[2] = M->matIndicesP[i+1][j-1];
        colunas[3] = M->matIndicesP[i+1][j];
    }
    else if( QUINA_BAIXO_DIREITA(i, j) ) {
        nx = sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;

        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

        double delta = sqrt( (M->dx[0]*M->dx[0]) + (M->dy[0]*M->dy[0]) );

        valores[1] = - ( 2.0*curvatura/(M->Re*delta) );
        valores[0] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        valores[2] = (1.0/M->dt) + ( 4.0/(M->Re*M->dy[0]*M->dy[0]) ) + ( 2.0*curvatura/(M->Re*delta) );
        valores[3] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));

        colunas[0] = M->matIndicesP[i-1][j];
        colunas[1] = M->matIndicesP[i-1][j+1];
        colunas[2] = M->matIndicesP[i][j];
        colunas[3] = M->matIndicesP[i][j+1];
    }
    else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
        nx = -sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;

        curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);
        double delta = sqrt( (M->dx[0]*M->dx[0]) + (M->dy[0]*M->dy[0]) );

        valores[3] = - ( 2.0*curvatura/(M->Re*delta) );
        valores[1] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));
        valores[0] = (1.0/M->dt) + ( 4.0/(M->Re*M->dy[0]*M->dy[0]) ) + ( 2.0*curvatura/(M->Re*delta) );
        valores[2] = -(2.0/(M->Re*M->dx[0]*M->dx[0]));

        colunas[0] = M->matIndicesP[i][j];
        colunas[1] = M->matIndicesP[i][j+1];
        colunas[2] = M->matIndicesP[i+1][j];
        colunas[3] = M->matIndicesP[i+1][j+1];
    }
    else DeuProblema("CASO INESPERADO LAP SURF.\n\n");

    MatSetValues(A, 1, &Linha, 4, colunas, valores, INSERT_VALUES);
    return;
}

void LaplacianoSupImplicitoRHS(Vec Vetor, MALHA *M, int i, int j, int Linha, double **U, double **V, double **P)
{
    double nx = 0, ny = 0;
    double u = 0, v = 0, um1 = 0, um2 = 0, del_um = 0, u_n = 0;


    if( LADO_VERTICAL_DIREITA(i, j) ) {
        nx = 1.0;
        ny = 0.0;

        u = 0.5*( U[i][j+1] + U[i+1][j+1] );
        v = 0.5*( V[i][j+1] + V[i][j+2] );
        um1 = - u*ny + v*nx;

        u = 0.5*( U[i][j-1] + U[i+1][j-1] );
        v = 0.5*( V[i][j] + V[i][j-1] );
        um2 = - u*ny + v*nx;

        del_um = (um1 - um2)/(2.0*M->dy[0]);
    }
    else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
        nx = -1.0;
        ny = 0.0;

        u = 0.5*( U[i][j+1] + U[i+1][j+1] );
        v = 0.5*( V[i][j+1] + V[i][j+2] );
        um1 = - u*ny + v*nx;

        u = 0.5*( U[i][j-1] + U[i+1][j-1] );
        v = 0.5*( V[i][j] + V[i][j-1] );
        um2 = - u*ny + v*nx;

        del_um = - (um1 - um2)/(2.0*M->dy[0]);
    }
    else if( LADO_HORIZONTAL_CIMA(i, j) ) {
        nx = 0.0;
        ny = 1.0;

        u = 0.5*( U[i+1][j] + U[i+2][j] );
        v = 0.5*( V[i+1][j] + V[i+1][j+1] );
        um1 = - u*ny + v*nx;

        u = 0.5*( U[i-1][j] + U[i][j] );
        v = 0.5*( V[i-1][j] + V[i-1][j+1] );
        um2 = - u*ny + v*nx;

        del_um = - (um1 - um2)/(2.0*M->dy[0]);
    }
    else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
        nx = 0.0;
        ny = -1.0;

        u = 0.5*( U[i+1][j] + U[i+2][j] );
        v = 0.5*( V[i+1][j] + V[i+1][j+1] );
        um1 = - u*ny + v*nx;

        u = 0.5*( U[i-1][j] + U[i][j] );
        v = 0.5*( V[i-1][j] + V[i-1][j+1] );
        um2 = - u*ny + v*nx;

        del_um = (um1 - um2)/(2.0*M->dy[0]);
    }
    else if( QUINA_CIMA_DIREITA(i, j) ) {
        nx = sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;


        u = 0.5*( U[i][j] + U[i-1][j] );
        v = 0.5*( V[i-1][j] + V[i-1][j+1] );
        um1 = - u*ny + v*nx;

        u = 0.5*( U[i][j-1] + U[i+1][j-1] );
        v = 0.5*( V[i][j] + V[i][j-1] );
        um2 = - u*ny + v*nx;

        del_um = (um1 - um2)/(2.0*M->dy[0]);
    }
    else if( QUINA_CIMA_ESQUERDA(i, j) ) {
        nx = -sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;

        u = 0.5*( U[i][j-1] + U[i+1][j-1] );
        v = 0.5*( V[i][j] + V[i][j-1] );
        um1 = - u*ny + v*nx;

        u = 0.5*( U[i+1][j] + U[i+2][j] );
        v = 0.5*( V[i+1][j] + V[i+1][j+1] );
        um2 = - u*ny + v*nx;

        del_um = (um1 - um2)/(2.0*M->dy[0]);
    }
    else if( QUINA_BAIXO_DIREITA(i, j) ) {
        nx = sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;

        u = 0.5*( U[i][j+1] + U[i+1][j+1] );
        v = 0.5*( V[i][j+1] + V[i][j+2] );
        um1 = - u*ny + v*nx;

        u = 0.5*( U[i-1][j] + U[i][j] );
        v = 0.5*( V[i-1][j] + V[i-1][j+1] );
        um2 = - u*ny + v*nx;

        del_um = (um1 - um2)/(2.0*M->dy[0]);
    }
    else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
        nx = -sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;

        u = 0.5*( U[i+1][j] + U[i+2][j] );
        v = 0.5*( V[i+1][j] + V[i+1][j+1] );
        um1 = - u*ny + v*nx;

        u = 0.5*( U[i][j+1] + U[i+1][j+1] );
        v = 0.5*( V[i][j+1] + V[i][j+2] );
        um2 = - u*ny + v*nx;

        del_um = (um1 - um2)/(2.0*M->dy[0]);
    }
    else DeuProblema("CASO INESPERADO LAP SURF.\n\n");


    u = 0.5*( U[i][j] + U[i+1][j] );
    v = 0.5*( V[i][j] + V[i][j+1] );
    u_n = u*nx + v*ny;

    // Fazendo a diagonal inferior
    double curvatura = CalculaCurvatura(&(M->interface), 0.5*(M->x[i] + M->x[i+1]), 0.5*(M->y[j] + M->y[j+1]), 1.5*M->minDxDy, 1.5*M->minDxDy, nx, ny);

    double valorRHS = (2.0/M->Re)*( curvatura*u_n - del_um ) + (curvatura/M->weber) - P[i][j];
    VecSetValues(Vetor, 1, &Linha, &valorRHS, INSERT_VALUES);

    return;
}
