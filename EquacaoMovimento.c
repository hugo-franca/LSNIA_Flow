#include "EquacaoMovimento.h"

PetscErrorCode InicializaMatrizCoeficientesU(Mat A, int NumLinhas, MALHA M, int iStart, int iEnd)
{
    int i, j, linha, indiceVetor = 0, coluna;
    int d_nnz[iEnd - iStart], o_nnz[iEnd - iStart], qtd[iEnd - iStart];
    PetscErrorCode ierr;

    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesU[linha].i;
        j = M.indicesU[linha].j;

        if( (M.pontosU[i][j].tipo==EMPTY) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==INFLOW) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==NOSLIP) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==SLIP) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==SIMETRIA) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA) ||
            (M.pontosU[i][j].tipo==SURFACE) )
        {
            qtd[indiceVetor] = 1;
            d_nnz[indiceVetor] = 1;
            o_nnz[indiceVetor] = 0;
            indiceVetor++;
            continue;
        }
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==NEUMANN ) {
            qtd[indiceVetor] = 2;
            d_nnz[indiceVetor] = 1;
            o_nnz[indiceVetor] = 0;

            coluna = (M.pontosU[i][j].extremo==ESQUERDA) ? M.matIndicesU[i+1][j] : M.matIndicesU[i-1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;

            indiceVetor++;
            continue;
        }
//        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==PERIODICIDADE ) {
//            qtd[indiceVetor] = 2;
//            d_nnz[indiceVetor] = 1;
//            o_nnz[indiceVetor] = 0;
//
//            coluna = (M.pontosU[i][j].extremo==ESQUERDA) ? M.matIndicesU[i+1][j] : M.matIndicesU[i-1][j];
//            if( coluna>=iStart && coluna<iEnd )
//                d_nnz[indiceVetor]++;
//            else
//                o_nnz[indiceVetor]++;
//
//            indiceVetor++;
//            continue;
//        }



        o_nnz[indiceVetor] = 0;

        //Diagonal principal sempre vai ter
        d_nnz[indiceVetor] = 1;
        qtd[indiceVetor] = 1;

        //Verificando se o U(i-1, j)
        if( (M.pontosU[i][j].tipo!=BOUNDARY) || (M.pontosU[i][j].extremo==DIREITA) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesU[i-1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE ) {
            qtd[indiceVetor]++;
            o_nnz[indiceVetor]++;
        }



        //Verificando se o U(i, j-1) entra na matriz ou nao
        if( i!=M.Nx && !BAIXO_U(i, j) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesU[i][j-1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }
        else if( i==M.Nx && M.pontosV[i-1][j].tipo!=BOUNDARY ) { //Usando o de tras, caso esteja no canto do dominio
            qtd[indiceVetor]++;
            coluna = M.matIndicesU[i][j-1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }


        //Verificando U(i, j+1)
        if( i!=M.Nx && !CIMA_U(i, j) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesU[i][j+1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }
        else if( i==M.Nx && M.pontosV[i-1][j+1].tipo!=BOUNDARY ) {//Usando o de tras, caso esteja no canto do dominio
            qtd[indiceVetor]++;
            coluna = M.matIndicesU[i][j+1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }

        //Verificando U(i+1, j)
        if( (M.pontosU[i][j].tipo!=BOUNDARY) || (M.pontosU[i][j].extremo==ESQUERDA) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesU[i+1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE ) {
            qtd[indiceVetor]++;
            o_nnz[indiceVetor]++;
        }

        indiceVetor++;
    }


    ierr = MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A, 0, qtd); CHKERRQ(ierr);
    ierr = MatSeqSBAIJSetPreallocation(A, 1, 0, qtd); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode MontaMatrizCoeficientesU(Mat A, MALHA M, double **ViscosityNewtGen)
{
    int i, j, iV; //indices da malha
    int linha=0, colunas[5];
    int qtd = 1; //Qtd de elementos nao nulos em uma linha
    int centro;
    double valores[5], theta, acrescentaCentro;
    double rImenosMeio, rImaisMeio, rImenos3Meios, cilindrico;
    double hx1, hx2, hy1, hy2;
    double coefXesq, coefXdir, coefXcentro, coefYbaixo, coefYcima, coefYcentro;
    double beta;
    int iStart, iEnd;
    PetscErrorCode ierr;

    

    theta = 1.0;
    cilindrico = (M.tipoCoord==AXI_CILINDRICO) ? 1.0 : 0.0;

    ierr = MatGetOwnershipRange(A, &iStart, &iEnd); CHKERRQ(ierr);

    for( linha=iStart; linha<iEnd; linha++ ) {

        i = M.indicesU[linha].i;
        j = M.indicesU[linha].j;

        if( i==M.Nx )
            iV = i-1;
        else
            iV = i;

        //Colocando casos triviais: EMPTY, INFLOW, NOSLIP, NEUMANN
        if( M.pontosU[i][j].tipo == EMPTY ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==INFLOW) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==NOSLIP) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==SLIP) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==SIMETRIA) ||
            (M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA) ||
            (M.pontosU[i][j].tipo==SURFACE) )
        {
            valores[0] = 1.0;
            colunas[0] = linha; //Diagonal principal
            qtd = 1;

            ierr = MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES); CHKERRQ(ierr);
            continue;
        }
        else if(M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==NEUMANN) { //TESTE: Colocando uma equacao quase-trivial nos contornos NEUMANN

            if( M.pontosU[i][j].extremo==DIREITA ) {
                //Esquerda
                valores[0] = -1.0/M.dx[i-1];
                colunas[0] = M.matIndicesU[i-1][j];

                //Central
                valores[1] = 1.0/M.dx[i-1];
                colunas[1] = linha;
            }
            else if( M.pontosU[i][j].extremo==ESQUERDA ) {
                //Central
                valores[0] = -1.0/M.dx[i];
                colunas[0] = linha;

                //Direita
                valores[1] = 1.0/M.dx[i];
                colunas[1] = M.matIndicesU[i+1][j];
            }

            qtd = 2;
            ierr = MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES); CHKERRQ(ierr);
            continue;
        }
//        else if(M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==PERIODICIDADE) { //TESTE: Colocando uma equacao quase-trivial nos contornos NEUMANN
//
//            if( M.pontosU[i][j].extremo==DIREITA ) {
//                // Vamos adicionar a equacao "U(0, j) - U(Nx, j) = 0"
//                // Esquerda
//                valores[0] = 1.0;
//                colunas[0] = M.matIndicesU[1][j];
////                colunas[0] = M.matIndicesU[0][j];
//
//                // Direita
//                valores[1] = -1.0;
//                colunas[1] = linha;
//            }
//            else if( M.pontosU[i][j].extremo==ESQUERDA ) {
//                // Vamos adicionar a equacao "U(0, j) - U(Nx, j) = 0"
//                // Esquerda do dominio
//                valores[0] = 1.0;
//                colunas[0] = linha;
//
//                // Direita
//                valores[1] = -1.0;
//                colunas[1] = M.matIndicesU[M.Nx-1][j];
////                colunas[1] = M.matIndicesU[M.Nx][j];
//            }
//
//            qtd = 2;
//            ierr = MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES); CHKERRQ(ierr);
//            continue;
//        }

        if( ((i==0) || (i==M.Nx)) && (M.pontosU[i][j].tipoBoundary!=PERIODICIDADE) )
            DeuProblema("\n\n MatrizCoefU: NAO ERA PRA TER ACONTECIDO ISSO X \n\n");


        //Valores da coordenada r usados no caso axisimetrico
        //Corrigir os nomes dessasv ariaveis depois. mto ruim assim
        if( M.tipoCoord==AXI_CILINDRICO ) {
            rImenos3Meios = M.r[i-1] - 0.5*M.dx[i-1];
            rImenosMeio = M.r[i] - 0.5*M.dx[i];
            rImaisMeio = M.r[i] + 0.5*M.dx[i];
        }
        else
            rImenosMeio = rImaisMeio = rImenos3Meios = 1.0;

        /// We adapt the "beta" that multiplies the diffusion term depending on the model
        beta = ( M.tipo_modelo==MODELO_EVPT ) ? M.eta_inf : M.beta;

        hx1 = (i==0) ? M.dx[i] : M.dx[i-1];
        hx2 = (i==M.Nx) ? M.dx[i-1] : M.dx[i];
        hy1 = (j==0) ? M.dy[j] : 0.5*(M.dy[j-1] + M.dy[j]);
        hy2 = (j==M.Ny-1) ? M.dy[j] : 0.5*(M.dy[j] + M.dy[j+1]);

	/*
	//coeficientes com minha 'correcao' em 30/05/2025: troquei sinais nos coefXesq, coefXdir e coefXcentro, troquei rImenos3Meios e rImaisMeio por rImenosMeio.
        coefXesq = -beta*theta*(M.dt/M.Re)*( (2.0 - (cilindrico*hx2/rImenosMeio))/(hx1*(hx1+hx2)) );
        coefXdir = -beta*theta*(M.dt/M.Re)*( (2.0 + (cilindrico*hx1/rImenosMeio))/(hx2*(hx1+hx2)) );
        coefXcentro = -beta*theta*(M.dt/M.Re)*( -((2.0 - (cilindrico*(hx2-hx1)/rImenosMeio))/(hx1*hx2)) - (1.0/(rImenosMeio*rImenosMeio)) );
        coefYbaixo = -beta*theta*(M.dt/M.Re)*( 2.0/(hy1*(hy1+hy2)) );
        coefYcima = -beta*theta*(M.dt/M.Re)*( 2.0/(hy2*(hy1+hy2)) );
        coefYcentro = -beta*theta*(M.dt/M.Re)*( -2.0/(hy1*hy2) );
	*/
	//coeficientes da versao original do Hugo - conferi a discretização e está ok 16/06/2025. O hugo utiliza uma forma mais compacta de representar nabla^2 u.
        coefXesq = -beta*theta*(M.dt/M.Re)*( (2.0 + (cilindrico*hx2/rImenosMeio))*rImenos3Meios/(rImenosMeio*hx1*(hx1+hx2)) );
        coefXdir = -beta*theta*(M.dt/M.Re)*( (2.0 - (cilindrico*hx1/rImenosMeio))*rImaisMeio/(rImenosMeio*hx2*(hx1+hx2)) );
        coefXcentro = -beta*theta*(M.dt/M.Re)*( -(2.0 + (cilindrico*(hx2-hx1)/rImenosMeio))/(hx1*hx2) );
        coefYbaixo = -beta*theta*(M.dt/M.Re)*( 2.0/(hy1*(hy1+hy2)) );
        coefYcima = -beta*theta*(M.dt/M.Re)*( 2.0/(hy2*(hy1+hy2)) );
        coefYcentro = -beta*theta*(M.dt/M.Re)*( -2.0/(hy1*hy2) );
	
        qtd = 0;
        valores[0]=valores[1]=valores[2]=valores[3]=valores[4]=0.0;
        acrescentaCentro = 0.0;

        /// === Se estiver em uma boundary de periodicidade na direita, vai adicionar agora um ponto la na esquerda do dominio
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j].extremo==DIREITA ) {
            valores[qtd] = coefXdir;
            colunas[qtd++] = M.matIndicesU[1][j];
        }

        //Conferindo o da esquerda (i-1, j)
        // Soh eh pra entrar aqui no caso de PERIODICIDADE na esquerda
        if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
            if( M.pontosU[i][j].tipoBoundary!=PERIODICIDADE )
                DeuProblema("MontaMatrizU: PROBLEMA ESQUERDA\n");
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow) [esses casos entram la em cima ja]
        }
        else {
            valores[qtd] = coefXesq;
            colunas[qtd++] = M.matIndicesU[i-1][j];
        }

        //Conferindo o de baixo (i, j-1)
        if( BAIXO_U(iV, j) ) {
            if( M.pontosV[iV][j].tipoBoundary == NOSLIP )
                acrescentaCentro += -coefYbaixo; //Dirichlet
            else if( M.pontosV[iV][j].tipoBoundary == SLIP )
                acrescentaCentro += -coefYbaixo; //Dirichlet
            else if( M.pontosV[iV][j].tipoBoundary == INFLOW )
                acrescentaCentro += -coefYbaixo; //Dirichlet
            else if( M.pontosV[iV][j].tipoBoundary == NEUMANN )
                acrescentaCentro += coefYbaixo;
            else if( M.pontosV[iV][j].tipoBoundary == SIMETRIA )
                acrescentaCentro += coefYbaixo;
            else DeuProblema("MontaMatriz U: PROBLEMA BAIXO %d %d %d\n", i, j, M.pontosV[iV][j].tipoBoundary);
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
        }
        else {
            valores[qtd] = coefYbaixo;
            colunas[qtd++] = M.matIndicesU[i][j-1];
        }

        //Ponto central sempre entra na matriz
        centro = qtd++;
        valores[centro] = 1.0 + coefXcentro + coefYcentro;
        colunas[centro] = linha;

        //Conferindo o de cima (i, j+1)
        if( CIMA_U(iV, j) ) {
            if( (M.pontosV[iV][j+1].tipoBoundary == NOSLIP) )
                acrescentaCentro += -coefYcima;
            else if( (M.pontosV[iV][j+1].tipoBoundary == SLIP) )
                acrescentaCentro += -coefYcima;
            else if( M.pontosV[iV][j+1].tipoBoundary == INFLOW )
                acrescentaCentro += -coefYcima;
            else if( M.pontosV[iV][j+1].tipoBoundary == NEUMANN )
                acrescentaCentro += coefYcima;
            else if( M.pontosV[iV][j+1].tipoBoundary == SIMETRIA )
                acrescentaCentro += coefYcima;
            else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j+1].tipoBoundary==PERIODICIDADE ) {
                // Vou supor que eh zero esse elemento acima....
            }
            else DeuProblema("MontaMatriz U: PROBLEMA CIMA %d %d\n", i, j);
        }
        else {
            valores[qtd] = coefYcima;
            colunas[qtd++] = M.matIndicesU[i][j+1];
        }

        /// Soh deve entrar aqui se for caso de periodicidade
        //Conferindo o da direita (i+1, j)
        if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
            if( M.pontosU[i][j].tipoBoundary!=PERIODICIDADE )
                DeuProblema("MontaMatrizU: Problema Direita\n");
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow) [esses casos entram la em cima ja]
        }
        else {
            valores[qtd] = coefXdir;
            colunas[qtd++] = M.matIndicesU[i+1][j];
        }

        /// === Se estiver em uma boundary de periodicidade na esquerda, vai adicionar agora um ponto la na direita do dominio
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j].extremo==ESQUERDA ) {
            valores[qtd] = coefXesq;
            colunas[qtd++] = M.matIndicesU[M.Nx-1][j];
        }

        valores[centro] += acrescentaCentro;
        ierr = MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return 0;
}

PetscErrorCode VetorIndependenteU(Vec VetVetor,
                                  double **UVelho, double **VVelho, double **WVelho,
                                  double **UNovo, double **VNovo, double **WNovo,
                                  double **PVelho,
                                  double **Txx, double **Txy, double **Ttt,
                                  double **ViscosityFunction,
                                  int n, MALHA M)
{
    int linha = 0;
    int i, j, iV;
    double theta = 1.0, valor;
    double hy1, hy2;
//    int iStart, iEnd;
//    PetscErrorCode ierr;
//    ierr = VecGetOwnershipRange(VetVetor, &iStart, &iEnd); CHKERRQ(ierr);

//    for( linha=iStart; linha<iEnd; linha++ ) {
    for( linha=0; linha<M.qtdIncognitasU; linha++ ) {
        i = M.indicesU[linha].i;
        j = M.indicesU[linha].j;

        if( M.pontosU[i][j].tipo == EMPTY )
            printf("\n\nPROBLEMA\n\n");
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==INFLOW ) {
            valor = UVelho[i][j];
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==NOSLIP ) {
            valor = UVelho[i][j];
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==SLIP ) {
            valor = UVelho[i][j];
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==SIMETRIA ) {
            valor = UVelho[i][j]; //Simplesmente mantem o valor que ja esta la (vai ser zero)
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosU[i][j].tipo==SURFACE ) {
            valor = UNovo[i][j];
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA ) {
            valor = 0.0;
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==NEUMANN ) {
            valor = 0.0;
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else { //FULL!!

            /// If we are using a generalized newtonian model then the "div(T)" term will be a bit different from the viscoelastic case
            // double div_stress = (M.tipo_modelo==MODELO_VP) ? DivT_u_NewtGen(M, UVelho, VVelho, ViscosityFunction, i, j) : DivT_u(Txx, Txy, Ttt, i, j, M);
            double div_stress = DivT_u(Txx, Txy, Ttt, i, j, M);

            valor = UVelho[i][j] - M.dt*( Cx(UVelho, VVelho, WVelho, i, j, M.dt*n, M) + GradPx(PVelho, i, j, M) - (1.0/M.Re)*div_stress - ((1.0-theta)/M.Re)*Lx(UVelho, i, j, M.dt*n, M) );


            /// == Adicionando tensao superficial nas celulas FULL_SURFACE (em fase de testes...)
    //         M.tensaoSup_full = 1;
    //         if( M.tensaoSup_full ) {
    //             if( (i!=M.Nx && M.celulasParede[i][j]==SURFACE) || (i!=0 && M.celulasParede[i-1][j]==SURFACE) ) {
    //                 double normalX, normalY, curvatura;
    //                 CalculaNormalCurvaturaEquacaoMovimento(&M, M.x[i], 0.5*(M.y[j] + M.y[j+1]), &normalX, &normalY, &curvatura);
    //                 valor += M.dt*(1.0/M.weber)*curvatura*normalX;
    // //                    PrintDebug("VALOR X: %lf %lf %lf\n", M.dt*(1.0/M.weber)*curvatura*normalX, curvatura, normalX);
    //             }
    //         }

            //Condicao de contorno NOSLIP quando a parede horizontal superior esta se movimentando
            iV = (i==M.Nx) ? (i-1) : i;
            if( CIMA_U(iV, j) ) {

                /// We adapt the "beta" that multiplies the diffusion term depending on the model
                double beta = ( M.tipo_modelo==MODELO_EVPT ) ? M.eta_inf : M.beta;

                /// Determinando o tipo deste contorno
                /// Precisa olhar dois pontos pra nao bugar casos de quina tipo na contracao.... (gambiarra)
                int tipoDir = M.pontosV[iV][j+1].tipoBoundary;
                int tipoEsq = (iV==0) ? NAO_BOUNDARY : M.pontosV[iV-1][j+1].tipoBoundary;
                int movimentoDir = M.pontosV[iV][j+1].movimento;
                int movimentoEsq = (iV==0) ? SEM_MOVIMENTO : M.pontosV[iV-1][j+1].movimento;
                int tipo = NAO_BOUNDARY, movimento = SEM_MOVIMENTO;

                if( tipoEsq==NOSLIP || tipoEsq==SLIP ) {
                    tipo = tipoEsq;
                    movimento = movimentoEsq;
                }
                if( tipoDir==NOSLIP || tipoDir==SLIP ) {
                    tipo = tipoDir;
                    movimento = movimentoDir;
                }

                if( (tipo==NOSLIP || tipo==SLIP) && (movimento==MOVIMENTO_X) ) {
                    hy1 = (j==0) ? M.dy[j] : 0.5*(M.dy[j-1] + M.dy[j]);
                    hy2 = (j==M.Ny-1) ? M.dy[j] : 0.5*(M.dy[j] + M.dy[j+1]);

                    double coefYcima = -beta*theta*(M.dt/M.Re)*( 2.0/(hy2*(hy1+hy2)) );
                    valor += -coefYcima*2.0*(FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], (n+1)*M.dt, UVelho ));
                }
            }

            //Condicao de contorno NOSLIP quando a parede horizontal inferior esta se movimentando
            if( BAIXO_U(iV, j) ) {

                /// We adapt the "beta" that multiplies the diffusion term depending on the model
                double beta = ( M.tipo_modelo==MODELO_EVPT ) ? M.eta_inf : M.beta;

                /// Determinando o tipo deste contorno
                /// Precisa olhar dois pontos pra nao bugar casos de quina tipo na contracao.... (gambiarra)
                int tipoDir = M.pontosV[iV][j].tipoBoundary;
                int tipoEsq = (iV==0) ? NAO_BOUNDARY : M.pontosV[iV-1][j].tipoBoundary;
                int movimentoDir = M.pontosV[iV][j].movimento;
                int movimentoEsq = (iV==0) ? SEM_MOVIMENTO : M.pontosV[iV-1][j].movimento;
                int tipo = NAO_BOUNDARY, movimento = SEM_MOVIMENTO;

                if( tipoEsq==NOSLIP || tipoEsq==SLIP ) {
                    tipo = tipoEsq;
                    movimento = movimentoEsq;
                }
                if( tipoDir==NOSLIP || tipoDir==SLIP ) {
                    tipo = tipoDir;
                    movimento = movimentoDir;
                }

                if( (tipo==NOSLIP || tipo==SLIP) && (movimento==MOVIMENTO_X) ) {
                    hy1 = (j==0) ? M.dy[j] : 0.5*(M.dy[j-1] + M.dy[j]);
                    hy2 = (j==M.Ny-1) ? M.dy[j] : 0.5*(M.dy[j] + M.dy[j+1]);

                    double coefYbaixo = -beta*theta*(M.dt/M.Re)*( 2.0/(hy1*(hy1+hy2)) );
                    valor += -coefYbaixo*2.0*(FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], (n+1)*M.dt, UVelho ));
                }
            }


            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
            //Vet[linha] = UVelho[i][j] - M.dt*( Cx(UVelho, VVelho, i, j, M) + ((p-PVelho[i-1][j])/M.dx[0]) - ((1.0-theta)/M.Re)*Lx(UVelho, i, j, M) );
        }
    }

    VecAssemblyBegin(VetVetor);
    VecAssemblyEnd(VetVetor);
    return 0;
}

PetscErrorCode InicializaMatrizCoeficientesV(Mat A, int NumLinhas, MALHA M, int iStart, int iEnd)
{
    int i, j, linha, indiceVetor = 0, coluna;
    int d_nnz[iEnd - iStart], o_nnz[iEnd - iStart], qtd[iEnd - iStart];
    PetscErrorCode ierr;

    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesV[linha].i;
        j = M.indicesV[linha].j;

        if( M.pontosV[i][j].tipo == EMPTY ||
            (M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==INFLOW) ||
            (M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==NOSLIP) ||
            (M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==SLIP) ||
            (M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==SIMETRIA) ||
            (M.pontosV[i][j].tipo==SURFACE) )
        {
            qtd[indiceVetor] = 1;
            d_nnz[indiceVetor] = 1;
            o_nnz[indiceVetor] = 0;
            indiceVetor++;
            continue;
        }
        else if( M.pontosV[i][j].tipo==BOUNDARY &&
                (M.pontosV[i][j].tipoBoundary==NEUMANN) ) {
            qtd[indiceVetor] = 2;
            d_nnz[indiceVetor] = 1;
            o_nnz[indiceVetor] = 0;

            coluna = (M.pontosV[i][j].extremo==BAIXO) ? M.matIndicesV[i][j+1] : M.matIndicesV[i][j-1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;

            indiceVetor++;
            continue;
        }

        o_nnz[indiceVetor] = 0;

        //Diagonal principal sempre vai ter
        qtd[indiceVetor] = 1;
        d_nnz[indiceVetor] = 1;

        /// Adicionando um valor pro caso de periodicidade na esquerda
        if( j!=M.Ny && ESQUERDA_V(i, j) && M.pontosU[i][j].tipoBoundary==PERIODICIDADE ) {
            qtd[indiceVetor]++;
            o_nnz[indiceVetor]++;
        }
        /// Adicionando um valor pro caso de periodicidade na direita
        if( j!=M.Ny && DIREITA_V(i, j) && M.pontosU[i+1][j].tipoBoundary==PERIODICIDADE ) {
            qtd[indiceVetor]++;
            o_nnz[indiceVetor]++;
        }


        //Conferindo V(i-1, j)
        if( j!=M.Ny && !ESQUERDA_V(i, j) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesV[i-1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }
        else if( j==M.Ny && M.pontosU[i][j-1].tipo!=BOUNDARY ) { //Usando o de tras, caso esteja no canto do dominio
            qtd[indiceVetor]++;
            coluna = M.matIndicesV[i-1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }



        //O V(i, j-1) sempre entra (MUDAR PQ TA ERRADO SE FOR NEUMANN EMBAIXO)
        qtd[indiceVetor]++;
        coluna = M.matIndicesV[i][j-1];
        if( coluna>=iStart && coluna<iEnd )
            d_nnz[indiceVetor]++;
        else
            o_nnz[indiceVetor]++;



        //Conferindo V(i, j+1)
        if( (M.pontosV[i][j].tipo!=BOUNDARY) ) { //Vai ter que mudar isso pra considerar neumann embaixo
            qtd[indiceVetor]++;
            coluna = M.matIndicesV[i][j+1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }



        //Conferindo V(i+1, j)
        if( j!=M.Ny && !DIREITA_V(i, j) ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesV[i+1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }
        else if( j==M.Ny && M.pontosU[i+1][j-1].tipo!=BOUNDARY ) { //Usando o de tras, caso esteja no canto do dominio
            qtd[indiceVetor]++;
            coluna = M.matIndicesV[i+1][j];
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

PetscErrorCode MontaMatrizCoeficientesV(Mat A, MALHA M, double **ViscosityNewtGen)
{
    int i, j, jU; //indices da malha
    int linha=0, colunas[5];
    int qtd = 1; //Qtd de elementos nao nulos em uma linha
    int centro;
    double valores[5], acrescentaCentro, theta, cilindrico;
    double hx1, hx2, hy1, hy2;
    double coefXesq, coefXcentro, coefXdir;
    double coefYbaixo, coefYcentro, coefYcima;
    double beta;
    int iStart, iEnd;
    PetscErrorCode ierr;

    theta = 1.0;
    cilindrico = (M.tipoCoord==AXI_CILINDRICO) ? 1.0 : 0.0;

    ierr = MatGetOwnershipRange(A, &iStart, &iEnd); CHKERRQ(ierr);

    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesV[linha].i;
        j = M.indicesV[linha].j;

        if( j==M.Ny )
            jU = j-1;
        else
            jU = j;

        if( M.pontosV[i][j].tipo == EMPTY ||
            (M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==INFLOW) ||
            (M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==NOSLIP) ||
            (M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==SLIP) ||
            (M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==SIMETRIA) ||
            (M.pontosV[i][j].tipo==SURFACE) )
        {
            valores[0] = 1.0;
            colunas[0] = linha; //Diagonal principal
            qtd = 1;

            ierr = MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES); CHKERRQ(ierr);
            continue;
        }
        else if( M.pontosV[i][j].tipo==BOUNDARY &&
                (M.pontosV[i][j].tipoBoundary==NEUMANN) ) { //TESTE: Colocando uma equacao quase-trivial nos contornos NEUMANN

            //Cima ou baixo, depende de onde estamos no dominio
            if( M.pontosV[i][j].extremo==BAIXO ) {

                //Central
                valores[0] = - 1.0/M.dy[j];
                colunas[0] = linha;

                //Cima
                valores[1] = 1.0/M.dy[j];
                colunas[1] = M.matIndicesV[i][j+1];
            }
            else if( M.pontosV[i][j].extremo==CIMA ) {

                //Baixo
                valores[0] = - 1.0/M.dy[j-1];
                colunas[0] = M.matIndicesV[i][j-1];

                //Central
                valores[1] = 1.0/M.dy[j-1];
                colunas[1] = linha;
            }
            else
                DeuProblema("MatrizV: Problema \n\n");

            qtd = 2;
            ierr = MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES); CHKERRQ(ierr);
            continue;
        }

        if( (j==0) || (j==M.Ny) )
            DeuProblema("\n\n Matriz V: NAO ERA PRA TER ACONTECIDO ISSO!!!! %d %d %d\n\n", i, j, M.pontosV[i][j].tipo);

        /// We adapt the "beta" that multiplies the diffusion term depending on the model
        beta = ( M.tipo_modelo==MODELO_EVPT ) ? M.eta_inf : M.beta;

        hx1 = (i==0) ? M.dx[i] : 0.5*( M.dx[i-1] + M.dx[i] );
        hx2 = (i==M.Nx-1) ? M.dx[i] : 0.5*( M.dx[i] + M.dx[i+1] );
        hy1 = M.dy[j-1];
        hy2 = M.dy[j];

        coefXesq = -beta*theta*(M.dt/M.Re)*((-cilindrico*hx2/(M.r[i]*hx1*(hx1+hx2))) + (2.0/(hx1*(hx1+hx2))));
        coefXdir = -beta*theta*(M.dt/M.Re)*((cilindrico*hx1/(M.r[i]*hx2*(hx1+hx2))) + (2.0/(hx2*(hx1+hx2))));
        coefXcentro = -beta*theta*(M.dt/M.Re)*((cilindrico*(hx2-hx1)/(M.r[i]*hx1*hx2)) - (2.0/(hx1*hx2)));
        coefYbaixo = -beta*theta*(M.dt/M.Re)*(2.0/(hy1*(hy1+hy2)));
        coefYcima = -beta*theta*(M.dt/M.Re)*(2.0/(hy2*(hy1+hy2)));
        coefYcentro = -beta*theta*(M.dt/M.Re)*(-2.0/(hy1*hy2));


        qtd = 0;
        valores[0] = valores[1] = valores[2] = valores[3] = valores[4] = 0.0;
        acrescentaCentro = 0.0;


        /// === Primeiramente eu vou verificar se estamos em um contorno de PERIODICIDADE
        /// === Preciso fazer isso antes de construir os outros elementos devido a estrutura do petsc...
        if( DIREITA_V(i, jU) && M.pontosU[i+1][jU].tipoBoundary==PERIODICIDADE ) {
            // Adiciona um ponto no canto esquerdo do dominio, ou seja, V(0, j)
            valores[qtd] = coefXdir;
            colunas[qtd++] = M.matIndicesV[0][j];
        }



        //Conferindo o V(i-1, j)
        if( ESQUERDA_V(i, jU) ) {

            if( M.pontosU[i][jU].tipoBoundary == NOSLIP )
                acrescentaCentro += -coefXesq;
            else if( M.pontosU[i][jU].tipoBoundary == INFLOW )
                acrescentaCentro += -coefXesq;
            else if( M.pontosU[i][jU].tipoBoundary == SLIP )
                acrescentaCentro += -coefXesq;
            else if( M.pontosU[i][jU].tipoBoundary == NEUMANN )
                acrescentaCentro += coefXesq;
            else if( M.pontosU[i][jU].tipoBoundary == SIMETRIA )
                acrescentaCentro += coefXesq;
            else if( M.pontosU[i][jU].tipoBoundary == PERIODICIDADE ) {
                // Nao faz nada, pois vou fazer la pra baixo
            }
            else if( M.pontosU[i][jU].tipoBoundary == EIXO_AXISSIMETRIA )
                acrescentaCentro += coefXesq;
            else
                DeuProblema("\n\n MatrizV: Problema Esquerda\n\n");
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
        }
        else {
            valores[qtd] = coefXesq;
            colunas[qtd++] = M.matIndicesV[i-1][j];
        }

        //Verificando se o V(i, j-1) (ponto de baixo)
        if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].extremo==BAIXO ) {
            DeuProblema("MatrizV: NAO ERA PRA ENTRAR AQUI. \n");
        }
        else {
            valores[qtd] = coefYbaixo;
            colunas[qtd++] = M.matIndicesV[i][j-1];
        }

        // O V(i, j) sempre entra na matriz
        centro = qtd++;
        valores[centro] += 1.0 + coefXcentro + coefYcentro;
        colunas[centro] = M.matIndicesV[i][j];


        //Verificando o V(i, j+1)
        if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].extremo==CIMA ) { //Colocando contorno no U la em cima do dominio
            DeuProblema("MatrizV: NAO ERA PRA ENTRAR AQUI. \n");
        }
        else {
            valores[qtd] = coefYcima;
            colunas[qtd++] = M.matIndicesV[i][j+1];
        }

        //Vericicando o V(i+1, j)
        if( DIREITA_V(i, jU) ) {
            if( M.pontosU[i+1][jU].tipoBoundary==NEUMANN )
                acrescentaCentro += coefXdir; //Neumann
            else if( M.pontosU[i+1][jU].tipoBoundary==SIMETRIA )
                acrescentaCentro += coefXdir;
            else if( M.pontosU[i+1][jU].tipoBoundary==NOSLIP )
                acrescentaCentro += -coefXdir;
            else if( M.pontosU[i+1][jU].tipoBoundary==SLIP )
                acrescentaCentro += -coefXdir;
            else if( M.pontosU[i+1][jU].tipoBoundary==INFLOW )
                acrescentaCentro += -coefXdir;
            else if( M.pontosU[i+1][jU].tipoBoundary==PERIODICIDADE ) {
                //Nao faz nada, pois ja fiz la pra cima
            }
            else
                DeuProblema("\n\n MatrizV: Problema direita\n\n");
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
        }
        else {
            valores[qtd] = coefXdir;
            colunas[qtd++] = M.matIndicesV[i+1][j];
        }

        /// === Finalizando eu vou verificar se estamos em um contorno de PERIODICIDADE na esquerda
        /// === Preciso fazer isso antes de construir os outros elementos devido a estrutura do petsc...
        if( ESQUERDA_V(i, jU) && M.pontosU[i+1][jU].tipoBoundary==PERIODICIDADE ) {
            // Adiciona um ponto no canto direito do dominio, ou seja, V(Nx-1, j)
            valores[qtd] = coefXesq;
            colunas[qtd++] = M.matIndicesV[M.Nx-1][j];
        }

        valores[centro] += acrescentaCentro;
        ierr = MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode VetorIndependenteV(Vec VetVetor, 
                                    double **UVelho, double **VVelho, double **UNovo, double **VNovo,
                                    double **PVelho, double **Txy, double **Tyy, 
                                    double **ViscosityFunction,
                                    MALHA M)
{
    int linha = 0;
    int i, j;
    double theta = 1.0, valor, g;
    double hx1, hx2;
    int iStart, iEnd;
    PetscErrorCode ierr;

    ierr = VecGetOwnershipRange(VetVetor, &iStart, &iEnd); CHKERRQ(ierr);

    g = M.gravY; //Gravidade adimensional (geralmente = -1)

    for( linha=iStart; linha<iEnd; linha++ ) {
        //PrintDebug("LINHA: %d\n", linha);
        i = M.indicesV[linha].i;
        j = M.indicesV[linha].j;

        if( M.pontosV[i][j].tipo==EMPTY ) {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("\n\n VetorIndependenteV: PROBLEMA %d %d\n\n", i, j);
        }
        else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==INFLOW ) {
            valor = VVelho[i][j];
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==NOSLIP ) {
            valor = VVelho[i][j];
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==SLIP ) {
            valor = VVelho[i][j];
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosV[i][j].tipo==SURFACE ) {
            valor = VNovo[i][j];
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==NEUMANN ) {
            valor = 0.0;
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==SIMETRIA ) {
            valor = 0.0;
            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
        /// Caso FULL
        else {
            /// If we are using a generalized newtonian model then the "div(T)" term will be a bit different from the viscoelastic case
            // double div_stress = (M.tipo_modelo==MODELO_VP) ? DivT_v_NewtGen(M, UVelho, VVelho, ViscosityFunction, i, j) : DivT_v(Txy, Tyy, i, j, M);
            double div_stress = DivT_v(Txy, Tyy, i, j, M);

            valor = VVelho[i][j] - M.dt*( Cy(UVelho, VVelho, i, j, M) + GradPy(PVelho, i, j, M) - (1.0/M.Re)*div_stress - (((1.0-theta)/M.Re)*Ly(VVelho, i, j, M)) - (g/(M.Fr*M.Fr)) );


            /// === Em fase de testes: colocando tensao superficial em algumas celulas FULL
            // M.tensaoSup_full = 1;
            // if( M.tensaoSup_full ) {
            //     if( (j!=M.Ny && M.celulasParede[i][j]==SURFACE) || (j!=0 && M.celulasParede[i][j-1]==SURFACE) ) {
            //         double normalX, normalY, curvatura;
            //         CalculaNormalCurvaturaEquacaoMovimento(&M, 0.5*(M.x[i] + M.x[i+1]), M.y[j], &normalX, &normalY, &curvatura);
            //         valor += M.dt*(1.0/M.weber)*curvatura*normalY;
            //         // PrintDebug("VALOR Y: %lf %lf\n", M.dt*(1.0/M.weber)*curvatura*normalY, curvatura);
            //     }
            // }


            //Condicao de contorno NOSLIP quando a parede vertical a direita esta se movimentando
            int jU = (j==M.Ny) ? (j-1) : j;
            if( DIREITA_V(i, jU) ) {

                /// Determinando o tipo deste contorno
                /// Precisa olhar dois pontos pra nao bugar casos de quina tipo na contracao.... (gambiarra)
                int tipoCima = M.pontosU[i+1][jU].tipoBoundary;
                int tipoBaixo = (jU==0) ? NAO_BOUNDARY : M.pontosU[i+1][jU-1].tipoBoundary;
                int movimentoCima = M.pontosU[i+1][jU].movimento;
                int movimentoBaixo = (jU==0) ? SEM_MOVIMENTO : M.pontosU[i+1][jU-1].movimento;
                int tipo = NAO_BOUNDARY, movimento = SEM_MOVIMENTO;

                if( tipoCima==NOSLIP || tipoCima==SLIP ) {
                    tipo = tipoCima;
                    movimento = movimentoCima;
                }
                if( tipoBaixo==NOSLIP || tipoBaixo==SLIP ) {
                    tipo = tipoBaixo;
                    movimento = movimentoBaixo;
                }

                if( (tipo==NOSLIP || tipo==SLIP) && (movimento==MOVIMENTO_Y) ) {
                    hx1 = (i==0) ? M.dx[i] : 0.5*(M.dx[i-1] + M.dx[i]);
                    hx2 = (i==M.Nx-1) ? M.dx[i] : 0.5*(M.dx[i] + M.dx[i+1]);

                    double coefXdir = -M.beta*theta*(M.dt/M.Re)*( 2.0/(hx2*(hx1+hx2)) );
                    valor += -coefXdir*2.0*(FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, VVelho ));
                }
            }

            //Condicao de contorno NOSLIP quando a parede vertical a esquerda esta se movimentando
            if( ESQUERDA_V(i, jU) ) {

                /// Determinando o tipo deste contorno
                /// Precisa olhar dois pontos pra nao bugar casos de quina tipo na contracao.... (gambiarra)
                int tipoCima = M.pontosU[i][jU].tipoBoundary;
                int tipoBaixo = (jU==0) ? NAO_BOUNDARY : M.pontosU[i][jU-1].tipoBoundary;
                int movimentoCima = M.pontosU[i][jU].movimento;
                int movimentoBaixo = (jU==0) ? SEM_MOVIMENTO : M.pontosU[i][jU-1].movimento;
                int tipo = NAO_BOUNDARY, movimento = SEM_MOVIMENTO;

                if( tipoCima==NOSLIP || tipoCima==SLIP ) {
                    tipo = tipoCima;
                    movimento = movimentoCima;
                }
                if( tipoBaixo==NOSLIP || tipoBaixo==SLIP ) {
                    tipo = tipoBaixo;
                    movimento = movimentoBaixo;
                }

                if( (tipo==NOSLIP || tipo==SLIP) && (movimento==MOVIMENTO_Y) ) {
                    hx1 = (i==0) ? M.dx[i] : 0.5*(M.dx[i-1] + M.dx[i]);
                    hx2 = (i==M.Nx-1) ? M.dx[i] : 0.5*(M.dx[i] + M.dx[i+1]);

                    double coefXesq = -M.beta*theta*(M.dt/M.Re)*( 2.0/(hx1*(hx1+hx2)) );
                    valor += -coefXesq*2.0*(FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, VVelho ));
                }

            }


            VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES);
        }
    }

    VecAssemblyBegin(VetVetor);
    VecAssemblyEnd(VetVetor);
    return 0;
}

PetscErrorCode InicializaMatrizCoeficientesW(Mat A, int NumLinhas, MALHA M, int iStart, int iEnd)
{
    int linha, indiceVetor = 0, coluna;
    int i, j;
    int d_nnz[iEnd - iStart], o_nnz[iEnd - iStart], qtd[iEnd - iStart];
    PetscErrorCode ierr;

    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesW[linha].i;
        j = M.indicesW[linha].j;


        o_nnz[indiceVetor] = 0;

        //Diagonal principal sempre vai ter
        qtd[indiceVetor] = 1;
        d_nnz[indiceVetor] = 1;


        //Verificando o da esquerda
        if( i!=0 && M.pontosW[i-1][j].tipo==FULL ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesW[i-1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }

        //Verificando o de baixo
        if( j!=0 && M.pontosW[i][j-1].tipo==FULL ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesW[i][j-1];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }

        //Verificando o da direita
        if( i!=M.Nx-1 && M.pontosW[i+1][j].tipo==FULL ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesW[i+1][j];
            if( coluna>=iStart && coluna<iEnd )
                d_nnz[indiceVetor]++;
            else
                o_nnz[indiceVetor]++;
        }

        //Verificando o de cima
        if( j!=M.Ny-1 && M.pontosW[i][j+1].tipo==FULL ) {
            qtd[indiceVetor]++;
            coluna = M.matIndicesW[i][j+1];
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

PetscErrorCode MontaMatrizCoeficientesW(Mat A, MALHA M)
{
    int i, j;
    int qtd;
    int linha=0, colunas[5], centro;
    double valores[5], acrescentaCentro;
    double rImenos1, rImais1, rI;
    double coefXesq, coefXdir, coefXcentro, coefYbaixo, coefYcima, coefYcentro;
    double hx1, hx2, hy1, hy2;
    double theta = 1.0, beta;
    int iStart, iEnd;

    beta = M.beta;
    theta = 1.0;

    MatGetOwnershipRange(A, &iStart, &iEnd);

    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesW[linha].i;
        j = M.indicesW[linha].j;

        rImais1 = (i!=M.Nx-1) ? M.r[i+1] : (M.r[i] + M.dx[i]);
        rImenos1 = (i!=0) ? M.r[i-1] : (M.r[i] - M.dx[i]);
        rI = M.r[i];


        hx1 = (i==0) ? M.dx[i] : 0.5*( M.dx[i-1] + M.dx[i] );
        hx2 = (i==M.Nx-1) ? M.dx[i] : 0.5*( M.dx[i] + M.dx[i+1] );
        hy1 = (j==0) ? M.dy[j] : 0.5*(M.dy[j-1] + M.dy[j]);
        hy2 = (j==M.Ny-1) ? M.dy[j] : 0.5*(M.dy[j] + M.dy[j+1]);

        coefXesq = -beta*theta*(M.dt/M.Re)*( (2.0 + (hx2/rI))*rImenos1/(rI*hx1*(hx1+hx2)) );
        coefXdir = -beta*theta*(M.dt/M.Re)*( (2.0 - (hx1/rI))*rImais1/(rI*hx2*(hx1+hx2)) );
        coefXcentro = -beta*theta*(M.dt/M.Re)*( -(2.0 + ((hx2-hx1)/rI))/(hx1*hx2) );
        coefYbaixo = -beta*theta*(M.dt/M.Re)*( 2.0/(hy1*(hy1+hy2)) );
        coefYcima = -beta*theta*(M.dt/M.Re)*( 2.0/(hy2*(hy1+hy2)) );
        coefYcentro = -beta*theta*(M.dt/M.Re)*( -2.0/(hy1*hy2) );


        qtd = 0;
        valores[0]=valores[1]=valores[2]=valores[3]=valores[4]=0.0;
        acrescentaCentro = 0.0;

        //Esquerda: Colocando os coeficientes do W(i-1, j)
        if( i!=0 && M.pontosW[i-1][j].tipo==FULL )
        {
            valores[qtd] = coefXesq;
            colunas[qtd++] = M.matIndicesW[i-1][j];
        }
        else {
            if( M.pontosU[i][j].tipoBoundary == NOSLIP )
                acrescentaCentro += -coefXesq;
            else if( M.pontosU[i][j].tipoBoundary == INFLOW )
                acrescentaCentro += -coefXesq;
            else if( M.pontosU[i][j].tipoBoundary == EIXO_AXISSIMETRIA )
                acrescentaCentro += coefXesq;
            else if( M.pontosU[i][j].tipoBoundary == NEUMANN )
                acrescentaCentro += coefXesq;
            else {
                DeuProblema("\n\n MatrizCoeficientesW: Problema Esquerda\n\n");
            }
        }

        //Baixo: Colocando o coeficiente do W(i, j-1)
        if( j!=0 && M.pontosW[i][j-1].tipo==FULL ) {
            valores[qtd] = coefYbaixo;
            colunas[qtd++] = M.matIndicesW[i][j-1];
        }
        else {
            if( M.pontosV[i][j].tipoBoundary == NOSLIP )
                acrescentaCentro += -coefYbaixo;
            else if( M.pontosV[i][j].tipoBoundary == INFLOW )
                acrescentaCentro += -coefYbaixo;
            else if( M.pontosV[i][j].tipoBoundary == NEUMANN )
                acrescentaCentro += coefYbaixo;
            else {
                DeuProblema("\n\n MatrizCoeficientesW: Problema Baixo\n\n");
            }
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
        }

        //Centro: Colocando o coeficiente do W(i, j)
        centro = qtd;
        //valores[qtd] += 1 + theta*rx*M.r[i]*((1.0/rImaisMeio) + (1.0/rImenosMeio)) + 2*theta*ry; //Sem regra da cadeia
        valores[qtd] += 1.0 + coefXcentro + coefYcentro;
        colunas[qtd++] = M.matIndicesW[i][j];

        //Cima: Colocando o coeficiente do W(i, j+1)
        if( j!=M.Ny-1 && M.pontosW[i][j+1].tipo==FULL ) {
            valores[qtd] = coefYcima;
            colunas[qtd++] = M.matIndicesW[i][j+1];
        }
        else {
            if( M.pontosV[i][j+1].tipoBoundary == NOSLIP )
                acrescentaCentro += -coefYcima;
            else if( M.pontosV[i][j+1].tipoBoundary == INFLOW )
                acrescentaCentro += -coefYcima;
            else if( M.pontosV[i][j+1].tipoBoundary == NEUMANN )
                acrescentaCentro += coefYcima;
            else {
                DeuProblema("\n\n MatrizCoeficientesW:  problema Cima\n");
            }
        }

        //Direita: Colocando o coeficiente do W(i+1, j)
        if( i!=M.Nx-1 && M.pontosW[i+1][j].tipo==FULL ) {
            valores[qtd] = coefXdir;
            colunas[qtd++] = M.matIndicesW[i+1][j];
        }
        else {
            if( M.pontosU[i+1][j].tipoBoundary == NOSLIP )
                acrescentaCentro += -coefXdir;
            else if( M.pontosU[i+1][j].tipoBoundary == INFLOW )
                acrescentaCentro += -coefXdir;
            else if( M.pontosU[i+1][j].tipoBoundary == NEUMANN )
                acrescentaCentro += coefXdir;
            else
                DeuProblema("\n\n MatrizCoeficientesW: problema direita\n");
            //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
        }

        valores[centro] += acrescentaCentro;
        MatSetValues(A, 1, &linha, qtd, colunas, valores, INSERT_VALUES);

//        if( j==0 || j==M.Ny-1 )
//            PrintDebug("vv (%2d, %2d): %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf\n", i, j, valores[0], valores[1], valores[2], valores[3], valores[4]);
    }

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    return 0;
}

PetscErrorCode VetorIndependenteW(Vec VetVetor, double **UVelho, double **VVelho, double **WVelho, double **PVelho,
                                  double **Txt, double **Tyt, MALHA M)
{
    int linha=0;
    int i, j;
    double valor, theta = 1.0;
    double rImais1, rImenos1, rI;
    double hx1, hx2;
    int iStart, iEnd;
    PetscErrorCode ierr;

    ierr = VecGetOwnershipRange(VetVetor, &iStart, &iEnd); CHKERRQ(ierr);


    for( linha=iStart; linha<iEnd; linha++ ) {
        i = M.indicesW[linha].i;
        j = M.indicesW[linha].j;

        if( M.pontosW[i][j].tipo==EMPTY )
            printf("\n PROBLEMA \n");
        else {
            valor = WVelho[i][j] - M.dt*( Ct(UVelho, VVelho, WVelho, i, j, M) - (1.0/M.Re)*DivT_t(Txt, Tyt, i, j, M)  );

            //Condicao de contorno quando a parede NOSLIP esta se movendo na direcao THETA
            if( ((i==M.Nx-1) || (M.pontosW[i+1][j].tipo!=FULL)) && (M.pontosU[i+1][j].tipoBoundary==NOSLIP) ) {
                rI = M.r[i];
                rImais1 = (i!=M.Nx-1) ? M.r[i+1] : (M.r[i] + M.dx[i]);
                hx1 = (i==0) ? M.dx[i] : 0.5*( M.dx[i-1] + M.dx[i] );
                hx2 = (i==M.Nx-1) ? M.dx[i] : 0.5*( M.dx[i] + M.dx[i+1] );
                double coefXdir = -M.beta*theta*(M.dt/M.Re)*( (2.0 - (hx1/rI))*rImais1/(rI*hx2*(hx1+hx2)) );

                valor += - 2.0*coefXdir*(M.pontosU[i+1][j].valorDirichlet);
            }
            else if( ((i==0) || (M.pontosW[i-1][j].tipo!=FULL)) && (M.pontosU[i][j].tipoBoundary==NOSLIP) ) {
                rI = M.r[i];
                rImenos1 = (i!=0) ? M.r[i-1] : (M.r[i] - M.dx[i]);
                hx1 = (i==0) ? M.dx[i] : 0.5*( M.dx[i-1] + M.dx[i] );
                hx2 = (i==M.Nx-1) ? M.dx[i] : 0.5*( M.dx[i] + M.dx[i+1] );
                double coefXesq = -M.beta*theta*(M.dt/M.Re)*( (2.0 + (hx2/rI))*rImenos1/(rI*hx1*(hx1+hx2)) );

                valor += - 2.0*coefXesq*(M.pontosU[i][j].valorDirichlet);
            }

//            if( j==0 || j==M.Ny-1 )
//                PrintDebug("(%2d, %2d): %20.15lf %20.15lf %20.15lf %20.15lf\n", i, j, Txt[i][j], Tyt[i][j], WVelho[i][j], valor);


            ierr = VecSetValues(VetVetor, 1, &linha, &valor, INSERT_VALUES); CHKERRQ(ierr);
        }
    }


    ierr = VecAssemblyBegin(VetVetor); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(VetVetor); CHKERRQ(ierr);
    return 0;
}


void CalculaNormalCurvaturaEquacaoMovimento(MALHA *M, double X, double Y, double *ResultNx, double *ResultNy, double *ResultCurv)
{
    PONTO *p, *pEscolhido;
    CURVA *c, *cEscolhida;
    double minDistancia, distancia, xi, eta;
    int i, j;

    int qtdParticulas = 15;

    /// Encontrando qual ponto da interface esta mais proxima de (x, Y)
    minDistancia = 1e+10;
    pEscolhido = NULL;
    cEscolhida = NULL;
    for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
//        c->medioX = 0.0;
//        c->medioY = 0.0;
//        c->qtdPontosCurva = 0;

        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            distancia = ((X - p->x)*(X - p->x)) + ((Y - p->y)*(Y - p->y));

            if( distancia<minDistancia ) {
                pEscolhido = p;
                cEscolhida = c;
                minDistancia = distancia;
            }

//            c->medioX += p->x;
//            c->medioY += p->y;
//            c->qtdPontosCurva++;
        }

//        c->medioX /= c->qtdPontosCurva;
//        c->medioY /= c->qtdPontosCurva;
    }


    /// == Verificando em qual celula esta este ponto
    /// == Depois encontrando um vetor normal ruim estilo freeflow apenas pra saber qual
    /// === qual eh o sentido apontando pra "fora" do bloco de fluido
    /// === Sao oito casos possiveis de vetor normal...
//    double direcaoX, direcaoY;
//    EncontraCelulaDaParticula(M, pEscolhido->x, pEscolhido->y, &i, &j);
//    CalculaDirecaoExternaSurfaceFull(M, i, j, &direcaoX, &direcaoY);




    static double A[4][4];
    static double det;
    for( i=0; i<4; i++ )
        for( j=0; j<4; j++ )
            A[i][j] = 0.0;

    /// == Construindo o sistema para calcular o vetor normal
    PONTO *pTras = pEscolhido;
    PONTO *pFrente = pEscolhido;
    for( i=0; i<qtdParticulas; i++ ) {
        pTras = pTras->ant;
        pFrente = pFrente->prox;

        p = pTras;
        //lhs
        A[0][0] += (p->x)*(p->x);
        A[0][1] += p->x;
        A[1][0] += p->x;
        A[1][1] += 1.0;
        //rhs
        A[0][2] += (p->x)*(p->y);
        A[1][2] += p->y;
    }

    /// Resolvendo o sistema para calcular o vetor normal
    det = A[0][0]*A[1][1] - A[0][1]*A[1][0]; //Determinante
    if( fabs(det)<1e-9 ) {
        *ResultNx = 1.0;
        *ResultNy = 0.0;
    }
    else {
        double a = (A[0][2]*A[1][1] - A[1][2]*A[0][1])/det;
        double norma = sqrt(a*a + 1.0);
        *ResultNx = -a/norma;
        *ResultNy = 1.0/norma;
    }

    /// Corrigindo o sentido da normal, se necessario
    /// Apenas uma versao inicial, talvez melhorar isso....
    ConfereDirecaoVetorNormal(M, pEscolhido, cEscolhida, ResultNx, ResultNy);
//    if( ( (*ResultNx)*direcaoX + (*ResultNy)*direcaoY ) < 0 ) { //Angulo em modulo eh maior que 90 graus? muda o sentido
//        *ResultNx = - *ResultNx;
//        *ResultNy = - *ResultNy;
//    }

//    PrintDebug("%lf %lf --- %lf %lf\n", X, Y, *ResultNx, *ResultNy);


    for( i=0; i<4; i++ )
        for( j=0; j<4; j++ )
            A[i][j] = 0.0;

    /// ==== Construindo o sistema pra calcular a curvatura
    pTras = pEscolhido;
    pFrente = pEscolhido;
    for( i=0; i<qtdParticulas; i++ ) {
        pTras = pTras->ant;
        pFrente = pFrente->prox;

        //Freeflow
        p = pTras;
        xi = (p->x-X)*(*ResultNy)- (p->y-Y)*(*ResultNx);
        eta = (p->x-X)*(*ResultNx) + (p->y-Y)*(*ResultNy);

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

        //Freeflow
        p = pFrente;
        xi = (p->x-X)*(*ResultNy)- (p->y-Y)*(*ResultNx);
        eta = (p->x-X)*(*ResultNx) + (p->y-Y)*(*ResultNy);

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
    }

    /// == Resolvendo o sistema por cramer pra calcular a curvatura
    det = A[0][0]*A[1][1]*A[2][2] + A[2][0]*A[0][1]*A[1][2] + A[1][0]*A[2][1]*A[0][2]
            - A[2][0]*A[1][1]*A[0][2] - A[0][0]*A[2][1]*A[1][2] - A[2][2]*A[1][0]*A[0][1]; //Determinante

	if( fabs(det)<1e-9 ) {
		*ResultCurv = 0.0;
	}

    //Usando regra de cramer
	double d = A[0][3]*A[1][1]*A[2][2] + A[2][3]*A[0][1]*A[1][2] + A[1][3]*A[2][1]*A[0][2]
            - A[2][3]*A[1][1]*A[0][2] - A[0][3]*A[2][1]*A[1][2] - A[2][2]*A[1][3]*A[0][1]; //Determinante
    d = d/det;

    *ResultCurv = -2.0*d;


//    PrintDebug("%lf %lf --- %lf\n", X, Y, *ResultCurv);

    return;
}

void ConfereDirecaoVetorNormal(MALHA *M, PONTO *P, CURVA *C, double *NormalX, double *NormalY)
{
    PONTO *p1, *p2;
    double h = 0.5*M->minDxDy;
    int resultado;

    double x= P->x + h*(*NormalX);
    double y = P->y + h*(*NormalY);

    //Verificando se p2 esta dentro ou fora da interface. even-odd rule
    p2 = C->pontoInicial;
    resultado = 0;
    LOOP_CURVA(C) {
        p1 = loop_stuff.ponto;

        if( ((p1->y > y) != (p2->y > y)) && (x < (p1->x + (p2->x - p1->x) * (y - p1->y) / (p2->y - p1->y)) ) ) {
            resultado = !resultado;
        }
        p2 = p1;
    }

    /// == Se a normal estiver apontando pra dentro, muda o sentido dela
    if( resultado ) {
        *NormalX = - *NormalX;
        *NormalY = - *NormalY;
    }

    return;
}

//void CalculaDirecaoExternaSurfaceFull(MALHA *M, int i, int j, double *DirecaoX, double *DirecaoY)
//{
//    if( M->celulas[i][j]==FULL ) {
//
//        int criterioEsq = (M->celulas[i-1][j]==FULL) && (M->celulasParede[i-1][j]!=SURFACE);
//        int criterioDir = (M->celulas[i+1][j]==FULL) && (M->celulasParede[i+1][j]!=SURFACE);
//        int criterioBaixo = (M->celulas[i][j-1]==FULL) && (M->celulasParede[i][j-1]!=SURFACE);
//        int criterioCima = (M->celulas[i][j+1]==FULL) && (M->celulasParede[i][j+1]!=SURFACE);
//
//        if( criterioEsq && !criterioDir && !criterioBaixo && !criterioCima ) {
//            *DirecaoX = 1.0;
//            *DirecaoY = 0.0;
//        }
//        else if( !criterioEsq && criterioDir && !criterioBaixo && !criterioCima ) {
//            *DirecaoX = -1.0;
//            *DirecaoY = 0.0;
//        }
//        else if( !criterioEsq && !criterioDir && criterioBaixo && !criterioCima ) {
//            *DirecaoX = 0.0;
//            *DirecaoY = 1.0;
//        }
//        else if( !criterioEsq && !criterioDir && !criterioBaixo && criterioCima ) {
//            *DirecaoX = 0.0;
//            *DirecaoY = -1.0;
//        }
//        else if( criterioEsq && !criterioDir && criterioBaixo && !criterioCima ) {
//            *DirecaoX = sqrt(2.0)/2.0;
//            *DirecaoY = sqrt(2.0)/2.0;
//        }
//        else if( !criterioEsq && criterioDir && criterioBaixo && !criterioCima ) {
//            *DirecaoX = -sqrt(2.0)/2.0;
//            *DirecaoY = sqrt(2.0)/2.0;
//        }
//        else if( criterioEsq && !criterioDir && !criterioBaixo && criterioCima ) {
//            *DirecaoX = sqrt(2.0)/2.0;
//            *DirecaoY = -sqrt(2.0)/2.0;
//        }
//        else if( !criterioEsq && criterioDir && !criterioBaixo && criterioCima ) {
//            *DirecaoX = -sqrt(2.0)/2.0;
//            *DirecaoY = -sqrt(2.0)/2.0;
//        }
//        // Caso especial quando tem FULL_SURFACE em tres cantos
//        else if( M->celulas[i][j+1]==SURFACE && M->celulas[i][j-1]==FULL && M->celulas[i-1][j]==FULL && M->celulas[i+1][j]==FULL ) {
//            *DirecaoX = 0.0;
//            *DirecaoY = 1.0;
//        }
//        else if( M->celulas[i][j+1]==FULL && M->celulas[i][j-1]==SURFACE && M->celulas[i-1][j]==FULL && M->celulas[i+1][j]==FULL ) {
//            *DirecaoX = 0.0;
//            *DirecaoY = -1.0;
//        }
//        else {
//            DesenhaMalhaVTK(*M, 0);
//            ImprimeInterfaceVTK(*M, 100000);
//            DeuProblema("1. CASO INESPERADO DO VETOR NORMAL NA FULL_SURFACE %d %d\n\n", i, j);
//        }
//    }
//    else if( M->celulas[i][j]==SURFACE ) {
//        //Superficie vertical: empty na direita
//        if( LADO_VERTICAL_DIREITA(i, j) ) {
//            *DirecaoX = 1.0;
//            *DirecaoY = 0.0;
//        }
//        //Superficie vertical: empty na esquerda
//        else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
//            *DirecaoX = -1.0;
//            *DirecaoY = 0.0;
//        }
//        //Superficie horizontal: empty em cima
//        else if( LADO_HORIZONTAL_CIMA(i, j) ) {
//            *DirecaoX = 0.0;
//            *DirecaoY = 1.0;
//        }
//        //Superficie horizontal: empty embaixo
//        else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
//            *DirecaoX = 0.0;
//            *DirecaoY = -1.0;
//        }
//        //Superficie quina: empty em cima e na direita
//        else if( QUINA_CIMA_DIREITA(i, j) ) {
//            *DirecaoX = sqrt(2.0)/2.0;
//            *DirecaoY = sqrt(2.0)/2.0;
//        }
//        //Superficie quina: empty em cima e na esquerda
//        else if( QUINA_CIMA_ESQUERDA(i, j) ) {
//            *DirecaoX = -sqrt(2.0)/2.0;
//            *DirecaoY = sqrt(2.0)/2.0;
//        }
//        //Superficie quina: empty em embaixo e na direita
//        else if( QUINA_BAIXO_DIREITA(i, j) ) {
//            *DirecaoX = sqrt(2.0)/2.0;
//            *DirecaoY = -sqrt(2.0)/2.0;
//        }
//        //Superficie quina: empty em embaixo e na esquerda
//        else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
//            *DirecaoX = -sqrt(2.0)/2.0;
//            *DirecaoY = -sqrt(2.0)/2.0;
//        }
//        else {
//            DesenhaMalhaVTK(*M, 0);
//            DeuProblema("2. CASO INESPERADO DO VETOR NORMAL NA FULL_SURFACE %d %d\n\n", i, j);
//        }
//    }
//    else {
//        DesenhaMalhaVTK(*M, 0);
//        DeuProblema("3. CASO INESPERADO DO VETOR NORMAL NA FULL_SURFACE %d %d\n\n", i, j);
//    }
//
//    return;
//}

void CalculaResiduoMomentoU(MALHA *M, double **UVelho, double **UMeio, double **UNovo,
                            double **VVelho, double **VMeio, double **VNovo,
                            double **PVelho, double **PNovo,
                            double **TxxVelho, double **TxxNovo, double **TxyVelho, double **TxyNovo, double **TttVelho, double **TttNovo,
                            double **Residuo)
{
    double parteTempo, parteConv, partePressao, parteLaplac, parteDivT;
    int i, j;

    for( i=0; i<=M->Nx; i++) {
        for( j=0; j<M->Ny; j++ ) {

            //Vou calcular apenas em pontos FULL
            if( M->pontosU[i][j].tipo!=FULL )
                continue;

            parteTempo = (UNovo[i][j] - UVelho[i][j])/(M->dt);
            parteConv = Cx(UNovo, VNovo, NULL, i, j, 0.0, *M);
            partePressao = GradPx(PNovo, i, j, *M);
            parteLaplac = -(M->beta/M->Re)*Lx(UNovo, i, j, 0.0, *M);
            parteDivT = -(1.0/M->Re)*DivT_u(TxxNovo, TxyNovo, TttNovo, i, j, *M);

            /// Apenas testando no tempo N pra ver se da zero
//            parteTempo = (UMeio[i][j] - UVelho[i][j])/(M->dt);
//            parteConv = Cx(UVelho, VVelho, NULL, i, j, 0.0, *M);
//            partePressao = GradPx(PVelho, i, j, *M);
//            parteLaplac = -(M->beta/M->Re)*Lx(UMeio, i, j, *M);
//            parteDivT = -(1.0/M->Re)*DivT_u(TxxVelho, TxyVelho, i, j, *M);

            Residuo[i][j] = parteTempo + parteConv + partePressao + parteLaplac + parteDivT;
        }
    }

    return;
}



