#include "MetodoProjecao.h"

extern MALHA malha_temp;

extern double **debugsCell;

double FuncaoParabolica(double PontoMin, double PontoMax, double ValorMax, double Ponto, double t, MALHA Malha)
{
    ///Caso comum
    double centro = 0.5*(PontoMin + PontoMax);
    return ValorMax*( ( (Ponto-PontoMin)*(Ponto-PontoMax) )/( (centro-PontoMin)*(centro-PontoMax) ) );
    
    /// Case: yield-stress Bingham model
    // double y = (2*Ponto - PontoMax - PontoMin)/(PontoMax - PontoMin);
    // double Bi = Malha.Bi;
    // double Re = Malha.Re;
    // double P_grad = 10;
    // double valor_inflow;
    // if( y<-Bi/(Re*P_grad) )
	// 	valor_inflow = 0.5*Re*P_grad*(1.0 - y*y) - Bi*(y + 1.0);
	// else if( y<Bi/(Re*P_grad) )
	// 	valor_inflow = 0.5*Re*P_grad - Bi + 0.5*Bi*Bi/(Re*P_grad);
	// else
	//     valor_inflow = 0.5*Re*P_grad*(1.0 - y*y) + Bi*(y - 1.0);
    
    // // PrintDebug("inflow: %lf %lf\n", y, valor_inflow);
    // return valor_inflow;
}

void CondicaoInicialProjecao(SOLVER_PROJECAO Solver, MALHA M, double T)
{
    int i, j;

//    double escala_random = 1e-3;

    for(i=0; i<=M.Nx; i++) {
        for(j=0; j<M.Ny; j++) {
            if( M.pontosU[i][j].tipo==FULL  ) {
                Solver.UVelho[i][j] = Solver.UMeio[i][j] = Solver.UNovo[i][j] = M.velInicialX;
//
//                double random_value = (double)rand()/RAND_MAX*2.0-1.0;
//                random_value *= escala_random;
//                Solver.UVelho[i][j] = Solver.UMeio[i][j] = Solver.UNovo[i][j] = random_value;
//                PrintDebug("rand %lf\n", random_value);

                // Verifica se esta dentro de alguma regiao free surface com velocidade (tipo uma gota)
                BLOCO *bloco = BuscaBlocoFreeSurface(&M, M.x[i], 0.5*(M.y[j]+M.y[j+1]));

                if( bloco ) {
                    Solver.UVelho[i][j] = Solver.UMeio[i][j] = Solver.UNovo[i][j] = bloco->freeVelX;
                }
            }

            ///APENAS PARA O CASO ENCONTRO DE GOTAS
//            if( M.pontosU[i][j].tipo==FULL && M.x[i]<1.5 )
//                Solver.UVelho[i][j] = Solver.UMeio[i][j] = Solver.UNovo[i][j] = 2.0;
//            else if( M.pontosU[i][j].tipo==FULL && M.x[i]>1.5 )
//                Solver.UVelho[i][j] = Solver.UMeio[i][j] = Solver.UNovo[i][j] = -2.0;
        }
    }

    for(i=0; i<M.Nx; i++) {
        for(j=0; j<=M.Ny; j++) {
            Solver.VVelho[i][j] = Solver.VMeio[i][j] = Solver.VNovo[i][j] = 0.0;
            if( M.pontosV[i][j].tipo==FULL ) {
                Solver.VVelho[i][j] = Solver.VMeio[i][j] = Solver.VNovo[i][j] = M.velInicialY;

//                double random_value = (double)rand()/RAND_MAX*2.0-1.0;
//                random_value *= escala_random;
//                Solver.VVelho[i][j] = Solver.VMeio[i][j] = Solver.VNovo[i][j] = random_value;

                // Verifica se esta dentro de alguma regiao free surface com velocidade (tipo uma gota)
                BLOCO *bloco = BuscaBlocoFreeSurface(&M, 0.5*(M.x[i]+M.x[i+1]), M.y[j]);

                if( bloco ) {
                    Solver.VVelho[i][j] = Solver.VMeio[i][j] = Solver.VNovo[i][j] = bloco->freeVelY;
                }
            }
        }
    }

    for(i=0; i<M.Nx; i++) {
        for(j=0; j<M.Ny; j++) {
            Solver.PVelho[i][j] = Solver.Psi[i][j] = Solver.PNovo[i][j] = 0.0;
            Solver.TxxVelho[i][j] = Solver.TxyVelho[i][j] = Solver.TyyVelho[i][j] = 0.0;
            Solver.TxxNovo[i][j] = Solver.TxyNovo[i][j] = Solver.TyyNovo[i][j] = 0.0;

            Solver.lambdaVelho[i][j] = Solver.lambdaNovo[i][j] = 1.0;
            Solver.nuVelho[i][j] = Solver.nuNovo[i][j] = 1.0;
            Solver.muVelho[i][j] = Solver.muNovo[i][j] = 0.0;

            if( M.pontosP[i][j].tipo!=EMPTY )
                Solver.evptStructureVelho[i][j] = Solver.evptStructureNovo[i][j] = 1.0;
        

            /// CANAL: Colocando pressao inicial apenas no caso do canal poisseuile
            if( M.pontosP[i][j].tipo!=EMPTY ) {
                // double x = 0.5*( M.x[i] + M.x[i+1] );//, y = 0.5*(M.y[j] + M.y[j+1]);
                // double p_in = 100.0, p_out = 0;
                // Solver.PVelho[i][j] = Solver.PNovo[i][j] = p_in - (x + 5.0)*(p_in - p_out)/10.0;

                // Solver.TxyVelho[i][j] = Solver.TxyNovo[i][j] = y*(p_in - p_out)/10.0;
                // Solver.PVelho[i][j] = Solver.PNovo[i][j] = p_in - (x + 0.5)*(p_in - p_out)/1.0;
                // Solver.TxyVelho[i][j] = Solver.TxyNovo[i][j] = y*(p_in - p_out)/1.0;
            }
        }
    }

    if( M.tipoCoord==AXI_CILINDRICO ) {
        for(i=0; i<M.Nx; i++)
            for(j=0; j<M.Ny; j++)
                Solver.WVelho[i][j] = Solver.WMeio[i][j] = Solver.WNovo[i][j] = 0.0;

        /// === Jogando um residuo inicial em T
//        double residuo = 1e-5;
////        double residuo = 0.0;
//        for(i=0; i<M.Nx; i++) {
//            for(j=0; j<M.Ny; j++) {
//                Solver.TxxVelho[i][j] = Solver.TxxNovo[i][j] = 10.0*residuo;
//                Solver.TyyVelho[i][j] = Solver.TyyNovo[i][j] = 0.5*residuo;
//                Solver.TttVelho[i][j] = Solver.TttNovo[i][j] = 0.142*residuo;
//            }
//        }
    }

    CondicaoContornoProjecao(Solver.UVelho, Solver.VVelho, Solver.WVelho, M, T);
    CondicaoContornoProjecao(Solver.UNovo, Solver.VNovo, Solver.WVelho, M, T);
    return;
}

void CondicaoContornoProjecao(double **U, double **V, double **W, MALHA M, double T)
{
    int i, j;
    TIPO_BOUNDARY contorno;
    BLOCO *bloco;
    double pMin, pMax;

	for(bloco=M.primBloco; bloco!=NULL; bloco=bloco->prox) {

        if( bloco->tipoU==BOUNDARY ) {
            contorno = bloco->tipoBoundaryU;
            if( (contorno==INFLOW) && (bloco->iMin==bloco->iMax) ) {
                pMin = pMax = M.yMin;
                for( j=0; j<bloco->jMin; j++ )
                    pMin += M.dy[j];
                for( j=0; j<bloco->jMax+1; j++ )
                    pMax += M.dy[j];

                // Verificando se tem condicao de simetria embaixo
                i = (bloco->iMin==M.Nx) ? bloco->iMin-1 : bloco->iMin;
                j = bloco->jMin;
                if( M.pontosV[i][j].tipoBoundary==SIMETRIA )
                    pMin = pMin - (pMax - pMin);

                // Verificando se tem condicao de simetria acima
                i = (bloco->iMin==M.Nx) ? bloco->iMin-1 : bloco->iMin;
                j = bloco->jMax;
                if( M.pontosV[i][j+1].tipoBoundary==SIMETRIA )
                    pMax = pMax + (pMax - pMin);

                for( j=bloco->jMin; j<=bloco->jMax; j++ ) {
                    //printf("HEHE j %d %lf \n", j, FuncaoParabolica(pMin, pMax, bloco->valorDirichlet, M.yMin + (j+0.5)*M.dy));
                    if( bloco->perfilInflow==PARABOLICO )
                        U[bloco->iMin][j] = FuncaoParabolica(pMin, pMax, bloco->valorDirichlet, M.y[j] + 0.5*M.dy[j], T, M);
                    else if( bloco->perfilInflow==RETO )
                        U[bloco->iMin][j] = bloco->valorDirichlet;
                }
            }
        }

        if( bloco->tipoV==BOUNDARY ) {
            contorno = bloco->tipoBoundaryV;
            if( (contorno==INFLOW) && (bloco->jMin==bloco->jMax) ) {
                pMin = M.xMin + (bloco->iMin)*M.dx[0];
                pMax = M.xMin + (bloco->iMax+1)*M.dx[0];

                if( bloco->perfilInflow==PARABOLICO_SIMETRICO_INICIO )
                    pMin = pMin - (pMax - pMin);
                else if( bloco->perfilInflow==PARABOLICO_SIMETRICO_FIM )
                    pMax = pMax + (pMax - pMin);

                for( i=bloco->iMin; i<=bloco->iMax; i++ ) {
                    if( bloco->perfilInflow==PARABOLICO || bloco->perfilInflow==PARABOLICO_SIMETRICO_INICIO || bloco->perfilInflow==PARABOLICO_SIMETRICO_FIM )
                        V[i][bloco->jMin] = FuncaoParabolica(pMin, pMax, bloco->valorDirichlet, M.xMin + (i+0.5)*M.dx[0], T, M); //ACHO QUE TA ERRADO CONFERIR
                    else if( bloco->perfilInflow==RETO )
                        V[i][bloco->jMin] = bloco->valorDirichlet;
                }
            }
        }

	}
}

PetscErrorCode InicializaEstruturasProjecao(SOLVER_PROJECAO *Solver, MALHA M)
{
    //Inicializando as matrizez de valores U, V, P, Psi
    double valorInicial = 0.0;
    Solver->UVelho = (double **)AlocaMatriz(M.Nx+1, M.Ny, sizeof(double), &valorInicial);
    Solver->UMeio = (double **)AlocaMatriz(M.Nx+1, M.Ny, sizeof(double), &valorInicial);
    Solver->UNovo = (double **)AlocaMatriz(M.Nx+1, M.Ny, sizeof(double), &valorInicial);
    Solver->VVelho = (double **)AlocaMatriz(M.Nx, M.Ny+1, sizeof(double), &valorInicial);
    Solver->VMeio = (double **)AlocaMatriz(M.Nx, M.Ny+1, sizeof(double), &valorInicial);
    Solver->VNovo = (double **)AlocaMatriz(M.Nx, M.Ny+1, sizeof(double), &valorInicial);
    Solver->PVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->Psi = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->PNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->WVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->WMeio = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->WNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);

    //Tensores viscoelastico (formulacao cartesiana)
    Solver->TxxVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TxyVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TyyVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TxxNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TxyNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TyyNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TxtVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TytVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TttVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TxtNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TytNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->TttNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);

    //Tensores viscoelastico (formulacao NSF)
    Solver->lambdaVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->muVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->nuVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->lambdaNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->muNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->nuNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->lambdaAux = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->nuAux = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);

    // Structure parameter for the EVPT model
    Solver->evptStructureVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->evptStructureNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);

    // Apparent viscosity for generalized newtonian models
    // Will always be 1 for purely newtonian/viscoelastic models
    double viscosidade_inicial = 1.0;
    Solver->MuNewtGenVelho = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &viscosidade_inicial);
    Solver->MuNewtGenNovo = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &viscosidade_inicial);


    Solver->norm_tau_dev = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);

    //Formulacao viscoelastica log conformation
    if( M.formulacaoVisco==VISCO_LOG ) {
        Solver->variaveis_log_conf.PsiVelho_xx = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiVelho_xy = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiVelho_yy = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiVelho_xt = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiVelho_yt = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiVelho_tt = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);

        Solver->variaveis_log_conf.PsiNovo_xx = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiNovo_xy = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiNovo_yy = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiNovo_xt = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiNovo_yt = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.PsiNovo_tt = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);

        Solver->variaveis_log_conf.O_11 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.O_12 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.O_13 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.O_21 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.O_22 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.O_23 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.O_31 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.O_32 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.O_33 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.lambda1 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.lambda2 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        Solver->variaveis_log_conf.lambda3 = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    }

    // Apenas pra imprimir no arquivo vtk e visualizar
    Solver->termoDominanteXX = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->termoDominanteXY = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
    Solver->termoDominanteYY = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);


    /// APenas pra imprimir no vtk e visualizar
    Solver->residualFaceU = (double **)AlocaMatriz(M.Nx+1, M.Ny, sizeof(double), &valorInicial);
    Solver->residualFaceV = (double **)AlocaMatriz(M.Nx, M.Ny+1, sizeof(double), &valorInicial);
    Solver->residualCell = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);

    //Armazenando a condicao inicial do problema
    //Estou colocando ZERO na pressao pq precisa de um valor pra usar no metodo da projecao
    CondicaoInicialProjecao(*Solver, M, 0.0);

    Solver->sistemasDestruidos = 1;

    return 0;
}

PetscErrorCode ConstroiSistemasLineares(SOLVER_PROJECAO *Solver, MALHA M)
{
    int qtd;
    int iStart, iEnd;
    PetscErrorCode ierr;

    Solver->sistemasDestruidos = 0;


    //Estruturas que armazenam as solucoes U e o vetor independente da equacao em x
    ierr = VecCreate(PETSC_COMM_WORLD, &(Solver->solU)); CHKERRQ(ierr);
	ierr = VecSetSizes(Solver->solU, PETSC_DECIDE, M.qtdIncognitasU); CHKERRQ(ierr);
    ierr = VecSetFromOptions(Solver->solU); CHKERRQ(ierr);
    ierr = VecDuplicate(Solver->solU, &(Solver->vetU)); CHKERRQ(ierr);

    //Solucoes V e vetor independente da eq em y
	VecCreate(PETSC_COMM_WORLD, &(Solver->solV));
	VecSetSizes(Solver->solV, PETSC_DECIDE, M.qtdIncognitasV);
    VecSetFromOptions(Solver->solV);
    VecDuplicate(Solver->solV, &(Solver->vetV));


    //Solucoes Psi e vet independente da poisson
	VecCreate(PETSC_COMM_WORLD, &(Solver->solPsi));
	VecSetSizes(Solver->solPsi, PETSC_DECIDE, M.qtdIncognitasP);
    VecSetFromOptions(Solver->solPsi);
    VecDuplicate(Solver->solPsi, &(Solver->vetP));


    //Solucoes W e vetor independente da equacao em theta
    if( (M.tipoCoord == AXI_CILINDRICO) && (M.interface.curvaInicial==NULL)  ) {
        VecCreate(PETSC_COMM_WORLD, &(Solver->solW));
        VecSetSizes(Solver->solW, PETSC_DECIDE, M.qtdIncognitasW);
        VecSetFromOptions(Solver->solW);
        VecDuplicate(Solver->solW, &(Solver->vetW));
    }



	//Estruturas que sao utilizadas no sistema linear da primeira equacao de movimento
	qtd = M.qtdIncognitasU;
	ierr = MatCreate(PETSC_COMM_WORLD, &(Solver->A_u)); CHKERRQ(ierr);//Cria a estrutura
    ierr = MatSetSizes(Solver->A_u, PETSC_DECIDE, PETSC_DECIDE, qtd, qtd); CHKERRQ(ierr); //Define as dimensoes
    ierr = MatSetFromOptions(Solver->A_u); CHKERRQ(ierr); //Finaliza a criacao da matriz
    ierr = VecGetOwnershipRange(Solver->solU, &iStart, &iEnd); CHKERRQ(ierr);
    ierr = InicializaMatrizCoeficientesU(Solver->A_u, qtd, M, iStart, iEnd); CHKERRQ(ierr);
    ierr = MontaMatrizCoeficientesU(Solver->A_u, M, Solver->MuNewtGenNovo); CHKERRQ(ierr);
    InicializaSolverKSP(&(Solver->kspU), Solver->A_u); CHKERRQ(ierr);



    //Estruturas que sao utilizadas no sistema linar da segunda equacao de movimento
    qtd = M.qtdIncognitasV;
	ierr = MatCreate(PETSC_COMM_WORLD, &(Solver->A_v)); CHKERRQ(ierr);//Cria a estrutura
    ierr = MatSetSizes(Solver->A_v, PETSC_DECIDE, PETSC_DECIDE, qtd, qtd); CHKERRQ(ierr);//Define as dimensoes
    ierr = MatSetFromOptions(Solver->A_v); CHKERRQ(ierr);//Finaliza a criacao da matriz
    ierr = VecGetOwnershipRange(Solver->solV, &iStart, &iEnd); CHKERRQ(ierr);
    ierr = InicializaMatrizCoeficientesV(Solver->A_v, qtd, M, iStart, iEnd); CHKERRQ(ierr);
    ierr = MontaMatrizCoeficientesV(Solver->A_v, M, Solver->MuNewtGenNovo); CHKERRQ(ierr);
    InicializaSolverKSP(&(Solver->kspV), Solver->A_v);

//    Estruturas utilizadas na equacao de Poisson
    qtd = M.qtdIncognitasP;
	ierr = MatCreate(PETSC_COMM_WORLD, &(Solver->A_p)); CHKERRQ(ierr);//Cria a estrutura
    ierr = MatSetSizes(Solver->A_p, PETSC_DECIDE, PETSC_DECIDE, qtd, qtd); CHKERRQ(ierr);//Define as dimensoes
    ierr = MatSetFromOptions(Solver->A_p); CHKERRQ(ierr);//Finaliza a criacao da matriz
    ierr = VecGetOwnershipRange(Solver->solPsi, &iStart, &iEnd); CHKERRQ(ierr);
    ierr = InicializaMatrizPoisson(Solver->A_p, qtd, M, iStart, iEnd); CHKERRQ(ierr);
    ierr = MontaMatrizPoisson(Solver->A_p, M, Solver->MuNewtGenNovo); CHKERRQ(ierr);
    InicializaSolverKSP(&(Solver->kspP), Solver->A_p);

    //Estruturas utilizadas na equacao de movimento em theta (soh resolve isso em escoamentos confinados...)
    if( (M.tipoCoord == AXI_CILINDRICO) && (M.interface.curvaInicial==NULL) ) {
        qtd = M.qtdIncognitasW;
        ierr = MatCreate(PETSC_COMM_WORLD, &(Solver->A_w)); CHKERRQ(ierr);//Cria a estrutura
        ierr = MatSetSizes(Solver->A_w, PETSC_DECIDE, PETSC_DECIDE, qtd, qtd); CHKERRQ(ierr);//Define as dimensoes
        ierr = MatSetFromOptions(Solver->A_w); CHKERRQ(ierr);//Finaliza a criacao da matriz
        ierr = VecGetOwnershipRange(Solver->solW, &iStart, &iEnd); CHKERRQ(ierr);
        ierr = InicializaMatrizCoeficientesW(Solver->A_w, qtd, M, iStart, iEnd); CHKERRQ(ierr);
        ierr = MontaMatrizCoeficientesW(Solver->A_w, M); CHKERRQ(ierr);
        InicializaSolverKSP(&(Solver->kspW), Solver->A_w);
    }

    return 0;
}

PetscErrorCode ExecutaPassoTemporalProjecao(SOLVER_PROJECAO *Solver, MALHA *M, int n, int rank, LISTA_CELULAS_SURFACE *ListaSurf, double *Residuo)
{
    PetscErrorCode ierr;
    LISTA_CELULAS_SURFACE *listaCelulasSurf;
    int qtdCurvas, i;
    CURVA *curva;
    static TIPO_CELULA **copiaMalha = NULL;
    static int reclassificaCelulasClean = 0;
    //static int quebrouJetting = 0, reclassificaCelulasClean = 0, frequenciaSuavizacao = 1000, mudaFrequenciaSuav = 0;
    //static int realizouUniao = 0;

    //printf("M->PontosBlocoParaview[0][0]=%lf\n", M->PontosBlocoParaview[0][0]);getchar();

    //usar o parametro abaixo se, por algum motivo, vc quiser fixar a matriz em problemas com superficie livre
    //se vc nao souber o que esta fazendo, deixe o valor ZERO neste parametro
    static int fixarMatriz = 0;

    //Vai ser usado em casos de superficie livre pra determinar se precisa recriar as matrizes do sistema linear ou nao
    if( !copiaMalha )
        copiaMalha = (TIPO_CELULA **)AlocaMatriz(M->Nx, M->Ny, sizeof(TIPO_CELULA), NULL);

    CondicaoContornoProjecao(Solver->UMeio, Solver->VMeio, Solver->WMeio, *M, (n+1)*M->dt); //Colocando as condicoes de contorno
    CondicaoContornoProjecao(Solver->UNovo, Solver->VNovo, Solver->WNovo, *M, (n+1)*M->dt);

//    IniciaTimer();
    /// ==== Primeiramente vou fazer a classificacao das celulas, ou seja, colocar as tags FULL, EMPTY e SURFACE
    /// ==== Em escoamentos confinados (sem superficie livre), isso eh feito na funcao LerModeloDoArquivo
    if( M->interface.curvaInicial!=NULL ) {
        /// ==== Fazendo a classificacao inicial das celulas (somente no primeiro passo)
        if( n==0 || reclassificaCelulasClean ) {
            InicializaMalhaComEmpty(M); //Limpando a malha. Colocando EMPTY em todas celulas

            // Contando quantas curvas tem
            qtdCurvas = 0;
            for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva )
                qtdCurvas++;
            listaCelulasSurf = (LISTA_CELULAS_SURFACE *)malloc(qtdCurvas*sizeof(LISTA_CELULAS_SURFACE));

            //Chamando a funcao que faz o preenchimento das celulas INICIAL
            i = 0;
            for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva )
                listaCelulasSurf[i++] = ClassificaCelulasEPontos(M, curva, Solver->UVelho, Solver->VVelho, Solver->PVelho);
            *ListaSurf = UniaoListasCelulaSurface(M, listaCelulasSurf, i);
            free(listaCelulasSurf);

            reclassificaCelulasClean = 0;
            DesenhaMalhaVTK(*M, 0);
            // DeuProblema("inicializou malha\n");
        }
        else if( !fixarMatriz ) { /// Fazendo apenas o REPREENCHIMENTO nos passos seguintes


            int j, mudou = 0;
            for( i=0; i<M->Nx; i++ )
                for( j=0; j<M->Ny; j++ )
                    copiaMalha[i][j] = M->celulas[i][j];

	    //VerificaInterface(M, Solver->UVelho, Solver->VVelho, Solver->UNovo, Solver->VNovo);
            ReclassificaCelulasEPontos(M, ListaSurf, n);
            for( i=0; i<M->Nx; i++ ) {
                for( j=0; j<M->Ny; j++ ) {
                    if( copiaMalha[i][j]!=M->celulas[i][j] )
                        mudou = 1;
                }
            }

            if( mudou && !(Solver->sistemasDestruidos) ) {
                DestroiSistemasLineares(Solver, *M);
            }
        }


        //// === Esse trecho faz mudancas topologicas de UNIAO. Esta em fase de teste isso ainda, entao deixei desativado aqui
        //static int embaracou = 0;
        double tempo = n*(M->dt);
        //static double merge_time = -1.0;
        //// merge_time = (merge_time>=0.0) ? merge_time : BuscaParametro(M, "merge_time");
        //merge_time = 0.0;
        //if( (tempo>merge_time) && (n%5==0) && !embaracou ) {
        //    embaracou = 1;

        //    // MudaNodeDaCurva(&M->interface, 1000.0, 0.0);

        //    DeletaCurvasFULL(M);
        //    if( RealizaMudancaTopologica_Uniao(M) ) {
        //        CheckUpInterface(M, "CheckUpInterface: apos realizar uniao\n");
        //        PrintDebug("SEGUNDA UNIAO!!! \n\n\n\n");
        //        RealizaMudancaTopologica_Uniao(M);
        //        CheckUpInterface(M, "CheckUpInterface: apos realizar uniao\n");
        //        RemoveParticulasFULL(M);
        //        CheckUpInterface(M, "CheckUpInterface: apos realizar uniao\n");
        //       AdicionaERemoveParticulas(M, 0.05*M->minDxDy, 0);

        //        // ImprimeInterfaceVTK(*M, 100000);
        //        // DeuProblema("parar\n");
        //    }
            
            
            
        //    // ImprimeInterfaceVTK(*M, 100000);
        //    if( RealizaEmbaracamentoForcado(M) ) {
        //        // ImprimeInterfaceVTK(*M, 200000);

        //        PrintDebug("Comecou uniao.\n");
        //        PONTO nodes[2];
        //        DesenhaMalhaVTK(*M, 0);
        //        RealizaDesembaracamento(&(M->interface), nodes);
        //        PrintDebug("Realiozu uniao.\n");


        //        // ImprimeInterfaceVTK(*M, 300000);
        //        // getchar();

        //        CheckUpInterface(M, "CheckUpInterface: apos realizar uniao\n");
        //        ArrumaMudancaTopologicaAxissimetrica(M);
        //        ArrumaMudancaTopologicaAxissimetrica(M);
        //        // ImprimeInterfaceVTK(*M, 400000);
        //        CheckUpInterface(M, "CheckUpInterface: apos realizar uniao\n");
        //    }
        //}


        ///// === Este trecho faz mudancas topologicas de QUEBRA.
        static int quebrouRecentemente = 0;
        static double tempo_minimo_quebra = 0.0;
        static int frequenciaQuebra = 39;

        frequenciaQuebra = (quebrouRecentemente) ? 1 : 39;
        quebrouRecentemente = (quebrouRecentemente>0) ? quebrouRecentemente-1 : 0;

        if( (!quebrouRecentemente) && (n%frequenciaQuebra==0) && (tempo>tempo_minimo_quebra) ) {

            // IRINEU, IMPORTANTE: 
            // A variavel dxQuebra determina a menor espessura aceitavel de um filamento antes de ele ser quebrado em dois pedacos
            // Se o codigo encontrar um filamento com espessura menor que dxQuebra, a funcao "RealizaEmbaracamentoForcado_Quebra2" quebra o filamento em dois pedacos
            // Talvez seja bom fazer testes ajustando este parametro
            static double dxQuebra = 0.05;
            dxQuebra = 1.25*M->minDxDy;

            
            // IRINEU, IMPORTANTE: Dentro da funcao "RealizaEmbaracamentoForcado_Quebra2" tem um parametro tambem que vc precisa olhar (qtdPontosEmbaraca).
            PONTO *pRemover = NULL;
            int quebrouSim = RealizaEmbaracamentoForcado_Quebra2(M, 1.0*dxQuebra, 5.0*dxQuebra, &pRemover, NULL);

            if( quebrouSim ) {
                PrintDebug("# Realizou mudanca topologica\n");
                
                // IRINEU, IMPORTANTE:
                // A variabel "quebrouRecentemente" impede que duas mudancas topologicas sejam realizadas uma logo apos a outra. Isto evita muitas goticulas sendo criadas
                // Esse numero enorme que eu coloquei faz com que somente 1 mudanca topologica aconteca durante toda a simulacao
                // Se voce colocar o numero "1000", por exemplo, vai ter um intervalo de 1000 passos temporais entre duas mudancas topologicas consecutivas
                // Eu tenho certeza que vc vai precisar mexer nesse parametro. Principalmente quando for necessario gerar multiplas goticulas na simulacao
                quebrouRecentemente = 1000000000;
            }

           
        }
   }

    

    // if( n%500==0 )
    //     DeletaCurvasFULL(M);

    // CheckUpInterface(M, "teste 1");

    /// ==== Colocando condicao inicial se estiver no primeiro passo
    if( n==0 ) {
        CondicaoInicialProjecao(*Solver, *M, 0.0);
        if( M->intervaloVTK_prop ) {
            ImprimeArquivoVTK(Solver->UVelho, Solver->VVelho, Solver->WVelho, Solver->PVelho,
                              Solver->TxxVelho, Solver->TxyVelho, Solver->TyyVelho, Solver->TxtVelho, Solver->TytVelho, Solver->TttVelho,
                              Solver->lambdaVelho, Solver->muVelho, Solver->nuVelho,
                              Solver->evptStructureVelho, Solver->norm_tau_dev,
                              Solver->MuNewtGenNovo,
                              Solver->termoDominanteXX, Solver->termoDominanteXY, Solver->termoDominanteYY,
                              Solver->residualFaceU, Solver->residualFaceV, Solver->residualCell,
                              *M, 0);
        }

        DesenhaMalhaVTK(*M, 0);
    }
    
//    IniciaTimer();
    /// ==== Agora aplicando a condicao tangencial de sup.livre pra calcular U e V no contorno
    if( M->interface.curvaInicial!=NULL )
        CondicaoContornoTangencial(n, M, *ListaSurf, Solver->UVelho, Solver->VVelho, Solver->UNovo, Solver->VNovo, Solver->PNovo, Solver->TxxVelho, Solver->TxyVelho, Solver->TyyVelho, &(Solver->kspS), &(Solver->A_s), &(Solver->vetS), &(Solver->solS));
//    PrintDebug("Contorno Tangencial: %lf segundos\n", FinalizaTimer());


//    PrintDebug("chegou 2\n");
//    IniciaTimer();
    /// ==== Constroi a ordenacao das incognitas nos sistemas lineares que virao a seguir
    if( M->interface.curvaInicial!=NULL && (n==0 || !fixarMatriz) && (Solver->sistemasDestruidos==1) )
        InicializaIndiceIncognitas(M);
//    PrintDebug("Indice incognitas: %lf segundos\n", FinalizaTimer());

    // Quando uma celula muda de EMPTY pra SURFACE
    // InicializaNovasPropriedadesSURFACE(M, ListaSurf, Solver->evptStructureVelho);

    // In the case of generalized Newtonian models, we update the apparent viscosity in the whole mesh
    // GeneralizedNewtonian_Viscosity(*M, n, Solver->UVelho, Solver->VVelho, Solver->MuNewtGenNovo);

    /// === The viscoplastic model is implemented as a generalized Newtonian model.
    /// === We calculate the non-newtonian part of the stess tensor in this case
    if( M->tipo_modelo==MODELO_VP ) {
        GeneralizedNewtonian_Tensor(*M, n, Solver->UVelho, Solver->VVelho, Solver->MuNewtGenNovo,
                                    Solver->TxxVelho, Solver->TxyVelho, Solver->TyyNovo);
        GeneralizedNewtonian_Tensor(*M, n, Solver->UVelho, Solver->VVelho, Solver->MuNewtGenNovo,
                                    Solver->TxxNovo, Solver->TxyNovo, Solver->TyyNovo);
    }

    // debug printing
    // int j;
    // for( j=0; j<M->Ny; j++ ) {
    //     double y = 0.5*( M->y[j] + M->y[j+1] );
    //     // double esq = Solver->MuNewtGenNovo[1][j];
    //     // double dir = Solver->MuNewtGenNovo[M->Nx-2][j];
    //     double esq = debugsCell[1][j];
    //     double dir = debugsCell[M->Nx-2][j];
    //     PrintDebug("surfs %10lf: %20.15lf %20.15lf %e\n", y, esq, dir, fabs(esq-dir));
    // }

    // In the case of EVPT model, we also solve an equation to update the structure parameter
    if( M->tipo_modelo==MODELO_EVPT ) {
        Evolucao_StructureParameter(*M, Solver->UVelho, Solver->VVelho, Solver->PVelho,
                                    Solver->TxxVelho, Solver->TxyVelho, Solver->TyyVelho,
                                    Solver->evptStructureVelho, Solver->evptStructureNovo);
    }

//    IniciaTimer();
    /// ==== Constroi os tres sistemas lineares (Mov. Direcao X, Mov. Direcao Y, e Poisson
    if( (n==0) || ((M->interface.curvaInicial!=NULL) && !fixarMatriz && Solver->sistemasDestruidos==1) ) {
        ConstroiSistemasLineares(Solver, *M);
//    PrintDebug("Constroi sistemas: %lf segundos\n", FinalizaTimer());
    }


//    IniciaTimer();
    /// ==== Resolvendo a primeira equacao de qtd. de movimento para calcular UMeio
    ierr = VetorIndependenteU(Solver->vetU, Solver->UVelho, Solver->VVelho, Solver->WVelho, Solver->UNovo, Solver->VNovo, Solver->WNovo,
                              Solver->PVelho, Solver->TxxVelho, Solver->TxyVelho, Solver->TttVelho, Solver->MuNewtGenNovo, n, *M); CHKERRQ(ierr);  
    ierr = KSPSolve(Solver->kspU, Solver->vetU, Solver->solU); CHKERRQ(ierr);
    ierr = CopiaVecUPraMatriz(Solver->solU, Solver->UMeio, *M, rank); CHKERRQ(ierr);
//    PrintDebug("SolveU: %lf segundos\n", FinalizaTimer());


    /// ==== Resolvendo a segunda equacao de qtd. de movimento para calcular VMeio
//    IniciaTimer();
    ierr = VetorIndependenteV(Solver->vetV, Solver->UVelho, Solver->VVelho, Solver->UNovo, Solver->VNovo, Solver->PVelho, Solver->TxyVelho, Solver->TyyVelho, Solver->MuNewtGenNovo, *M); CHKERRQ(ierr);
    ierr = KSPSolve(Solver->kspV, Solver->vetV, Solver->solV); CHKERRQ(ierr);
    ierr = CopiaVecVPraMatriz(Solver->solV, Solver->VMeio, *M, rank); CHKERRQ(ierr);
//    PrintDebug("SolveV: %lf segundos\n", FinalizaTimer());


    /// ==== Axissimetrico: Resolvendo a terceira equacao de qtd. de movimento para calcular WMeio
    /// ==== Soh resolve esta equacao em escoamentos confinados. Em free-surface, considero u_theta = 0
    if( (M->tipoCoord==AXI_CILINDRICO) && (M->interface.curvaInicial==NULL) ) {
        ierr = VetorIndependenteW(Solver->vetW, Solver->UVelho, Solver->VVelho, Solver->WVelho, Solver->PVelho, Solver->TxtVelho, Solver->TytVelho, *M); CHKERRQ(ierr);
        ierr = KSPSolve(Solver->kspW, Solver->vetW, Solver->solW); CHKERRQ(ierr);
        ierr = CopiaVecPsiPraMatriz(Solver->solW, Solver->WMeio, *M, rank); CHKERRQ(ierr);
    }

    

//    IniciaTimer();
    /// ==== Resolvendo a equacao de Poisson para calcular Psi
    ierr = VetorIndependentePoisson(Solver->vetP, Solver->UMeio, Solver->VMeio, Solver->PVelho, Solver->PNovo, Solver->TxxVelho, Solver->TxyVelho, Solver->TyyVelho, Solver->Psi, Solver->MuNewtGenNovo, *M, ListaSurf); CHKERRQ(ierr);
    ierr = KSPSolve(Solver->kspP, Solver->vetP, Solver->solPsi); CHKERRQ(ierr);
    ierr = CopiaVecPsiPraMatriz(Solver->solPsi, Solver->Psi, *M, rank); CHKERRQ(ierr);
//    PrintDebug("SolvePoisson: %lf segundos\n", FinalizaTimer());


    /// ==== Usando as formulas de atualizacao para calcular UNovo, VNovo, PNovo
//    IniciaTimer();
    CalculaVelocidadeNova(Solver->UNovo, Solver->VNovo, Solver->WNovo, Solver->UMeio, Solver->VMeio, Solver->WMeio, Solver->Psi, *M);
    CalculaPressaoNova(Solver->PNovo, Solver->PVelho, Solver->Psi, *M);
//    PrintDebug("PressaoNova: %lf segundos\n", FinalizaTimer());


    /// ==== Resolve as equacoes constitutivas para atualizar os tensores (calcula TxxNovo, TxyNovo, TyyNovo)
//    IniciaTimer();
    EquacoesConstitutivas(Solver->lambdaVelho, Solver->muVelho, Solver->nuVelho, Solver->lambdaNovo, Solver->muNovo, Solver->nuNovo, Solver->lambdaAux, Solver->nuAux,
                            Solver->TxxVelho, Solver->TxyVelho, Solver->TyyVelho, Solver->TxtVelho, Solver->TytVelho, Solver->TttVelho,
                            Solver->TxxNovo, Solver->TxyNovo, Solver->TyyNovo, Solver->TxtNovo, Solver->TytNovo, Solver->TttNovo,
                            Solver->evptStructureNovo, Solver->norm_tau_dev,
                            &(Solver->variaveis_log_conf),
                            Solver->UVelho, Solver->VVelho, Solver->WVelho, Solver->UNovo, Solver->VNovo, Solver->WNovo, Solver->PNovo,
                            Solver->termoDominanteXX, Solver->termoDominanteXY, Solver->termoDominanteYY,
                            n, *M);
//    PrintDebug("Constitutivas: %lf segundos\n", FinalizaTimer());

//    static double **residuoEq = NULL, valorInicial = 0.0;
//    if( residuoEq == NULL )
//        residuoEq = (double **)AlocaMatriz(M->Nx + 1, M->Ny, sizeof(double), &valorInicial);
//
//    CalculaResiduoMomentoU(M, Solver->UVelho, Solver->UMeio, Solver->UNovo,
//                            Solver->VVelho, Solver->VMeio, Solver->VNovo,
//                            Solver->PVelho, Solver->PNovo,
//                            Solver->TxxVelho, Solver->TxxNovo, Solver->TxyVelho, Solver->TxyNovo, Solver->TttVelho, Solver->TttNovo,
//                            residuoEq);

//    PrintDebug("RESIDUO: %lf\n", residuoEq[M->Nx/2][M->Ny/2]);


//    IniciaTimer();
    /// ==== Faz uma diferenca entre a nova solucao e a solucao do passo anterior p/ medir se estacionou
    *Residuo = NormaEntreSolucoesU(*Solver, *M);
//    PrintDebug("NormaSolucoes: %lf segundos\n", FinalizaTimer());

//    IniciaTimer();
    /// ==== Realizando a adveccao das particulas
    if( M->interface.curvaInicial!=NULL && !fixarMatriz ) {
        //Copia todas as solucoes das matrizes Novo para as matrizes Velho
        //CopiaSolucaoNovaParaVelha(Solver, *M);
        //Aplica novamente a condicao tangencial no contorno
        CondicaoContornoTangencial(n, M, *ListaSurf, Solver->UNovo, Solver->VNovo, Solver->UNovo, Solver->VNovo, Solver->PNovo, Solver->TxxVelho, Solver->TxyVelho, Solver->TyyVelho, &(Solver->kspS), &(Solver->A_s), &(Solver->vetS), &(Solver->solS));
        // CheckUpInterface(M, "teste 1.5");


        //Advecta as particulas
        AdvectaInterface(M, Solver->UVelho, Solver->VVelho, Solver->UNovo, Solver->VNovo);
        // CheckUpInterface(M, "teste 2");


        if( n>10 && (n%45==0) ) {
            // RemoveParticulasFULL(M);
            // CheckUpInterface(M, "teste 3");
            // RemoveParticulasNaSimetriaSURFACE(M);
            // CheckUpInterface(M, "teste 4");
        }

//        PrintDebug("chegou 1\n");

        // Suavizando
        if( (M->freqTSUR>0) && (n>5) && !(n%(M->freqTSUR)) )
            SuavizaInterfaceTSUR(M);
	
	//Teste July 23, 2025
	//printf("\n\nVerifica ANTES de Adicionar Particula\n");
	//VerificaInterface(M, Solver->UVelho, Solver->VVelho, Solver->UNovo, Solver->VNovo);


        // CheckUpInterface(M, "teste 5");

//        PrintDebug("chegou 2\n");

//        RemoveCurvasMuitoPequenas(M, 10);

        //Adiciona e remove particulas, se necessario
	//printf("minDxDy=%lf\n", M->minDxDy);getchar();
        AdicionaERemoveParticulas(M, 0.05*M->minDxDy, 0);

	//Teste July 23, 2025
	//printf("\n\nVerifica APOS de Adicionar Particula\n");
	//VerificaInterface(M, Solver->UVelho, Solver->VVelho, Solver->UNovo, Solver->VNovo);


//        PrintDebug("chegou 3\n");

//        ImprimeInterfaceVTK(*M, 100000);
//        DeuProblema("parou\n");

//        PrintDebug("chegou 10\n");
    }
//    PrintDebug("Adveccao Particulas: %lf segundos\n", FinalizaTimer());


    /// ==== Copia todas as solucoes das matrizes Novo para as matrizes Velho
//    IniciaTimer();
    CopiaSolucaoNovaParaVelha(Solver, *M);
//    PrintDebug("CopiaSolucaoNovaPraVelha: %lf segundos\n", FinalizaTimer());


    /// === Se nao tiver mais nenhuma celula SURFACE, vou fixar a matriz
    /// === Pensa no canal de placas paralelas sendo preenchido. Uma hora ele fica todo cheio e nao vai mais ter celula SURFACE
    if( ListaSurf->prim==NULL )
        fixarMatriz = 1;

    /// === Se tiver superficie livre, destroi as estruturas de todos sistemas lineares
   if( M->interface.curvaInicial!=NULL && !fixarMatriz )
       DestroiSistemasLineares(Solver, *M);
    return 0;
}

void ProcuraParticulaNaN(MALHA *M, const char *Mensagem)
{
    CURVA *curva;
    PONTO *p;

    for( curva=M->interface.curvaInicial; curva; curva=curva->proxCurva ) {
        LOOP_CURVA(curva) {
            p = loop_stuff.ponto;

            if( isnan(p->x) || isinf(p->x) )
                DeuProblema(Mensagem);
            if( isnan(p->y) || isinf(p->y) )
                DeuProblema(Mensagem);
        }
    }

    return;
}

void CalculaVelocidadeNova(double **UNovo, double **VNovo, double **WNovo, double **UMeio, double **VMeio, double **WMeio, double **Psi, MALHA M)
{
    int i, j, linha;
    static double psiDir = 0.0, psiEsq = 0.0, psiBaixo = 0.0, psiCima = 0.0, psiCentro = 0.0;
    double h1, h2;
    double coefMenos, coefMais, coefCentro;


    //Primeira equacao (calculando U)
    for( linha=0; linha<M.qtdIncognitasU; linha++ ) {
        i = M.indicesU[linha].i;
        j = M.indicesU[linha].j;

        //Casos que nao precisa calcular pq eh dirichlet na velocidade U
        if( M.pontosU[i][j].tipo == EMPTY )
            DeuProblema("\n\nPROBLEMA\n\n");
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==INFLOW )
            continue;
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==NOSLIP )
            continue;
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==SLIP )
            continue;
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA )
            continue;
        else if( M.pontosU[i][j].tipo==SURFACE )
            continue;


        //Caso o Psi da direita esteja fora do dominio
        if( M.pontosU[i][j].tipo == BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
            if( M.pontosU[i][j].tipoBoundary==NEUMANN )
                psiDir = - Psi[i-1][j]; //Neumann na velocidade vira dirichlet no psi
            else if( M.pontosU[i][j].tipoBoundary==SIMETRIA )
                psiDir = Psi[i-1][j];
            else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
                psiDir = Psi[0][j];
            else if( M.pontosU[i][j].tipoBoundary==KNOWN_PRESSURE )
                psiDir = 0.0;
            else
                DeuProblema("PROBLEMA DIREITA ATUALIZACAO VEL\n");
        }
        else
            psiDir = Psi[i][j];



        //Caso o Psi da esquerda esteja fora
        if( M.pontosU[i][j].tipo == BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
            if( M.pontosU[i][j].tipoBoundary==NEUMANN )
                psiEsq = - Psi[i][j];
            else if( M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA )
                psiEsq = Psi[i][j];
            else if( M.pontosU[i][j].tipoBoundary==SIMETRIA )
                psiEsq = Psi[i][j];
            else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
                psiEsq = Psi[M.Nx-1][j];
            else if( M.pontosU[i][j].tipoBoundary==KNOWN_PRESSURE )
                psiEsq = 0.0;
            else
                DeuProblema("PROBLEMA ESQUERDA ATUALIZACAO VEL %d %d %d\n", i, j, M.pontosU[i][j].tipoBoundary);
        }
        else
            psiEsq = Psi[i-1][j];


        h1 = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        h2 = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];
        coefMenos = -h2/(h1*(h1+h2));
        coefMais = h1/(h2*(h1+h2));
        coefCentro = (h2-h1)/(h1*h2);
        psiCentro = Interpolacao(h1, h2, psiEsq, psiDir);

        UNovo[i][j] = UMeio[i][j] - (coefMenos*psiEsq + coefCentro*psiCentro + coefMais*psiDir);
    }


    //Segunda equacao (calculando V)
    for( linha=0; linha<M.qtdIncognitasV; linha++ ) {
        i = M.indicesV[linha].i;
        j = M.indicesV[linha].j;

        if( M.pontosV[i][j].tipo == EMPTY )
            DeuProblema("\n\nPROBLEMA\n\n");
        else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==INFLOW )
            continue;
        else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==NOSLIP )
            continue;
        else if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].tipoBoundary==SLIP )
            continue;
        else if( M.pontosV[i][j].tipo==SURFACE )
            continue;

        //Colocando contorno la em cima do dominio
        if( M.pontosV[i][j].tipo == BOUNDARY && M.pontosV[i][j].extremo==CIMA ) {
            if( M.pontosV[i][j].tipoBoundary == NEUMANN  )
                psiCima = - Psi[i][j-1];
            else if( M.pontosV[i][j].tipoBoundary == SIMETRIA  )
                psiCima = Psi[i][j-1];
            else
                DeuProblema("PROBLEMA CIMA ATUALIZACAO VEL\n");
        }
        else
            psiCima = Psi[i][j];


        //Colocando contorno la em baixo do dominio
        if( M.pontosV[i][j].tipo == BOUNDARY && M.pontosV[i][j].extremo==BAIXO ) {
            if( M.pontosV[i][j].tipoBoundary == NEUMANN )
                psiBaixo = - Psi[i][j];
            else if( M.pontosV[i][j].tipoBoundary == SIMETRIA )
                psiBaixo = Psi[i][j+1];
            else
                DeuProblema("PROBLEMA BAIXO ATUALIZACAO VEL\n");
        }
        else
            psiBaixo = Psi[i][j-1];

        h1 = (j==0) ? 0.5*M.dy[j] : 0.5*M.dy[j-1];
        h2 = (j==M.Ny) ? 0.5*M.dy[j-1] : 0.5*M.dy[j];
        coefMenos = -h2/(h1*(h1+h2));
        coefMais = h1/(h2*(h1+h2));
        coefCentro = (h2-h1)/(h1*h2);
        psiCentro = Interpolacao(h1, h2, psiBaixo, psiCima);

        VNovo[i][j] = VMeio[i][j] - (coefMenos*psiBaixo + coefCentro*psiCentro + coefMais*psiCima);
    }

    //Terceira equacao (calculando W)
    if( M.tipoCoord==AXI_CILINDRICO ) {
        for( linha=0; linha<M.qtdIncognitasW; linha++ ) {
            i = M.indicesW[linha].i;
            j = M.indicesW[linha].j;

            WNovo[i][j] = WMeio[i][j];
        }
    }

    return;
}

void CalculaPressaoNova(double **PNovo, double **P, double **Psi, MALHA M)
{
    int i, j, linha;

    // VERIFICAR A MODIFICACAO DISSO AQUI DEPOIS PRA ACRESCENTAR
    for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
        i = M.indicesP[linha].i;
        j = M.indicesP[linha].j;

        if( M.pontosP[i][j].tipo==EMPTY )
            DeuProblema("\n PROBLEMA \n");
        else if( M.pontosP[i][j].tipo==FULL || M.pontosP[i][j].tipo==SURFACE )
            PNovo[i][j] = P[i][j] + Psi[i][j]/M.dt;// - (teta/Re)*Lp(Psi, i, j);
        else DeuProblema("\n Problema na funcao CalculaPressaoNova. \n");
    }

    return;
}

extern double *hVelho;

PetscErrorCode AdvectaInterface(MALHA *M, double **UVelho, double **VVelho, double **UNovo, double **VNovo)
{
    CURVA *curva, *proxCurva = NULL;
    PONTO *ponto, *proxPonto;
    double velU1, velV1;
    int iVelho, iNovo, jVelho, jNovo;
    double xVelho, yVelho;
    int sensor=0, sensor2=0;
    double epsilon=1.e-12;
    //double menor = 100.0, dif, velocidade=0, valor;

    for( curva=M->interface.curvaInicial; curva!=NULL; curva=proxCurva) {
        proxCurva = curva->proxCurva;

        ponto = curva->pontoInicial;
        do {
            proxPonto = ponto->prox;

            EncontraCelulaDaParticula(M, ponto->x, ponto->y, &iVelho, &jVelho);
	    xVelho = ponto->x;
	    yVelho = ponto->y;
            /// === Um teste: remover particulas que estao em celulas EMPTY
            /// === Isto acontece em alguns casos degenerados
            if( M->celulas[iVelho][jVelho]==EMPTY ) {
//                PrintDebug("LEMBRAR DE ARRUMAR A REMOCAO O CASO DENEGERADO NA FUNCAO ADVECTAINTERFACE\n\n");
//                PrintDebug("REMOVEU KKKKK %lf %lf\n\n", ponto->x, ponto->y);
//                getchar();

                RemovePontoDaCurva(curva, ponto);
                ponto = proxPonto;


                // Se por acaso removeu a curva inteira... tira a curva da interface
                // Isso acontece as vezes raramente em casos de splashing com gotinhas pequenas (horrivel)
                if( curva->pontoInicial==NULL ) {
                    /// Deletando essa curva da interface
                    if( curva->curvaAnt )
                        curva->curvaAnt->proxCurva = proxCurva;
                    if( proxCurva )
                        proxCurva->curvaAnt = curva->curvaAnt;
                    if( M->interface.curvaInicial==curva )
                        M->interface.curvaInicial = proxCurva;
                    LiberaMemoriaCurva(curva);
                }

                continue;
            }

            //Nao vai ser advectado
            if( ponto->fixo ) {
                // PrintDebug("FIXO!!!\n");
                ponto = proxPonto;
                continue;
            }

            /// ===== Euler modificado
            VelocidadeEmUmaParticula(*M, ponto, UNovo, VNovo, &velU1, &velV1);
            AdvectaParticula(ponto, velU1, velV1, M->dt, &(ponto->x), &(ponto->y));

            /// ===== RK21
//            PONTO pontoAux;
//            double velU2, velV2;
//            pontoAux.x = ponto->x;
//            pontoAux.y = ponto->y;
//            VelocidadeEmUmaParticula(*M, &pontoAux, UVelho, VVelho, &velU1, &velV1);
//            AdvectaParticula(&pontoAux, velU1, velV1, M->dt, &(pontoAux.x), &(pontoAux.y));
//            VelocidadeEmUmaParticula(*M, &pontoAux, UNovo, VNovo, &velU2, &velV2);
//            AdvectaParticula(ponto, 0.5*(velU1 + velU2), 0.5*(velV1 + velV2), M->dt, &(ponto->x), &(ponto->y));


            //Verifica se a particula saiu do dominio
            if( ponto->x >= M->xMax ) {
                ponto->x = M->xMax - 1e-8;
//                ponto->outflow = 1;
//                ponto->fixo = 1;
            }
            else if( ponto->x <= M->xMin )
                ponto->x = M->xMin + 1e-8;

            if( ponto->y >= M->yMax )
                ponto->y = M->yMax - 1e-8;
            else if( ponto->y <= M->yMin )
                ponto->y = M->yMin + 1e-8;

            EncontraCelulaDaParticula(M, ponto->x, ponto->y, &iNovo, &jNovo);

	    //Identifiquei que no problema do 'droplet orifice' os indices iVelho e jVelho haviam sido sobrescritos com o novo ponto.
            //EncontraCelulaDaParticula(M, ponto->ant->x, ponto->ant->y, &iVelho, &jVelho);

	    //if (iVelho == 5 && jVelho==119){
		    //printf("ponto->x=%f, ponto->y=%f\n", ponto->x, ponto->y);
		    //printf("VELHO: iNovo=%d, iVelho=%d, jNovo=%d, jVelho=%d\n", iNovo, iVelho, jNovo, jVelho);getchar();
	    //}
	    //if (iNovo == 5 && jVelho==119){
		    //printf("ponto->x=%f, ponto->y=%f\n", ponto->x, ponto->y);
		    //printf("NOVO: iNovo=%d, iVelho=%d, jNovo=%d, jVelho=%d\n", iNovo, iVelho, jNovo, jVelho);getchar();
	    //}

	    //printf("M->x1=%f, M->x2=%f, M->y1=%f, M->y3=%f\n", M->x1, M->x2, M->y1, M->y3);getchar();
	    sensor2=0;
	    if ( ((ponto->x >= (M->x1-epsilon)) && (ponto->x <= (M->x2+epsilon))) && ((ponto->y >= (M->y3-epsilon)) && (ponto->y <= (M->y1+epsilon))) ){
		    // printf("\n\nANTES: Particula marcadora dentro de celula BOUNDARY\n");
		    // printf("celula[%d][%d] = %d\n", iNovo, jNovo, M->celulas[iNovo][jNovo]);
		    // printf("iVelho=%d, jVelho=%d, xVelho=%f, yVelho=%f\n", iVelho, jVelho, xVelho, yVelho);
		    // printf("iNovo=%d, jNovo=%d, xNovo=%f, yNovo=%f\n", iNovo, jNovo, ponto->x, ponto->y);
		    sensor2=1;
	    }
	    //if (M->celulas[5][119] == 4){
		    //printf("\ncelula BOUNDARY foi transformada em SURFACE\n");
		    //printf("iVelho=%d, jVelho=%d, xVelho=%f, yVelho=%f\n", iVelho, jVelho, xVelho, yVelho);
		    //printf("iNovo=%d, jNovo=%d, xNovo=%f, yNovo=%f\n", iNovo, jNovo, ponto->x, ponto->y);
		    //getchar();
	    //}

	    sensor=0;
            // Se a particula passou por uma parede vertical
            if( (iNovo>iVelho) && M->pontosU[iNovo][jVelho].tipo==BOUNDARY && sensor2==1){
		//printf("particula atravessou uma BOUNDARY\n");getchar();
                ponto->x = M->x[iNovo] - 1e-8;
		sensor=1;
	    }else if( (iNovo<iVelho) && M->pontosU[iVelho][jVelho].tipo==BOUNDARY && sensor2==1){
                ponto->x = M->x[iVelho] + 1e-8;
		sensor=2;
	    }

            // Se a parede passou por uma parede horizontal
            if( (jNovo>jVelho) && M->pontosV[iVelho][jNovo].tipo==BOUNDARY && sensor2==1){
                ponto->y = M->y[jNovo] - 1e-8;
		sensor=3;
	    }
            else if( (jNovo<jVelho) && M->pontosV[iVelho][jVelho].tipo==BOUNDARY && sensor2==1){
                ponto->y = M->y[jVelho] + 1e-8;
		sensor=4;
	    }

            // Se a particula passou por uma parede na diagonal
            if( (jNovo<jVelho) && (iNovo<iVelho) && (M->pontosV[iNovo][jNovo+1].tipo==BOUNDARY) && (M->pontosU[iNovo+1][jNovo].tipo==BOUNDARY) && sensor2==1) {
                ponto->x = M->x[iNovo + 1] + 1e-8;
                ponto->y = M->y[jNovo + 1] + 1e-8;
		sensor=5;
            }
            else if( (jNovo<jVelho) && (iNovo>iVelho) && (M->pontosV[iNovo][jNovo+1].tipo==BOUNDARY) && (M->pontosU[iNovo][jNovo].tipo==BOUNDARY) && sensor2==1) {
                ponto->x = M->x[iNovo] - 1e-8;
                ponto->y = M->y[jNovo + 1] + 1e-8;
		sensor=6;
            }

	    //sensor2=0;
	    //if (ponto->x >= 0.25 && (ponto->y >= -0.2 && ponto->y <= 0.0)){
	    //if ( ((ponto->x >= (M->x1-epsilon)) && (ponto->x <= (M->x2+epsilon))) && ((ponto->y >= (M->y3-epsilon)) && (ponto->y <= (M->y1+epsilon))) ){
		    //printf("\n\nNAO CORRIGIU A CELULA: Particula marcadora dentro de celula BOUNDARY\n");
		    //printf("iVelho=%d, jVelho=%d, xVelho=%f, yVelho=%f\n", iVelho, jVelho, xVelho, yVelho);
		    //printf("iNovo=%d, jNovo=%d, xNovo=%f, yNovo=%f\n", iNovo, jNovo, ponto->x, ponto->y);
		    //sensor2=1;
		    //exit(1);
	    //}

	    if (sensor!=0 && sensor2==1){
		    // printf("\n\nMODIFICADA: Particula foi modificada corretamente\n");
		    // printf("iVelho=%d, jVelho=%d, xVelho=%f, yVelho=%f\n", iVelho, jVelho, xVelho, yVelho);
		    // printf("iNovo=%d, jNovo=%d, xNovo=%f, yNovo=%f\n", iNovo, jNovo, ponto->x, ponto->y);
		    // printf("sensor=%d\n", sensor);
		    // printf("sensor2=%d\n", sensor2);
		    if (M->celulas[iNovo][jNovo] == FULL || M->celulas[iNovo][jNovo] == SURFACE){//deveria ser boundary para estar certo
			    printf("Particula modificada erroneamente, pois celula eh FULL\n");			  
			    exit(1);
		    }
		    //getchar();
	    }

            ponto = proxPonto;
        } while( ponto!=curva->pontoInicial );

    }

    return 0;
}

PetscErrorCode VerificaInterface(MALHA *M, double **UVelho, double **VVelho, double **UNovo, double **VNovo)
{
    //printf("Dentro de verifica interface\n");getchar();

    CURVA *curva, *proxCurva = NULL;
    PONTO *ponto, *proxPonto;
    double velU1, velV1;
    int iVelho, iNovo, jVelho, jNovo;
    //double menor = 100.0, dif, velocidade=0, valor;

    //printf("M->celulasParede[%d][%d]=%d\n", 5, 119, M->celulasParede[5][119]);

    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva) {
	LOOP_CURVA(curva) {
        	ponto = loop_stuff.ponto;

        	//proxCurva = curva->proxCurva;

        	//ponto = curva->pontoInicial;
        //do {
            //proxPonto = ponto->prox;

            EncontraCelulaDaParticula(M, ponto->x, ponto->y, &iNovo, &jNovo);
	    //printf("VerificaInterface: iNovo=%d, jNovo=%d, xNovo=%f, yNovo=%f\n", iNovo, jNovo, ponto->x, ponto->y);

	    
	    if (ponto->x >= M->x1 && (ponto->y >= M->y3 && ponto->y <= M->y1)){
		    printf("\n\nVERIFICA INTERFACE\nParticula marcadora dentro de celula BOUNDARY\n");
		    printf("M->celulas[%d][%d] = %d\n", iNovo, jNovo, M->celulas[iNovo][jNovo]);
		    printf("M->celulasParede[%d][%d]=%d\n", iNovo, jNovo, M->celulasParede[iNovo][jNovo]);
		    printf("iNovo=%d, jNovo=%d, xNovo=%f, yNovo=%f\n", iNovo, jNovo, ponto->x, ponto->y);
		    exit(1);
	    }

            //ponto = proxPonto;
        } //while( ponto!=curva->pontoInicial );

    }

    return 0;
}


/* Funcoes auxiliares ao projecao */
PetscErrorCode CopiaVecUPraMatriz(Vec Vetor, double **Matriz, MALHA M, int rank)
{
    int linha;
//    PetscErrorCode ierr;
    PetscScalar *pVetor;
    //VecScatter scatter;
    //Vec vetorSeq;

//    ierr = VecScatterCreateToZero(Vetor, &scatter, &vetorSeq); CHKERRQ(ierr);
//    ierr = VecScatterBegin(scatter, Vetor, vetorSeq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecScatterEnd(scatter, Vetor, vetorSeq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);




    VecGetArray(Vetor, &pVetor);
    for( linha=0; linha<M.qtdIncognitasU; linha++ ) {
        Matriz[M.indicesU[linha].i][M.indicesU[linha].j] = pVetor[linha];
//        ierr = VecGetValues(Vetor, 1, &linha, &(Matriz[M.indicesU[linha].i][M.indicesU[linha].j])); CHKERRQ(ierr);
    }

    //BroadcastMatrizMPI(Matriz, M.Nx+1, M.Ny);

    //VecDestroy(&vetorSeq); CHKERRQ(ierr);
    VecRestoreArray(Vetor, &pVetor);
    return 0;
}

PetscErrorCode CopiaVecVPraMatriz(Vec Vetor, double **Matriz, MALHA M, int rank)
{
    int linha;
//    PetscErrorCode ierr;
    PetscScalar *pVetor;
    //VecScatter scatter;
    //Vec vetorSeq;

//    ierr = VecScatterCreateToZero(Vetor, &scatter, &vetorSeq); CHKERRQ(ierr);
//    ierr = VecScatterBegin(scatter, Vetor, vetorSeq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecScatterEnd(scatter, Vetor, vetorSeq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);

    VecGetArray(Vetor, &pVetor);
    if( rank==0 ) {
        for( linha=0; linha<M.qtdIncognitasV; linha++ ) {
            Matriz[M.indicesV[linha].i][M.indicesV[linha].j] = pVetor[linha];
//            ierr = VecGetValues(Vetor, 1, &linha, &(Matriz[M.indicesV[linha].i][M.indicesV[linha].j])); CHKERRQ(ierr);
        }
    }

    //BroadcastMatrizMPI(Matriz, M.Nx, M.Ny+1);

    //ierr = VecDestroy(&vetorSeq); CHKERRQ(ierr);
    VecRestoreArray(Vetor, &pVetor);
    return 0;
}

PetscErrorCode CopiaVecPsiPraMatriz(Vec Vetor, double **Matriz, MALHA M, int rank)
{
    int linha;
//    PetscErrorCode ierr;
    PetscScalar *pVetor;
//    VecScatter scatter;
//    Vec vetorSeq;

//    ierr = VecScatterCreateToZero(Vetor, &scatter, &vetorSeq); CHKERRQ(ierr);
//    ierr = VecScatterBegin(scatter, Vetor, vetorSeq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecScatterEnd(scatter, Vetor, vetorSeq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
//    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);

    int sizevec;
    VecGetSize(Vetor, &sizevec);

    VecGetArray(Vetor, &pVetor);
    if( rank==0 ) {
        for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
            Matriz[M.indicesP[linha].i][M.indicesP[linha].j] = pVetor[linha];
//            ierr = VecGetValues(Vetor, 1, &linha, &(Matriz[M.indicesP[linha].i][M.indicesP[linha].j])); CHKERRQ(ierr);
            //Matriz[M.indicesP[linha].i][M.indicesP[linha].j] = vetor[linha];
        }
    }

    //BroadcastMatrizMPI(Matriz, M.Nx, M.Ny+1);

//    ierr = VecDestroy(&vetorSeq); CHKERRQ(ierr);
    VecRestoreArray(Vetor, &pVetor);
    return 0;
}

double NormaEntreSolucoesU(SOLVER_PROJECAO Solver, MALHA M)
{
    int i, j, qtd, linha = 0, linhaNorma = 0;
    double valor;
    PetscErrorCode ierr;
    Vec vetor;

    if(M.tipoCoord==CARTESIANO)
        qtd = (M.Nx+1)*(M.Ny) + (M.Nx)*(M.Ny+1) + 3*(M.Nx)*(M.Ny);
    else
        qtd = (M.Nx+1)*(M.Ny) + (M.Nx)*(M.Ny+1) + 6*(M.Nx)*(M.Ny);

    ierr = VecCreate(PETSC_COMM_WORLD, &vetor); CHKERRQ(ierr);
	ierr = VecSetSizes(vetor, PETSC_DECIDE, qtd); CHKERRQ(ierr);
    ierr = VecSetFromOptions(vetor); CHKERRQ(ierr);


    for( linha=0; linha<M.qtdIncognitasU; linha++ ) {
        i = M.indicesU[linha].i;
        j = M.indicesU[linha].j;

        valor = Solver.UVelho[i][j] - Solver.UNovo[i][j];
        VecSetValues(vetor, 1, &linhaNorma, &valor, INSERT_VALUES);
        
        Solver.residualFaceU[i][j] = fabs(valor);
        linhaNorma++;
    }

    for( linha=0; linha<M.qtdIncognitasV; linha++ ) {
        i = M.indicesV[linha].i;
        j = M.indicesV[linha].j;

        valor = Solver.VVelho[i][j] - Solver.VNovo[i][j];
        VecSetValues(vetor, 1, &linhaNorma, &valor, INSERT_VALUES);

        Solver.residualFaceV[i][j] = fabs(valor);
        linhaNorma++;
    }

    for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
        i = M.indicesP[linha].i;
        j = M.indicesP[linha].j;


        valor = Solver.TxxVelho[i][j] - Solver.TxxNovo[i][j];
        VecSetValues(vetor, 1, &linhaNorma, &valor, INSERT_VALUES);
        linhaNorma++;

        valor = Solver.TxyVelho[i][j] - Solver.TxyNovo[i][j];
        VecSetValues(vetor, 1, &linhaNorma, &valor, INSERT_VALUES);
        linhaNorma++;

        valor = Solver.TyyVelho[i][j] - Solver.TyyNovo[i][j];
        VecSetValues(vetor, 1, &linhaNorma, &valor, INSERT_VALUES);

        valor = Solver.MuNewtGenVelho[i][j] - Solver.MuNewtGenNovo[i][j];
        Solver.residualCell[i][j] = fabs(valor);
        linhaNorma++;
    }

    if( M.tipoCoord==AXI_CILINDRICO ) {

        for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
            i = M.indicesP[linha].i;
            j = M.indicesP[linha].j;

            valor = Solver.WVelho[i][j] - Solver.WNovo[i][j];
            VecSetValues(vetor, 1, &linhaNorma, &valor, INSERT_VALUES);
            linhaNorma++;

            valor = Solver.TxtVelho[i][j] - Solver.TxtNovo[i][j];
            VecSetValues(vetor, 1, &linhaNorma, &valor, INSERT_VALUES);
            linhaNorma++;

            valor = Solver.TttVelho[i][j] - Solver.TttNovo[i][j];
            VecSetValues(vetor, 1, &linhaNorma, &valor, INSERT_VALUES);
            linhaNorma++;
        }

    }


//    VecAssemblyBegin(vetor);
//    VecAssemblyEnd(vetor);


    ierr = VecNorm(vetor, NORM_2, &valor); CHKERRQ(ierr);
    ierr = VecDestroy(&vetor); CHKERRQ(ierr);
    return valor;
}

void CopiaSolucaoNovaParaVelha(SOLVER_PROJECAO *Solver, MALHA M)
{
    int i, j;
    for( i=0; i<=M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            Solver->UVelho[i][j] = Solver->UNovo[i][j];

            // if( M.pontosU[i][j].tipoBoundary==NOSLIP ) {
            //     int iT, jT;
            //     for( iT=0; iT<4; iT++ )
            //         for( jT=0; jT<4; jT++ )
            //             M.pontosU[i][j].tensorParedeVelho[iT][jT] = M.pontosU[i][j].tensorParedeNovo[iT][jT];
            // }
        }
    }


    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<=M.Ny; j++ ) {
            Solver->VVelho[i][j] = Solver->VNovo[i][j];

            // if( M.pontosV[i][j].tipoBoundary==NOSLIP ) {
            //     int iT, jT;
            //     for( iT=0; iT<4; iT++ )
            //         for( jT=0; jT<4; jT++ )
            //             M.pontosV[i][j].tensorParedeVelho[iT][jT] = M.pontosV[i][j].tensorParedeNovo[iT][jT];
            // }
        }
    }


    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            Solver->PVelho[i][j] = Solver->PNovo[i][j];
            Solver->TxxVelho[i][j] = Solver->TxxNovo[i][j];
            Solver->TxyVelho[i][j] = Solver->TxyNovo[i][j];
            Solver->TyyVelho[i][j] = Solver->TyyNovo[i][j];
            Solver->lambdaVelho[i][j] = Solver->lambdaNovo[i][j];
            Solver->muVelho[i][j] = Solver->muNovo[i][j];
            Solver->nuVelho[i][j] = Solver->nuNovo[i][j];

            Solver->evptStructureVelho[i][j] = Solver->evptStructureNovo[i][j];

            Solver->MuNewtGenVelho[i][j] = Solver->MuNewtGenNovo[i][j];

            if( M.tipoCoord == AXI_CILINDRICO ) {
                Solver->WVelho[i][j] = Solver->WNovo[i][j];
                Solver->TxtVelho[i][j] = Solver->TxtNovo[i][j];
                Solver->TytVelho[i][j] = Solver->TytNovo[i][j];
                Solver->TttVelho[i][j] = Solver->TttNovo[i][j];
            }
        }
    }

    return;
}

PetscErrorCode DestroiSistemasLineares(SOLVER_PROJECAO *Solver, MALHA M)
{
    PetscErrorCode ierr;

    Solver->sistemasDestruidos = 1;

    ierr = VecDestroy(&(Solver->vetU)); CHKERRQ(ierr);
    ierr = VecDestroy(&(Solver->solU)); CHKERRQ(ierr);
    ierr = MatDestroy(&(Solver->A_u)); CHKERRQ(ierr);
    ierr = KSPDestroy(&(Solver->kspU)); CHKERRQ(ierr);

    ierr = VecDestroy(&(Solver->vetV)); CHKERRQ(ierr);
    ierr = VecDestroy(&(Solver->solV)); CHKERRQ(ierr);
    ierr = MatDestroy(&(Solver->A_v)); CHKERRQ(ierr);
    ierr = KSPDestroy(&(Solver->kspV)); CHKERRQ(ierr);

   if( (M.tipoCoord==AXI_CILINDRICO) && (M.interface.curvaInicial==NULL)  ) {
       ierr = VecDestroy(&(Solver->vetW)); CHKERRQ(ierr);
       ierr = VecDestroy(&(Solver->solW)); CHKERRQ(ierr);
       ierr = MatDestroy(&(Solver->A_w)); CHKERRQ(ierr);
       ierr = KSPDestroy(&(Solver->kspW)); CHKERRQ(ierr);
   }

    ierr = VecDestroy(&(Solver->vetP)); CHKERRQ(ierr);
    ierr = VecDestroy(&(Solver->solPsi)); CHKERRQ(ierr);
    ierr = MatDestroy(&(Solver->A_p)); CHKERRQ(ierr);
    ierr = KSPDestroy(&(Solver->kspP)); CHKERRQ(ierr);

    return 0;
}

void FinalizaPrograma(MALHA *M, LISTA_CELULAS_SURFACE *ListaSurf, SOLVER_PROJECAO *Solver)
{
    CURVA *curva, *curva2, *proxCurva;
    PONTO *ponto, *proxPonto;
    REGIAO *regiao;
    BLOCO *bloco, *proxBloco;
    CELULA_SURFACE *celulaSurface, *proxCelSurf;

    ///Liberando memoria dos indices do sistema linear
    free(M->indicesU);
    free(M->indicesV);
    free(M->indicesP);
    DesalocaMatriz((void**)M->matIndicesU, M->Nx+1, M->Ny);
    DesalocaMatriz((void**)M->matIndicesV, M->Nx, M->Ny+1);
    DesalocaMatriz((void**)M->matIndicesP, M->Nx, M->Ny);
    DesalocaMatriz((void**)M->matIndicesW, M->Nx, M->Ny);

    //Desaloca pontos bloco paraview
    //DesalocaMatriz((void**)M->PontosBlocoParaview, 5, 3);

    ///informacao das celulas/pontos
    DesalocaMatriz((void**)M->celulas, M->Nx, M->Ny);
    DesalocaMatriz((void**)M->celulasParede, M->Nx, M->Ny);
    DesalocaMatriz((void**)M->pontosU, M->Nx+1, M->Ny);
    DesalocaMatriz((void**)M->pontosV, M->Nx, M->Ny+1);
    DesalocaMatriz((void**)M->pontosP, M->Nx, M->Ny);
    DesalocaMatriz((void**)M->pontosW, M->Nx, M->Ny);
    free(M->dx);
    free(M->dy);
    free(M->r);
    free(M->x);
    free(M->y);

    DesalocaMatriz((void**)Solver->residualFaceU, M->Nx+1, M->Ny);
    DesalocaMatriz((void**)Solver->residualFaceV, M->Nx, M->Ny+1);
    DesalocaMatriz((void**)Solver->residualCell, M->Nx, M->Ny+1);

    DesalocaMatriz((void**)Solver->termoDominanteXX, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->termoDominanteXY,  M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->termoDominanteYY, M->Nx, M->Ny);

    DesalocaMatriz((void**)Solver->MuNewtGenVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->MuNewtGenNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->norm_tau_dev, M->Nx, M->Ny);

    DesalocaMatriz((void**)Solver->evptStructureVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->evptStructureNovo, M->Nx, M->Ny);

    DesalocaMatriz((void**)M->malhaRecursoes, M->Nx, M->Ny);

    if( M->interface.curvaInicial!=NULL ) {
        DesalocaMatriz((void**)M->pontosContinuidade, M->Nx, M->Ny);
        DesalocaMatriz((void**)M->pontosSupLivre, M->Nx+1, M->Ny+1);

        DesalocaMatriz((void**)M->curvaturaCelulas, M->Nx, M->Ny);
        DesalocaMatriz((void**)M->pontoCelula, M->Nx, M->Ny);
        DesalocaMatriz((void**)M->vetNx, M->Nx, M->Ny);
        DesalocaMatriz((void**)M->vetNy, M->Nx, M->Ny);
    }

    ///Blocos da leitura da malha a partir do arquivo .sim
    bloco = M->primBloco;
    while( bloco!=NULL ) {
        proxBloco = bloco->prox;
        free(bloco);
        bloco = proxBloco;
    }

    ///Liberando as estruturas do solver projecao
    DesalocaMatriz((void**)Solver->UVelho, M->Nx+1, M->Ny);
    DesalocaMatriz((void**)Solver->UMeio, M->Nx+1, M->Ny);
    DesalocaMatriz((void**)Solver->UNovo, M->Nx+1, M->Ny);
    DesalocaMatriz((void**)Solver->VVelho, M->Nx, M->Ny+1);
    DesalocaMatriz((void**)Solver->VMeio, M->Nx, M->Ny+1);
    DesalocaMatriz((void**)Solver->VNovo, M->Nx, M->Ny+1);
    DesalocaMatriz((void**)Solver->WVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->WMeio, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->WNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->PVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->Psi, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->PNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TxxVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TxxNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TxyVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TxyNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TyyVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TyyNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TytNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TttNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TxtNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TytVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TttVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->TxtVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->lambdaVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->lambdaNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->lambdaAux, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->muNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->muVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->nuNovo, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->nuVelho, M->Nx, M->Ny);
    DesalocaMatriz((void**)Solver->nuAux, M->Nx, M->Ny);



    if( M->formulacaoVisco==VISCO_LOG ) {
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiVelho_xx, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiVelho_xy, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiVelho_yy, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiVelho_xt, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiVelho_yt, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiVelho_tt, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiNovo_xx, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiNovo_xy, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiNovo_yy, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiNovo_xt, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiNovo_yt, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.PsiNovo_tt, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_11, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_12, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_13, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_21, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_22, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_23, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_31, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_32, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.O_33, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.lambda1, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.lambda2, M->Nx, M->Ny);
        DesalocaMatriz((void**)Solver->variaveis_log_conf.lambda3, M->Nx, M->Ny);
    }

    if( !Solver->sistemasDestruidos )
        DestroiSistemasLineares(Solver, *M);


    ///Deletando as regioes da interface
    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {

        //Regiao da esquerda primeiro
        regiao = curva->regiaoEsq;
        //Coloca NULL em todos os ponteiros (de todas curvas) que apontam pra essa regiao
        if( regiao!=NULL ) {
            for( curva2=M->interface.curvaInicial; curva2!=NULL; curva2=curva2->proxCurva ) {
                if( curva2->regiaoEsq==regiao )
                    curva2->regiaoEsq = NULL;
                if( curva2->regiaoDir==regiao )
                    curva2->regiaoDir = NULL;
            }

            free(regiao);
        }

        //Regiao da direita agora
        regiao = curva->regiaoDir;
        //Coloca NULL em todos os ponteiros (de todas curvas) que apontam pra essa regiao
        if( regiao!=NULL ) {
            for( curva2=M->interface.curvaInicial; curva2!=NULL; curva2=curva2->proxCurva ) {
                if( curva2->regiaoEsq==regiao )
                    curva2->regiaoEsq = NULL;
                if( curva2->regiaoDir==regiao )
                    curva2->regiaoDir = NULL;
            }

            free(regiao);
        }
    }

    ///Deletando toda a interface
    curva = M->interface.curvaInicial;
    while( curva!=NULL ) {
        proxCurva = curva->proxCurva;

        ponto = curva->pontoInicial;
        do {
            proxPonto = ponto->prox;
            free(ponto);
            ponto = proxPonto;
        } while(ponto!=curva->pontoInicial);
        curva->pontoInicial = NULL;

        free(curva);
        curva = proxCurva;
    }


    ///Lista de celulas surface
    celulaSurface = ListaSurf->prim;
    while( celulaSurface!=NULL ) {
        proxCelSurf = celulaSurface->prox;
        free(celulaSurface);
        celulaSurface = proxCelSurf;
    }




    return;
}



/// ==== Calculo de energias na simulacao
double CalculaEnergiaCinetica(MALHA *M, double **U, double **V)
{
    double energia, velU, velV, normaVel;
    int i, j, linha;

    double total_volume = 0.0;

    energia = 0.0;
    /// Percorrendo cada celula e calculando a energia na area dela (depois somando tudo)
    for( linha=0; linha<M->qtdIncognitasP; linha++ ) {
        i = M->indicesP[linha].i;
        j = M->indicesP[linha].j;

        if( M->pontosP[i][j].tipo!=FULL && M->pontosP[i][j].tipo!=SURFACE )
            DeuProblema("Algo deu errado aqui.... CalculaEnergiaCinetica\n");

        /// Vou ignorar celulas SURFACE, pois tambem ignorei na energia de pressao...
        if( M->pontosP[i][j].tipo==SURFACE )
            continue;

        // Interpolando velocidades U e V no centro desta celula
        velU = 0.5*(U[i][j] + U[i+1][j]);
        velV = 0.5*(V[i][j] + V[i][j+1]);
        normaVel = velU*velU + velV*velV; //norma ao quadrado

        // Volume desta celula (ou area no caso 2D)
        double volume;
        if( M->tipoCoord==AXI_CILINDRICO ) 
            volume = M_PI*(M->x[i+1] + M->x[i])*M->dx[i]*M->dy[j];
        else
            volume = M->dx[i]*M->dy[j];

        // Energia cinetica: E_c = 0.5*densidade*v^2
        energia += volume*normaVel;
        total_volume += volume;
    }

//    PrintDebug("Energia cinetica: %.15lf\n", total_volume);

    return 0.5*M->Re*M->weber*energia;
}

double CalculaEnergiaGrav(MALHA *M) {
    double energia, volume;
    int i, j, linha;

    energia = 0.0;
    for( linha=0; linha<M->qtdIncognitasP; linha++ ) {
        i = M->indicesP[linha].i;
        j = M->indicesP[linha].j;

        if( M->pontosP[i][j].tipo!=FULL && M->pontosP[i][j].tipo!=SURFACE )
            DeuProblema("Algo deu errado aqui.... CalculaEnergiaGravitacional\n");

        // Vou ignorar celulas SURFACE pra nao precisar usar contorno pra pressao...
        if( M->pontosP[i][j].tipo==SURFACE )
            continue;

        // Volume desta celula (ou area no caso 2D)
        if( M->tipoCoord==AXI_CILINDRICO ) 
            volume = M_PI*(M->x[i+1] + M->x[i])*M->dx[i]*M->dy[j];
        else
            volume = M->dx[i]*M->dy[j];

        double y = 0.5*( M->y[j] + M->y[j+1] );

        energia += volume * y;
    }

    return ((M->Re * M->weber) / (M->Fr * M->Fr)) * energia;
}

double CalculaEnergiaSuperficie(MALHA *M)
{
    double area, comprimento;
    CURVA *c;
    PONTO *p, *p2;

    /// Calculating the area of the surface (or just the length in 2D)
    area = 0.0;
    for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
        LOOP_CURVA(c) {
            p = loop_stuff.ponto;
            p2 = p->prox;

            comprimento = sqrt( (p->x-p2->x)*(p->x-p2->x) + (p->y-p2->y)*(p->y-p2->y) );
            area += (M->tipoCoord==AXI_CILINDRICO ) ? M_PI*(p->x + p2->x)*comprimento : comprimento;
        }
    }

    // PrintDebug("Energia superficie: %.15lf\n", area);
    return M->Re*area;
}

double CalculaEnergiaDissipativaVisc(MALHA *M, double **U, double **V, double EnergiaAnterior)
{
    double dudx, dudy, dvdx, dvdy;
    double funcao_dissip, vdr;
    int i, j, linha;

    vdr = 0.0;
    for( linha=0; linha<M->qtdIncognitasP; linha++ ) {
        i = M->indicesP[linha].i;
        j = M->indicesP[linha].j;

        if( M->pontosP[i][j].tipo!=FULL && M->pontosP[i][j].tipo!=SURFACE )
            DeuProblema("Algo deu errado aqui.... CalculaEnergiaCinetica\n");

        /// Vou ignorar celulas SURFACE, pois tambem ignorei na energia de pressao...
        if( M->pontosP[i][j].tipo==SURFACE )
            continue;

        dudx = DerivadaDUDX(U, i, j, *M);
        dudy = DerivadaDUDY(U, i, j, 0.0, *M);
        dvdx = DerivadaDVDX(V, i, j, *M);
        dvdy = DerivadaDVDY(V, i, j, *M);

        funcao_dissip = ( 2.0*dudx*dudx + 2.0*dvdy*dvdy + (dvdx + dudy)*(dvdx + dudy) );
        // funcao_dissip += - (2.0/3.0)*( dudx + dvdy );
//        parte_continuidade += -(2.0/3.0)*viscosidade*( dudx + dvdy );

        // Area/Volume desta celula
        double volume;
        if( M->tipoCoord==AXI_CILINDRICO ) 
            volume = M_PI*(M->x[i+1] + M->x[i])*M->dx[i]*M->dy[j];
        else
            volume = M->dx[i]*M->dy[j];

        // vdr: integral da funcao dissipacao no volume do fluido
        vdr += volume*funcao_dissip;
    }

    // Fazendo a integral no tempo (adicionando ao valor que veio do tempo anterior)
    double energia = EnergiaAnterior + (M->dt)*M->weber*vdr;

    // PrintDebug("Energia dissipada: %.15lf\n", energia);

    return energia;
}

void CalculaTodasEnergias(MALHA *M, double **U, double **V, double **P, int PassoTemporal)
{
    static double e_dissipada = 0.0;
    static double e_total = 0.0;
    static int print_energia = 1;

    e_dissipada = CalculaEnergiaDissipativaVisc(M, U, V, e_dissipada);

    static FILE *arq = NULL;
    if( arq==NULL ) {
        char nomeArq[500];
        sprintf(nomeArq, "ArquivosPlot/VTK/%s/energias.txt", M->nomeArquivo);
        arq = fopen(nomeArq, "wt");
    }

    if( (PassoTemporal-1)%print_energia==0 ) {
        double e_cinetica = CalculaEnergiaCinetica(M, U, V);
        double e_sup = CalculaEnergiaSuperficie(M);
        double e_grav = CalculaEnergiaGrav(M);

        //e_total = e_cinetica + e_sup + e_dissipada;
        e_total = e_cinetica + e_grav + e_sup + e_dissipada;
        fprintf(arq, "%lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", PassoTemporal*M->dt,
                e_cinetica, e_grav, e_sup, e_dissipada, e_total);
        fflush(arq);
    }

    return;
}

double CalculaFlowType(MALHA *M, double **U, double **V, int i, int j) {
    double dudx, dudy, dvdx, dvdy;
    double D, omega, flow_type;

    dudx = DerivadaDUDX(U, i, j, *M);
    dudy = DerivadaDUDY(U, i, j, 0.0, *M);
    dvdx = DerivadaDVDX(V, i, j, *M);
    dvdy = DerivadaDVDY(V, i, j, *M);

    D = (1.0/2.0) * sqrt(4.0*dudx*dudx + 2.0 * (dudy + dvdx) * (dudy + dvdx) + 4.0 * dvdy*dvdy);
    omega = (1.0/2.0) * sqrt((dudy - dvdx) * (dudy - dvdx) + (dvdx - dudy) * (dvdx - dudy));

    flow_type = (D - omega) / (D + omega);

    return flow_type;
}

double CalculaMagnitudeVelocidade(double **U, double **V, int i, int j) {
    return sqrt((U[i][j] * U[i][j]) + (V[i][j] * V[i][j]));
}

double CalculaEnergiaDissipativaViscCelula(MALHA *M, double **U, double **V, int i, int j, double EnergiaAnterior) {
    double dudx, dudy, dvdx, dvdy;
    double funcao_dissip, vdr, energia;

    dudx = DerivadaDUDX(U, i, j, *M);
    dudy = DerivadaDUDY(U, i, j, 0.0, *M);
    dvdx = DerivadaDVDX(V, i, j, *M);
    dvdy = DerivadaDVDY(V, i, j, *M);

    funcao_dissip = ( 2.0*dudx*dudx + 2.0*dvdy*dvdy + (dvdx + dudy)*(dvdx + dudy) );

    vdr = M_PI*(M->x[i+1] + M->x[i])*M->dx[i]*M->dy[j] * funcao_dissip;

    energia = EnergiaAnterior + (M->dt)*M->weber*vdr;

    return energia;
}
