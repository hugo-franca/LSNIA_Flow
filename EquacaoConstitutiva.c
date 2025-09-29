#include "EquacaoConstitutiva.h"

extern double **debugsCell;

/// Function used for Generalized Newtonian models to calculate the apparent viscosity in the entire mesh
/// For traditional Newtonian models, this viscosity should always be 1 everywhere
double GeneralizedNewtonian_Viscosity(MALHA M, int n, double **U, double **V, double **MuNewtGen)
{
    int linha;


    for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
        int i = M.indicesP[linha].i;
        int j = M.indicesP[linha].j;

        double dudx = DerivadaDUDX(U, i, j, M);
        double dudy = DerivadaDUDY(U, i, j, (n+1)*M.dt, M);
        double dvdx = DerivadaDVDX(V, i, j, M);
        double dvdy = DerivadaDVDY(V, i, j, M);

        // double strain_rate = sqrt( 2.0*dudx*dudx + 2.0*dvdy*dvdy + 2.0*0.5*(dudy + dvdx)*(dudy + dvdx) );
        double strain_rate = sqrt( dudx*dudx + dvdy*dvdy + 0.5*(dudy + dvdx)*(dudy + dvdx) );

        // double mumax = 10000, muTemp = 1.0;
        // if (strain_rate > 0. && M.Bi > 0.) {
        //     double temp = M.Bi/(sqrt(2.0)*strain_rate) + 1.0;
        //     muTemp = (temp<mumax) ? temp : mumax; //min(temp, mumax)
        // } 
        // else {
        //     if (M.Bi > 0.)
        //         muTemp = mumax;
        //     else
        //         muTemp = 1.0;
        // }
        // MuNewtGen[i][j] = muTemp;


        // Bingham viscoplastic
        double regularization = 1.0e-03;
        double yieldstress_viscosity = M.Bi*(1.0 - exp(-strain_rate/regularization))/(sqrt(2.0)*strain_rate + 1e-20);
        MuNewtGen[i][j] = 1.0 + yieldstress_viscosity;

        /// Basilisk
        // double mumax = 10000.0;
        // if(strain_rate > 0. && M.Bi > 0.) {
        //     double temp = M.Bi/(sqrt(2.0)*strain_rate) + 1.0;
        //     MuNewtGen[i][j] = (temp<mumax) ? temp : mumax;
        // }
        // else {
        //     if (M.Bi > 0.)
        //         MuNewtGen[i][j] = mumax;
        //     else
        //         MuNewtGen[i][j] = 1.0;
        // }



        if( debugsCell==NULL ) {
            double valorInicial = 0.0;
            debugsCell = (double **)AlocaMatriz(M.Nx, M.Ny, sizeof(double), &valorInicial);
        }
        debugsCell[i][j] = dudy;
        

        // if( MuNewtGen[i][j] )

        /// === Cross model
//        MuNewtGen[i][j] = mu_inf + ( mu_zero - mu_inf )/( 1.0 + pow(K*Gamma[i][j], m) );

        /// === Carreau-Yasuda model
        // double expoente = (m-1)/a;
        // double fator = 1.0 + pow(K*Gamma[i][j], a);
        // MuNewtGen[i][j] = mu_inf + ( mu_zero - mu_inf )*pow(fator, expoente);

        /// === Power law (paper kai-san)
//        MuNewtGen[i][j] = mu_zero/( 1.0 + (mu_zero*pow(Gamma[i][j], 1-m)/K) );

        /// ==== Adimensionalizando
//        MuNewtGen[i][j] /= mu_zero;
        // MuNewtGen[i][j] /= mu_apparent;

//        MuNewtGen[i][j] = 1.0;
    }

    return 1;
}


void GeneralizedNewtonian_Tensor(MALHA M, int n, double **U, double **V, double **MuNewtGen,
                                double **Txx, double **Txy, double **Tyy)
{
    int linha;

    for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
        int i = M.indicesP[linha].i;
        int j = M.indicesP[linha].j;

        double dudx = DerivadaDUDX(U, i, j, M);
        double dudy = DerivadaDUDY(U, i, j, (n+1)*M.dt, M);
        double dvdx = DerivadaDVDX(V, i, j, M);
        double dvdy = DerivadaDVDY(V, i, j, M);

        // double strain_rate = sqrt( 2.0*dudx*dudx + 2.0*dvdy*dvdy + 2.0*0.5*(dudy + dvdx)*(dudy + dvdx) );
        double strain_rate = sqrt( dudx*dudx + dvdy*dvdy + 0.5*(dudy + dvdx)*(dudy + dvdx) );

        // Bingham viscoplastic
        double regularization = 1.0e-03;
        double yieldstress_viscosity = M.Bi*(1.0 - exp(-strain_rate/regularization))/(sqrt(2.0)*strain_rate + 1e-20);
        MuNewtGen[i][j] = 1.0 + yieldstress_viscosity; // Calculating the apparent viscosity just for visualizations

        // Non-Newtonian part of the extra-stress tensor in the VP model
        Txx[i][j] = 2.0*yieldstress_viscosity*dudx;
        Txy[i][j] = yieldstress_viscosity*(dudy + dvdx);
        Tyy[i][j] = 2.0*yieldstress_viscosity*dvdy;
        // if( (M.passoTemporal%10==0) && (i==M.Nx/2) ) {
        //     PrintDebug("AAAA %d %d: %20lf %20e %20lf\n", i, j, yieldstress_viscosity, strain_rate, dudy);
        // }
    }

    return;
}

void EquacoesConstitutivas(double **LambdaVelho, double **MuVelho, double **NuVelho, double **LambdaNovo,  double **MuNovo, double **NuNovo, double **LambdaAux, double **NuAux,
                            double **TxxVelho, double **TxyVelho, double **TyyVelho, double **TxtVelho, double **TytVelho, double **TttVelho,
                            double **TxxNovo, double **TxyNovo, double **TyyNovo, double **TxtNovo, double **TytNovo, double **TttNovo,
                            double **EVPT_Structure, double **Norm_tau_dev,
                            VARIAVEIS_LOG_CONFORMATION *Variaveis_log_conf,
                            double **UVelho, double **VVelho, double **WVelho, double **UNovo, double **VNovo, double **WNovo, double **P,
                            double **TermoDominanteXX, double **TermoDominanteXY, double **TermoDominanteYY,
                            int n, MALHA M)
{
    int i, j, linha;
    double xi, tol, velU, velV, modV;

    double **PsiVelho_xx = Variaveis_log_conf->PsiVelho_xx;
    double **PsiVelho_xy = Variaveis_log_conf->PsiVelho_xy;
    double **PsiVelho_yy = Variaveis_log_conf->PsiVelho_yy;
    double **PsiVelho_xt = Variaveis_log_conf->PsiVelho_xt;
    double **PsiVelho_yt = Variaveis_log_conf->PsiVelho_yt;
    double **PsiVelho_tt = Variaveis_log_conf->PsiVelho_tt;
    double **PsiNovo_xx = Variaveis_log_conf->PsiNovo_xx;
    double **PsiNovo_xy = Variaveis_log_conf->PsiNovo_xy;
    double **PsiNovo_yy = Variaveis_log_conf->PsiNovo_yy;
    double **PsiNovo_xt = Variaveis_log_conf->PsiNovo_xt;
    double **PsiNovo_yt = Variaveis_log_conf->PsiNovo_yt;
    double **PsiNovo_tt = Variaveis_log_conf->PsiNovo_tt;
    double **O_11 = Variaveis_log_conf->O_11;
    double **O_12 = Variaveis_log_conf->O_12;
    double **O_13 = Variaveis_log_conf->O_13;
    double **O_21 = Variaveis_log_conf->O_21;
    double **O_22 = Variaveis_log_conf->O_22;
    double **O_23 = Variaveis_log_conf->O_23;
    double **O_31 = Variaveis_log_conf->O_31;
    double **O_32 = Variaveis_log_conf->O_32;
    double **O_33 = Variaveis_log_conf->O_33;
    double **Lambda1 = Variaveis_log_conf->lambda1;
    double **Lambda2 = Variaveis_log_conf->lambda2;
    double **Lambda3 = Variaveis_log_conf->lambda3;

    // VP model: no need for additional constitutive equations
    if( M.tipo_modelo==MODELO_VP )
        return;

    // No need to solve anything if beta==1
    if( M.beta==1.0 )
        return;

    xi = (1.0 - M.beta)/(/*M.Re**/M.We);
    tol = M.tolNSF;

    /// Fazendo uma passagem inicial pra atualizar as condicoes de contorno para o tensor
//    for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
//        break;
//        i = M.indicesP[linha].i;
//        j = M.indicesP[linha].j;
//
//        double **tVelho, **tNovo;
//        double velU, velV, velW, coordR;
//        double dwdr = 0.0, dwdz = 0.0, dvdr = 0.0, dvdz = 0.0, dudr = 0.0, dudz = 0.0;
//        double convT = 0.0;
//        double giesekus;
//        int parede = 0;
//
//
//        if( M.pontosU[i][j].tipoBoundary==NOSLIP || M.pontosU[i+1][j].tipoBoundary==NOSLIP ) {
//            int iParede = (M.pontosU[i][j].tipoBoundary==NOSLIP) ? i : i+1;
//            parede = 1;
//
//            tVelho = M.pontosU[iParede][j].tensorParedeVelho;
//            tNovo = M.pontosU[iParede][j].tensorParedeNovo;
//            velU = 0.0;
//            velV = 0.0;
//            velW = M.pontosU[iParede][j].valorDirichlet;
//            coordR = M.x[iParede];
//
//            if( M.pontosU[iParede][j].extremo==ESQUERDA ) {
//                dudr = (UNovo[i+1][j] - UNovo[i][j])/M.dx[i];
//                dvdr = (0.5*(VNovo[i][j] + VNovo[i][j+1]) - velV)/( 0.5*M.dx[i] );
//                dwdr = (WNovo[i][j] - velW)/( 0.5*M.dx[i] );
//            }
//            else if( M.pontosU[iParede][j].extremo==DIREITA ) {
//                dudr = (UNovo[i+1][j] - UNovo[i][j])/M.dx[i];
//                dvdr = (velV - 0.5*(VNovo[i][j] + VNovo[i][j+1]))/( 0.5*M.dx[i] );
//                dwdr = (velW - WNovo[i][j])/( 0.5*M.dx[i] );
//            }
//            else
//                DeuProblema("EquacoesConstitutivas: problema na atualizacao do contorno.\n\n");
//            dudz = dwdz = dvdz = 0.0;
//
//            convT = - (2.0*tVelho[1][3])*velW/coordR;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][1]*tVelho[1][1] + tVelho[1][2]*tVelho[1][2] + tVelho[1][3]*tVelho[1][3]);
//            tNovo[1][1] = tVelho[1][1] + M.dt*( - convT + 2.0*dudr*tVelho[1][1] + 2.0*dudz*tVelho[1][2] - 2.0*velW*tVelho[1][3]/coordR + 2.0*xi*dudr - (tVelho[1][1]/M.We) - giesekus );
//
//            convT = - (tVelho[2][3])*velW/coordR;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][1]*tVelho[1][2] + tVelho[1][2]*tVelho[2][2] + tVelho[1][3]*tVelho[2][3]);
//            tNovo[1][2] = tVelho[1][2] + M.dt*( - convT + dudr*tVelho[1][2] - velW*tVelho[2][3]/coordR + dudz*tVelho[2][2] + dvdr*tVelho[1][1] + dvdz*tVelho[1][2] + xi*(dudz + dvdr) - (tVelho[1][2]/M.We) - giesekus );
//
//            convT = (tVelho[1][1] - tVelho[3][3])*velW/M.r[i];
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][1]*tVelho[1][3] + tVelho[1][2]*tVelho[2][3] + tVelho[1][3]*tVelho[3][3]);
//            tNovo[1][3] = tVelho[1][3] + M.dt*( - convT + dudr*tVelho[1][3] - velW*tVelho[3][3]/coordR + dudz*tVelho[2][3] + dwdr*tVelho[1][1] + velU*tVelho[1][3]/coordR + dwdz*tVelho[1][2] + xi*(dwdr - velW/coordR) - (tVelho[1][3]/M.We) - giesekus );
//
//            convT = 0.0;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][2]*tVelho[1][2] + tVelho[2][2]*tVelho[2][2] + tVelho[2][3]*tVelho[2][3]);
//            tNovo[2][2] = tVelho[2][2] + M.dt*( - convT + 2.0*dvdr*tVelho[1][2] + 2.0*dvdz*tVelho[2][2] + 2.0*xi*dvdz - (tVelho[2][2]/M.We) - giesekus);
//
//            convT = (tVelho[1][2])*velW/coordR;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][2]*tVelho[1][3] + tVelho[2][2]*tVelho[2][3] + tVelho[2][3]*tVelho[3][3]);
//            tNovo[2][3] = tVelho[2][3] + M.dt*( - convT + dwdr*tVelho[1][2] + velU*tVelho[2][3]/coordR + dwdz*tVelho[2][2] + dvdr*tVelho[1][3] + dvdz*tVelho[2][3] + xi*dwdz - (tVelho[2][3]/M.We) - giesekus );
//
//            convT = (2.0*tVelho[1][1])*velW/coordR;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][3]*tVelho[1][3] + tVelho[2][3]*tVelho[2][3] + tVelho[3][3]*tVelho[3][3]);
//            tNovo[3][3] = tVelho[3][3] + M.dt*( - convT + 2.0*dwdr*tVelho[1][3] + 2.0*velU*tVelho[3][3]/coordR + 2.0*dwdz*tVelho[2][3] + 2.0*xi*velU/coordR - (tVelho[3][3]/M.We) - giesekus );
//
//            tNovo[2][1] = tNovo[1][2];
//            tNovo[3][1] = tNovo[1][3];
//            tNovo[3][2] = tNovo[2][3];
//        }
//
//        if( M.pontosV[i][j].tipoBoundary==NOSLIP || M.pontosV[i][j+1].tipoBoundary==NOSLIP ) {
//            int jParede = (M.pontosV[i][j].tipoBoundary==NOSLIP) ? j : j+1;
//            parede = 1;
//
//            tVelho = M.pontosV[i][jParede].tensorParedeVelho;
//            tNovo = M.pontosV[i][jParede].tensorParedeNovo;
//            velU = 0.0;
//            velV = 0.0;
//            velW = M.pontosV[i][jParede].valorDirichlet;
//            coordR = 0.5*( M.x[i] + M.x[i+1]  );
//
//            if( M.pontosV[i][jParede].extremo==BAIXO ) {
//                dudz = (0.5*(UNovo[i][j] + UNovo[i+1][j]) - velU)/( 0.5*M.dy[j] );
//                dvdz = (VNovo[i][j+1] - VNovo[i][j])/M.dy[j];
//                dwdz = (WNovo[i][j] - velW)/( 0.5*M.dy[j] );
//            }
//            else if( M.pontosV[i][jParede].extremo==CIMA ) {
//                dudz = (velU - 0.5*(UNovo[i][j] + UNovo[i+1][j]))/( 0.5*M.dy[j] );
//                dvdz = (VNovo[i][j+1] - VNovo[i][j])/M.dy[j];
//                dwdz = (velW - WNovo[i][j])/( 0.5*M.dy[j] );
//            }
//            else
//                DeuProblema("EquacoesConstitutivas: problema na atualizacao do contorno.\n\n");
//            dudr = dwdr = dvdr = 0.0;
//
//            convT = - (2.0*tVelho[1][3])*velW/coordR;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][1]*tVelho[1][1] + tVelho[1][2]*tVelho[1][2] + tVelho[1][3]*tVelho[1][3]);
//            tNovo[1][1] = tVelho[1][1] + M.dt*( - convT + 2.0*dudr*tVelho[1][1] + 2.0*dudz*tVelho[1][2] - 2.0*velW*tVelho[1][3]/coordR + 2.0*xi*dudr - (tVelho[1][1]/M.We) - giesekus );
//
//            convT = - (tVelho[2][3])*velW/coordR;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][1]*tVelho[1][2] + tVelho[1][2]*tVelho[2][2] + tVelho[1][3]*tVelho[2][3]);
//            tNovo[1][2] = tVelho[1][2] + M.dt*( - convT + dudr*tVelho[1][2] - velW*tVelho[2][3]/coordR + dudz*tVelho[2][2] + dvdr*tVelho[1][1] + dvdz*tVelho[1][2] + xi*(dudz + dvdr) - (tVelho[1][2]/M.We) - giesekus );
//
//            convT = (tVelho[1][1] - tVelho[3][3])*velW/M.r[i];
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][1]*tVelho[1][3] + tVelho[1][2]*tVelho[2][3] + tVelho[1][3]*tVelho[3][3]);
//            tNovo[1][3] = tVelho[1][3] + M.dt*( - convT + dudr*tVelho[1][3] - velW*tVelho[3][3]/coordR + dudz*tVelho[2][3] + dwdr*tVelho[1][1] + velU*tVelho[1][3]/coordR + dwdz*tVelho[1][2] + xi*(dwdr - velW/coordR) - (tVelho[1][3]/M.We) - giesekus );
//
//            convT = 0.0;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][2]*tVelho[1][2] + tVelho[2][2]*tVelho[2][2] + tVelho[2][3]*tVelho[2][3]);
//            tNovo[2][2] = tVelho[2][2] + M.dt*( - convT + 2.0*dvdr*tVelho[1][2] + 2.0*dvdz*tVelho[2][2] + 2.0*xi*dvdz - (tVelho[2][2]/M.We) - giesekus);
//
//            convT = (tVelho[1][2])*velW/coordR;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][2]*tVelho[1][3] + tVelho[2][2]*tVelho[2][3] + tVelho[2][3]*tVelho[3][3]);
//            tNovo[2][3] = tVelho[2][3] + M.dt*( - convT + dwdr*tVelho[1][2] + velU*tVelho[2][3]/coordR + dwdz*tVelho[2][2] + dvdr*tVelho[1][3] + dvdz*tVelho[2][3] + xi*dwdz - (tVelho[2][3]/M.We) - giesekus );
//
//            convT = (2.0*tVelho[1][1])*velW/coordR;
//            giesekus = (M.alpha/(1.0-M.beta))*(tVelho[1][3]*tVelho[1][3] + tVelho[2][3]*tVelho[2][3] + tVelho[3][3]*tVelho[3][3]);
//            tNovo[3][3] = tVelho[3][3] + M.dt*( - convT + 2.0*dwdr*tVelho[1][3] + 2.0*velU*tVelho[3][3]/coordR + 2.0*dwdz*tVelho[2][3] + 2.0*xi*velU/coordR - (tVelho[3][3]/M.We) - giesekus );
//
//            tNovo[2][1] = tNovo[1][2];
//            tNovo[3][1] = tNovo[1][3];
//            tNovo[3][2] = tNovo[2][3];
//        }

//        if( M.pontosU[i][j].tipoBoundary==NOSLIP && j==M.Ny/2 && (M.passoTemporal%M.intervaloPrint==0) ) {
//            int iT, jT;
//            PrintDebug("VELW: %lf\n", velW);
//            for( iT=1; iT<4; iT++ ) {
//                for( jT=1; jT<4; jT++ ) {
//                    PrintDebug("%e ", 0.5*M.pontosU[i][j].tensorParedeVelho[iT][jT]);
//                }
//                PrintDebug("\n");
//            }
//            PrintDebug("\n");
//        }

//    }




    ///Calculando as variaveis auxiliares pro termo convectivo do NSF
    if( M.formulacaoVisco==VISCO_NSF || M.formulacaoVisco==VISCO_HIBRIDO ) {
        for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
            i = M.indicesP[linha].i;
            j = M.indicesP[linha].j;

            velU = 0.5*( UNovo[i+1][j] + UNovo[i][j] );
            velV = 0.5*( VNovo[i][j+1] + VNovo[i][j] );
            //velU = 0.5*( UVelho[i+1][j] + UVelho[i][j] );
            //velV = 0.5*( VVelho[i][j+1] + VVelho[i][j] );
            modV = velU*velU + velV*velV;

            LambdaAux[i][j] = LambdaVelho[i][j]/( modV + tol );
            NuAux[i][j] = modV*NuVelho[i][j];
        }
    }

    if( M.formulacaoVisco==VISCO_LOG ) {
        for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
            i = M.indicesP[linha].i;
            j = M.indicesP[linha].j;

            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            A[1][1] = 1.0 + TxxVelho[i][j]/xi;
            A[1][2] = TxyVelho[i][j]/xi;
            A[2][1] = TxyVelho[i][j]/xi;
            A[2][2] = 1.0 + TyyVelho[i][j]/xi;
            if( M.tipoCoord==AXI_CILINDRICO ) {
                A[1][3] = A[3][1] = TxtVelho[i][j]/xi;
                A[2][3] = A[3][2] = TytVelho[i][j]/xi;
                A[3][3] = 1.0 + TttVelho[i][j]/xi;
            }

            CalculaAutovalores_Jacobi(A, (M.tipoCoord==AXI_CILINDRICO) ? 3 : 2, auto_valores, aux1, aux2, O, &nrot);

            PsiVelho_xx[i][j] = O[1][1]*O[1][1]*log(auto_valores[1]) + O[1][2]*O[1][2]*log(auto_valores[2]);

            PsiVelho_xy[i][j] = O[1][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[1][2]*log(auto_valores[2]);

            PsiVelho_yy[i][j] = O[2][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[2][2]*log(auto_valores[2]);
            if( M.tipoCoord==AXI_CILINDRICO ) {
                PsiVelho_xx[i][j] += O[1][3]*O[1][3]*log(auto_valores[3]);
                PsiVelho_xy[i][j] += O[1][3]*O[2][3]*log(auto_valores[3]);
                PsiVelho_yy[i][j] += O[2][3]*O[2][3]*log(auto_valores[3]);

                PsiVelho_xt[i][j] = O[1][1]*O[3][1]*log(auto_valores[1]) + O[1][2]*O[3][2]*log(auto_valores[2]) + O[1][3]*O[3][3]*log(auto_valores[3]);
                PsiVelho_yt[i][j] = O[2][1]*O[3][1]*log(auto_valores[1]) + O[2][2]*O[3][2]*log(auto_valores[2]) + O[2][3]*O[3][3]*log(auto_valores[3]);
                PsiVelho_tt[i][j] = O[3][1]*O[3][1]*log(auto_valores[1]) + O[3][2]*O[3][2]*log(auto_valores[2]) + O[3][3]*O[3][3]*log(auto_valores[3]);
            }


            O_11[i][j] = O[1][1];
            O_12[i][j] = O[1][2];
            O_13[i][j] = O[1][3];
            O_21[i][j] = O[2][1];
            O_22[i][j] = O[2][2];
            O_23[i][j] = O[2][3];
            O_31[i][j] = O[3][1];
            O_32[i][j] = O[3][2];
            O_33[i][j] = O[3][3];

            Lambda1[i][j] = auto_valores[1];
            Lambda2[i][j] = auto_valores[2];
            Lambda3[i][j] = auto_valores[3];
        }
    }

    ///Agora aplicando as equacoes constitutivas em cada celula
    for( linha=0; linha<M.qtdIncognitasP; linha++ ) {
        i = M.indicesP[linha].i;
        j = M.indicesP[linha].j;

        velU = 0.5*( UNovo[i+1][j] + UNovo[i][j] );
        velV = 0.5*( VNovo[i][j+1] + VNovo[i][j] );
        //velU = 0.5*( UVelho[i+1][j] + UVelho[i][j] );
        //velV = 0.5*( VVelho[i][j+1] + VVelho[i][j] );
        modV = velU*velU + velV*velV;

        if( M.formulacaoVisco==VISCO_CART )
            EquacoesConstitutivasCartesiano(i, j, TxxVelho, TxyVelho, TyyVelho, TxtVelho, TytVelho, TttVelho,
                                            TxxNovo, TxyNovo, TyyNovo, TxtNovo, TytNovo, TttNovo,
                                            EVPT_Structure, Norm_tau_dev,
                                            UNovo, VNovo, WNovo, P, n, M);
        else if( M.formulacaoVisco==VISCO_NSF ) {
            EquacoesConstitutivasNSF(i, j, LambdaVelho, MuVelho, NuVelho, LambdaNovo, MuNovo, NuNovo, LambdaAux, NuAux,
                                        UVelho, VVelho, UNovo, VNovo, n, M);

            //CONVERSAO
            tol = M.tolNSFconversao;
            TxxNovo[i][j] = xi*( -1.0 + ( ( LambdaNovo[i][j]*velU*velU - MuNovo[i][j]*2.0*velU*velV + NuNovo[i][j]*velV*velV )/( modV + tol ) ) );
            TxyNovo[i][j] = (xi/(modV + tol))*( LambdaNovo[i][j]*velU*velV + MuNovo[i][j]*(velU*velU - velV*velV) - NuNovo[i][j]*velU*velV );
            TyyNovo[i][j] = xi*( -1.0 + ( ( LambdaNovo[i][j]*velV*velV + MuNovo[i][j]*2.0*velU*velV + NuNovo[i][j]*velU*velU )/( modV + tol ) ) );
        }
        else if( M.formulacaoVisco==VISCO_HIBRIDO ) {

            if( ProximoDaParede_Hibrido(i, j, M) ) {

                //PrintDebug("(%d %d) ", i, j);

                EquacoesConstitutivasCartesiano(i, j, TxxVelho, TxyVelho, TyyVelho, TxtVelho, TytVelho, TttVelho,
                                            TxxNovo, TxyNovo, TyyNovo, TxtNovo, TytNovo, TttNovo,
                                            EVPT_Structure, Norm_tau_dev,
                                            UNovo, VNovo, WNovo, P, n, M);
                //EquacoesConstitutivasCartesiano(i, j, TxxVelho, TxyVelho, TyyVelho, TxxNovo, TxyNovo, TyyNovo, UVelho, VVelho, n, M);

                //Conversao
                tol = M.tolNSFconversao;
                LambdaNovo[i][j] = 1.0 + ( (TxxNovo[i][j]*velU*velU + TxyNovo[i][j]*2.0*velU*velV + TyyNovo[i][j]*velV*velV)/(xi*(modV+tol)) );
                MuNovo[i][j] = (-TxxNovo[i][j]*velU*velV + TxyNovo[i][j]*(velU*velU - velV*velV) + TyyNovo[i][j]*velU*velV)/(xi*(modV+tol));
                NuNovo[i][j] = 1.0 + ( (TxxNovo[i][j]*velV*velV - TxyNovo[i][j]*2.0*velU*velV + TyyNovo[i][j]*velU*velU)/(xi*(modV+tol)) );
            }
            else {

                EquacoesConstitutivasNSF(i, j, LambdaVelho, MuVelho, NuVelho, LambdaNovo, MuNovo, NuNovo, LambdaAux, NuAux,
                                        UVelho, VVelho, UNovo, VNovo, n, M);

                //Conversao
                tol = M.tolNSFconversao;
                TxxNovo[i][j] = xi*( -1.0 + ( ( LambdaNovo[i][j]*velU*velU - MuNovo[i][j]*2.0*velU*velV + NuNovo[i][j]*velV*velV )/( modV + tol ) ) );
                TxyNovo[i][j] = (xi/(modV + tol))*( LambdaNovo[i][j]*velU*velV + MuNovo[i][j]*(velU*velU - velV*velV) - NuNovo[i][j]*velU*velV );
                TyyNovo[i][j] = xi*( -1.0 + ( ( LambdaNovo[i][j]*velV*velV + MuNovo[i][j]*2.0*velU*velV + NuNovo[i][j]*velU*velU )/( modV + tol ) ) );
            }

        }
        else if( M.formulacaoVisco==VISCO_LOG ) {
            EquacoesConstitutivasLogConformation(i, j,
                                                 PsiVelho_xx, PsiVelho_xy, PsiVelho_yy, PsiVelho_xt, PsiVelho_yt, PsiVelho_tt,
                                                 PsiNovo_xx, PsiNovo_xy, PsiNovo_yy, PsiNovo_xt, PsiNovo_yt, PsiNovo_tt,
                                                 O_11, O_12, O_13, O_21, O_22, O_23, O_31, O_32, O_33,
                                                 Lambda1, Lambda2, Lambda3,
                                                 UNovo, VNovo, WNovo, n, M);

            double psi[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            psi[1][1] = PsiNovo_xx[i][j];
            psi[1][2] = PsiNovo_xy[i][j];
            psi[2][1] = PsiNovo_xy[i][j];
            psi[2][2] = PsiNovo_yy[i][j];
            if( M.tipoCoord==AXI_CILINDRICO ) {
                psi[1][3] = PsiNovo_xt[i][j];
                psi[2][3] = PsiNovo_yt[i][j];
                psi[3][1] = PsiNovo_xt[i][j];
                psi[3][2] = PsiNovo_yt[i][j];
                psi[3][3] = PsiNovo_tt[i][j];
            }

            CalculaAutovalores_Jacobi(psi, (M.tipoCoord==AXI_CILINDRICO)? 3 : 2, auto_valores, aux1, aux2, O, &nrot);

            double A11, A12, A13, A22, A23, A33;
            A11 = O[1][1]*O[1][1]*exp(auto_valores[1]) + O[1][2]*O[1][2]*exp(auto_valores[2]);
            A12 = O[1][1]*O[2][1]*exp(auto_valores[1]) + O[2][2]*O[1][2]*exp(auto_valores[2]);
            A22 = O[2][1]*O[2][1]*exp(auto_valores[1]) + O[2][2]*O[2][2]*exp(auto_valores[2]);
            if( M.tipoCoord==AXI_CILINDRICO ) {
                A11 += O[1][3]*O[1][3]*exp(auto_valores[3]);
                A12 += O[1][3]*O[2][3]*exp(auto_valores[3]);
                A22 += O[2][3]*O[2][3]*exp(auto_valores[3]);

                A13 = O[1][1]*O[3][1]*exp(auto_valores[1]) + O[1][2]*O[3][2]*exp(auto_valores[2]) + O[1][3]*O[3][3]*exp(auto_valores[3]);
                A23 = O[2][1]*O[3][1]*exp(auto_valores[1]) + O[2][2]*O[3][2]*exp(auto_valores[2]) + O[2][3]*O[3][3]*exp(auto_valores[3]);
                A33 = O[3][1]*O[3][1]*exp(auto_valores[1]) + O[3][2]*O[3][2]*exp(auto_valores[2]) + O[3][3]*O[3][3]*exp(auto_valores[3]);
            }


            TxxNovo[i][j] = xi*(A11 - 1.0);
            TxyNovo[i][j] = xi*(A12);
            TyyNovo[i][j] = xi*(A22 - 1.0);
            if( M.tipoCoord==AXI_CILINDRICO ) {
                TxtNovo[i][j] = xi*(A13);
                TytNovo[i][j] = xi*(A23);
                TttNovo[i][j] = xi*(A33 - 1.0);
            }
        }
        else DeuProblema("EquacoesConstitutivas: Opcao viscoelastico invalida.\n\n");

        CalculaTermoDominanteEqsConstituvas(i, j, TxxVelho, TxyVelho, TyyVelho, TxtVelho, TytVelho, TttVelho,
                                            TxxNovo, TxyNovo, TyyNovo, TxtNovo, TytNovo, TttNovo,
                                            UNovo, VNovo, WNovo, TermoDominanteXX, TermoDominanteXY, TermoDominanteYY,
                                            n, M);


    }

    //DeuProblema("ACABOU\n");

}

void EquacoesConstitutivasCartesiano(int i, int j,
                                     double **Txx_Velho, double **Txy_Velho, double **Tyy_Velho, double **Txt_Velho, double **Tyt_Velho, double **Ttt_Velho,
                                     double **Txx_Novo, double **Txy_Novo, double **Tyy_Novo, double **Txt_Novo, double **Tyt_Novo, double **Ttt_Novo,
                                     double **EVPT_Structure, double **Norm_tau_dev,
                                     double **U, double **V, double **W, double **P, int n, MALHA M)
{
    double dudx, dudy, dvdx, dvdy, dwdx = 0.0, dwdy = 0.0;
    double convecTxx, convecTxy, convecTyy, convecTtt, convecTxt, convecTyt;
    double xi;
    double parteWi, parteCilindrico;
    double evpt_elasticModulus, evpt_viscosity;

    // Ignora se for beta==1
    if( M.beta==1.0 )
        return;

    dudx = DerivadaDUDX(U, i, j, M);
    dudy = DerivadaDUDY(U, i, j, (n+1)*M.dt, M);
    dvdx = DerivadaDVDX(V, i, j, M);
    dvdy = DerivadaDVDY(V, i, j, M);

    /// In the case Wi = 0, we can completely ignore the Oldroyd Derivative
    /// We only need to solve the Bingham model ((|t_d| - Bi)/|td|)*T = 2(1-beta)D
    // if( Wi==0 ) {
    //     double stress_dev = EVPT_CalculaSigmaDev(M, i, j, U, V, P, Txx_Velho, Txy_Velho, Tyy_Velho);
    //     double bingham_criteria = (stress_dev - M.Bi)/stress_dev;
    //     if( bingham_criteria>0.0 ) {
    //         Txx_Novo[i][j] = 2.0*(1 - M.beta)*dudx/bingham_criteria;
    //         Txy_Novo[i][j] = (1 - M.beta)*(dudx + dvdy)/bingham_criteria;
    //         Tyy_Novo[i][j] = 2.0*(1 - M.beta)*dvdy/bingham_criteria;
    //     }
    //     else { // I dont update anything in this case
    //         Txx_Novo[i][j] = Txx_Velho[i][j];
    //         Txy_Novo[i][j] = Txy_Velho[i][j];
    //         Tyy_Novo[i][j] = Tyy_Velho[i][j];
    //     }
    //     return;
    // }


    convecTxx = Conv_Tensor(Txx_Velho, U, V, i, j, "txx", M) - (2.0*Txt_Velho[i][j])*W[i][j]/M.r[i];
    convecTxy = Conv_Tensor(Txy_Velho, U, V, i, j, "txy", M) - (Tyt_Velho[i][j])*W[i][j]/M.r[i];
    convecTyy = Conv_Tensor(Tyy_Velho, U, V, i, j, "tyy", M) - 0.0;
    

    /// === Simulacoes axissimetricas tem umas componentes a mais pra calcular
    if( (M.tipoCoord==AXI_CILINDRICO) ) {
        convecTtt = Conv_Tensor(Ttt_Velho, U, V, i, j, "ttt", M) + (2.0*Txx_Velho[i][j])*W[i][j]/M.r[i];

        /// === Apenas em casos confinados
        if( M.interface.curvaInicial==NULL ) {
            dwdx = DerivadaDWDX(W, i, j, M);
            dwdy = DerivadaDWDY(W, i, j, M);
            convecTxt = Conv_Tensor(Txt_Velho, U, V, i, j, "txt", M) + (Txx_Velho[i][j] - Ttt_Velho[i][j])*W[i][j]/M.r[i];
            convecTyt = Conv_Tensor(Tyt_Velho, U, V, i, j, "tyt", M) + (Txy_Velho[i][j])*W[i][j]/M.r[i];
        }
        
    }

    if( M.tipo_modelo==MODELO_VE || M.tipo_modelo==MODELO_EVP_SARAMITO ) {
        evpt_elasticModulus = 1.0;
        evpt_viscosity = 1.0;
    }
    else {
        evpt_elasticModulus = EVPT_ElasticModulus(M, EVPT_Structure[i][j]);
        evpt_viscosity = EVPT_Viscosity(M, EVPT_Structure[i][j]);

        // PrintDebug("ELASTIC MODULUS: %lf %lf %lf\n", evpt_elasticModulus, M.evpt_m, EVPT_Structure[i][j]);
    }

    xi = (M.tipo_modelo==MODELO_EVPT) ? evpt_elasticModulus/M.We : (1.0 - M.beta)/M.We;
    double evp_term = 1.0;
    if( M.tipo_modelo==MODELO_EVPT )
        evp_term = (evpt_elasticModulus/evpt_viscosity);
    else if( M.tipo_modelo==MODELO_EVP_SARAMITO ) {
        double stress_dev = EVPT_CalculaSigmaDev(M, i, j, U, V, P, Txx_Velho, Txy_Velho, Tyy_Velho, Ttt_Velho);
        double bingham_criteria = (stress_dev - M.Bi)/stress_dev;
        evp_term = ( bingham_criteria>0 ) ? bingham_criteria : 0.0;

        Norm_tau_dev[i][j] = stress_dev;
    }



    /// ===  Equacao Txx (ou Trr)
    parteWi = - (1.0/M.We)*evp_term*( Txx_Velho[i][j] )
             - ( (M.epsilon/(1.0-M.beta))*(Txx_Velho[i][j] + Tyy_Velho[i][j] + Ttt_Velho[i][j])*Txx_Velho[i][j] )
            - ( (M.alpha/(1.0-M.beta))*(Txx_Velho[i][j]*Txx_Velho[i][j] + Txy_Velho[i][j]*Txy_Velho[i][j] + Txt_Velho[i][j]*Txt_Velho[i][j]) )
            ;

    parteCilindrico = ( M.tipoCoord == AXI_CILINDRICO ) ? - 2.0*(W[i][j]/M.r[i])*Txt_Velho[i][j] : 0.0;

    Txx_Novo[i][j] = Txx_Velho[i][j] + M.dt*( - convecTxx + parteCilindrico + 2.0*dudx*Txx_Velho[i][j] + 2.0*dudy*Txy_Velho[i][j] + 2.0*xi*dudx + parteWi );
    
    if( (isnan(Txx_Novo[i][j]) && !isnan(Txx_Velho[i][j])) || (Txx_Novo[i][j]>1e+10) ) {
        // DesenhaMalhaVTK(M, 0);
        // ImprimeInterfaceVTK(M, 1000000);
        DeuProblema("(%d %d): txx nan..... %e %e %lf %lf %lf\n", i, j, evpt_elasticModulus, EVPT_Structure[i][j], parteWi, convecTxx, parteCilindrico);
    }

    /// ===  Equacao Txy (ou Trz)
    parteWi = - (1.0/M.We)*evp_term*( Txy_Velho[i][j] )
            - ( (M.epsilon/(1.0-M.beta))*(Txx_Velho[i][j] + Tyy_Velho[i][j] + Ttt_Velho[i][j])*Txy_Velho[i][j] )
            - ( (M.alpha/(1.0-M.beta))*( Txx_Velho[i][j]*Txy_Velho[i][j] + Tyy_Velho[i][j]*Txy_Velho[i][j] + Txt_Velho[i][j]*Tyt_Velho[i][j] ) )
            ;
    parteCilindrico = ( M.tipoCoord == AXI_CILINDRICO ) ? (dudx + dvdy)*Txy_Velho[i][j] - (W[i][j]/M.r[i])*Tyt_Velho[i][j] : 0.0;
//    parteCilindrico = ( M.tipoCoord == AXI_CILINDRICO ) ? (dudx + dvdy)*Txy_Velho[i][j] : 0.0;
    Txy_Novo[i][j] = Txy_Velho[i][j] + M.dt*( - convecTxy + parteCilindrico + dvdx*Txx_Velho[i][j] + dudy*Tyy_Velho[i][j] + xi*(dudy + dvdx) + parteWi );

    if( isnan(Txy_Novo[i][j]) && !isnan(Txy_Velho[i][j]) ) {
        // DesenhaMalhaVTK(M, 0);
        // ImprimeInterfaceVTK(M, 1000000);
        DeuProblema("(%d %d): txy nan..... %e %e %lf %lf %lf\n", i, j, evpt_elasticModulus, EVPT_Structure[i][j], parteWi, convecTxx, parteCilindrico);
    }

    /// === Equacao Tyy (ou Tzz)
    parteWi = - (1.0/M.We)*evp_term*( Tyy_Velho[i][j] )
            - ( (M.epsilon/(1.0-M.beta))*(Txx_Velho[i][j] + Tyy_Velho[i][j] + Ttt_Velho[i][j])*Tyy_Velho[i][j] )
            - ( (M.alpha/(1.0-M.beta))*( Txy_Velho[i][j]*Txy_Velho[i][j] + Tyt_Velho[i][j]*Tyt_Velho[i][j] + Tyy_Velho[i][j]*Tyy_Velho[i][j] ) )
            ;
    Tyy_Novo[i][j] = Tyy_Velho[i][j] + M.dt*( - convecTyy + 2.0*dvdx*Txy_Velho[i][j] + 2.0*dvdy*Tyy_Velho[i][j] + 2.0*xi*dvdy + parteWi );

    if( isnan(Tyy_Novo[i][j]) && !isnan(Tyy_Velho[i][j]) ) {
        // DesenhaMalhaVTK(M, 0);
        // ImprimeInterfaceVTK(M, 1000000);
        DeuProblema("(%d %d): tyy nan..... %e %e %lf %lf %lf\n", i, j, evpt_elasticModulus, EVPT_Structure[i][j], parteWi, convecTxx, parteCilindrico);
    }

    if( M.tipoCoord == AXI_CILINDRICO ) {

        /// === Equacao T_thetatheta
        parteWi = - (1.0/M.We)*evp_term*( Ttt_Velho[i][j] )
                - ( (M.epsilon/(1.0-M.beta))*(Txx_Velho[i][j] + Tyy_Velho[i][j] + Ttt_Velho[i][j])*Ttt_Velho[i][j] )
                - ( (M.alpha/(1.0-M.beta))*( Txt_Velho[i][j]*Txt_Velho[i][j] + Ttt_Velho[i][j]*Ttt_Velho[i][j] + Tyt_Velho[i][j]*Tyt_Velho[i][j] ) )
                ;
        parteCilindrico = 2.0*(U[i][j]/M.r[i])*Ttt_Velho[i][j];
//        Ttt_Novo[i][j] = Ttt_Velho[i][j] + M.dt*( - convecTtt + parteCilindrico + 2.0*xi*(U[i][j]/M.r[i]) + parteWi );
        Ttt_Novo[i][j] = Ttt_Velho[i][j] + M.dt*( - convecTtt + parteCilindrico + 2.0*dwdx*Txt_Velho[i][j] + 2.0*dwdy*Tyt_Velho[i][j] + 2.0*xi*(U[i][j]/M.r[i]) + parteWi );

        /// === As duas equacoes abaixo soh sao resolvidas em casos confinados. Pra free surface considera T_r_theta = T_z_theta = 0
        if( M.interface.curvaInicial==NULL ) {

            /// ==== Equacao T_r_theta
            parteWi = - (1.0/M.We)*(evpt_elasticModulus/evpt_viscosity)*( Txt_Velho[i][j] )
                    - ( (M.epsilon/(1.0-M.beta))*(Txx_Velho[i][j] + Tyy_Velho[i][j] + Ttt_Velho[i][j])*Txt_Velho[i][j] )
                    - ( (M.alpha/(1.0-M.beta))*( Txx_Velho[i][j]*Txt_Velho[i][j] + Txt_Velho[i][j]*Ttt_Velho[i][j] + Txy_Velho[i][j]*Tyt_Velho[i][j] ) )
                    ;
            parteCilindrico = (U[i][j]/M.r[i])*Txt_Velho[i][j] - (W[i][j]/M.r[i])*Ttt_Velho[i][j];
            Txt_Novo[i][j] = Txt_Velho[i][j] + M.dt*( - convecTxt + parteCilindrico + dudx*Txt_Velho[i][j] + dudy*Tyt_Velho[i][j] + dwdx*Txx_Velho[i][j] + dwdy*Txy_Velho[i][j] + xi*(dwdx - (W[i][j]/M.r[i])) + parteWi );

            /// ==== Equacao T_z_theta
            parteWi = - (1.0/M.We)*(evpt_elasticModulus/evpt_viscosity)*( Tyt_Velho[i][j] )
                    - ( (M.epsilon/(1.0-M.beta))*(Txx_Velho[i][j] + Tyy_Velho[i][j] + Ttt_Velho[i][j])*Tyt_Velho[i][j] )
                    - ( (M.alpha/(1.0-M.beta))*( Txt_Velho[i][j]*Txy_Velho[i][j] + Tyt_Velho[i][j]*Ttt_Velho[i][j] + Tyt_Velho[i][j]*Tyy_Velho[i][j] ) )
                    ;
            parteCilindrico = (U[i][j]/M.r[i])*Tyt_Velho[i][j];
            Tyt_Novo[i][j] = Tyt_Velho[i][j] + M.dt*( - convecTyt + parteCilindrico + dwdx*Txy_Velho[i][j] + dwdy*Tyy_Velho[i][j] + dvdx*Txt_Velho[i][j] + dvdy*Tyt_Velho[i][j] + xi*dwdy + parteWi );
        }

    }

//    if( j==0 || j==M.Ny-1 )
//        PrintDebug("(%2d, %2d): %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf\n", i, j, Txx_Novo[i][j], Txy_Novo[i][j], Tyy_Novo[i][j], Txt_Novo[i][j], Tyt_Novo[i][j], Ttt_Novo[i][j], W[i][j]);






    return;
}

/// FORMULACAO NATURAL-STRESS
void EquacoesConstitutivasNSF(int i, int j, double **LambdaVelho, double **MuVelho, double **NuVelho, double **LambdaNovo,  double **MuNovo, double **NuNovo, double **LambdaAux, double **NuAux,
                              double **UVelho, double **VVelho, double **UNovo, double **VNovo, int n, MALHA M)
{
    double dudx, dudy, dvdx, dudt, dvdt, velV, velU, modV, convecLambda, convecMu, convecNu;
    double parteTempo, parteConvec, parteDivW, parteWi;
    double tol;

    // Ignora se for beta==1, pois eh newtoniano
    if( M.beta==1.0 )
        return;

    tol = M.tolNSF;

    dudx = DerivadaDUDX(UNovo, i, j, M);
    dudy = DerivadaDUDY(UNovo, i, j, (n+1)*M.dt, M);
    dvdx = DerivadaDVDX(VNovo, i, j, M);
    dudt = DerivadaDUDT(UVelho, UNovo, i, j, M);
    dvdt = DerivadaDVDT(VVelho, VNovo, i, j, M);
    velU = 0.5*( UNovo[i+1][j] + UNovo[i][j] );
    velV = 0.5*( VNovo[i][j+1] + VNovo[i][j] );
    modV = velU*velU + velV*velV;


    convecLambda = Conv_Tensor(LambdaAux, UNovo, VNovo, i, j, "lambdaAux", M);
    convecMu = Conv_Tensor(MuVelho, UNovo, VNovo, i, j, "mu", M);
    convecNu = Conv_Tensor(NuAux, UNovo, VNovo, i, j, "nuAux", M);


    parteTempo = (-2.0*MuVelho[i][j]/(modV + tol))*(velV*dudt - velU*dvdt);
    parteConvec = - modV*convecLambda;
    parteDivW = -2.0*MuVelho[i][j]*( (velV*velV - velU*velU)*( dvdx + dudy ) + 4.0*velU*velV*dudx )/( modV + tol );
    parteWi = ( -(LambdaVelho[i][j] - 1.0)/M.We ) //Parte Oldroyd-B
            + ( (M.epsilon/M.We)*(LambdaVelho[i][j] + NuVelho[i][j] - 2.0)*(1.0 - LambdaVelho[i][j]) ) //Parte PTT
            - ( (M.alpha/M.We)*( ((LambdaVelho[i][j] - 1.0)*(LambdaVelho[i][j] - 1.0)) + (MuVelho[i][j]*MuVelho[i][j]) ) ); //Parte Geisekus
    LambdaNovo[i][j] = LambdaVelho[i][j] + M.dt*( parteTempo  + parteConvec + parteDivW + parteWi );


    parteTempo = -( (LambdaVelho[i][j] - NuVelho[i][j])/(modV + tol))*(velU*dvdt - velV*dudt);
    parteConvec = - convecMu;
    parteDivW = -NuVelho[i][j]*( (velV*velV - velU*velU)*( dvdx + dudy ) + 4.0*velU*velV*dudx )/( modV + tol );
    parteWi = - (MuVelho[i][j]/M.We)
            + ( (M.epsilon/M.We)*(LambdaVelho[i][j] + NuVelho[i][j] - 2.0)*(-MuVelho[i][j]) ) //Parte PTT
            + ( (M.alpha/M.We)*(LambdaVelho[i][j] + NuVelho[i][j] - 2.0)*(-MuVelho[i][j]) ); //Parte Geisekus
    MuNovo[i][j] = MuVelho[i][j] + M.dt*( parteTempo  + parteConvec + parteDivW + parteWi );

    parteTempo = (-2.0*MuVelho[i][j]/(modV + tol))*(velU*dvdt - velV*dudt);
    parteConvec = - convecNu/( modV + tol );
    parteDivW = 0.0;
    parteWi = (-(NuVelho[i][j] - 1.0)/M.We)
            + ( (M.epsilon/M.We)*(LambdaVelho[i][j] + NuVelho[i][j] - 2.0)*(1.0-NuVelho[i][j]) ) //Parte PTT
            - ( (M.alpha/M.We)*( ((NuVelho[i][j] - 1.0)*(NuVelho[i][j] - 1.0)) + (MuVelho[i][j]*MuVelho[i][j]) ) ); //Parte Geisekus
    NuNovo[i][j] = NuVelho[i][j] + M.dt*( parteTempo  + parteConvec + parteDivW + parteWi );

    return;
}

void EquacoesConstitutivasLogConformation(int i, int j,
                                          double **PsiVelho_xx, double **PsiVelho_xy, double **PsiVelho_yy, double **PsiVelho_xt, double **PsiVelho_yt, double **PsiVelho_tt,
                                          double **PsiNovo_xx, double **PsiNovo_xy, double **PsiNovo_yy, double **PsiNovo_xt, double **PsiNovo_yt, double **PsiNovo_tt,
                                          double **O_11, double **O_12, double **O_13, double **O_21, double **O_22, double **O_23, double **O_31, double **O_32, double **O_33,
                                          double **Lambda1, double **Lambda2, double **Lambda3,
                                          double **U, double **V, double **W, int n, MALHA M)
{
    double o11, o12, o13 = 0.0, o21, o22, o23 = 0.0, o31 = 0.0, o32 = 0.0, o33 = 0.0;
    double M11, M12, M13, M21, M22, M23, M31, M32, M33;
    double B11, B12, B13, B22, B23, B33;
    double S12 = 0.0, S13 = 0.0, S23 = 0.0, S_Tio12, S_Tio13, S_Tio23; // Matriz sigma
    double dudx, dudy, dvdx, dvdy, dwdx = 0.0, dwdy = 0.0, velW = 0.0, velU = 0.0, coordR = 1.0;
    double convecPsiXX, convecPsiXY, convecPsiYY, convecPsiXT, convecPsiYT, convecPsiTT;

    dudx = DerivadaDUDX(U, i, j, M);
    dudy = DerivadaDUDY(U, i, j, (n+1)*M.dt, M);
    dvdx = DerivadaDVDX(V, i, j, M);
    dvdy = DerivadaDVDY(V, i, j, M);
    convecPsiXX = Conv_Tensor(PsiVelho_xx, U, V, i, j, "PsiXX", M) - (2.0*PsiVelho_xt[i][j])*W[i][j]/M.r[i];
    convecPsiXY = Conv_Tensor(PsiVelho_xy, U, V, i, j, "PsiXY", M) - (PsiVelho_yt[i][j])*W[i][j]/M.r[i];
    convecPsiYY = Conv_Tensor(PsiVelho_yy, U, V, i, j, "PsiYY", M) - 0.0;
    if( M.tipoCoord==AXI_CILINDRICO ) {
        dwdx = DerivadaDWDX(W, i, j, M);
        dwdy = DerivadaDWDY(W, i, j, M);
        velU = 0.5*(U[i][j] + U[i+1][j]);
        velW = W[i][j];
        coordR = M.r[i];

        convecPsiXT = Conv_Tensor(PsiVelho_xt, U, V, i, j, "PsiXT", M) + (PsiVelho_xx[i][j] - PsiVelho_tt[i][j])*W[i][j]/M.r[i];
        convecPsiYT = Conv_Tensor(PsiVelho_yt, U, V, i, j, "PsiYT", M) + (PsiVelho_xy[i][j])*W[i][j]/M.r[i];
        convecPsiTT = Conv_Tensor(PsiVelho_tt, U, V, i, j, "PsiTT", M) + (2.0*PsiVelho_xx[i][j])*W[i][j]/M.r[i];
    }


    o11 = O_11[i][j];
    o12 = O_12[i][j];
    o21 = O_21[i][j];
    o22 = O_22[i][j];
    if( M.tipoCoord==AXI_CILINDRICO ) {
        o13 = O_13[i][j];
        o23 = O_23[i][j];
        o31 = O_31[i][j];
        o32 = O_32[i][j];
        o33 = O_33[i][j];
    }

    if( M.tipoCoord==CARTESIANO ) {
        M11 = o11*o11*dudx + o11*o21*dvdx + o11*o21*dudy + o21*o21*dvdy;
        M12 = o11*o12*dudx + o12*o21*dvdx + o11*o22*dudy + o21*o22*dvdy;
        M21 = o11*o12*dudx + o11*o22*dvdx + o12*o21*dudy + o21*o22*dvdy;
        M22 = o12*o12*dudx + o12*o22*dvdx + o12*o22*dudy + o22*o22*dvdy;
    }
    else {
        M11 = o11*(o11*dudx + o21*dvdx + o31*dwdx) + o21*(o11*dudy + o21*dvdy + o31*dwdy) + o31*(-o11*velW/coordR + o31*velU/coordR);
        M12 = o12*(o11*dudx + o21*dvdx + o31*dwdx) + o22*(o11*dudy + o21*dvdy + o31*dwdy) + o32*(-o11*velW/coordR + o31*velU/coordR);
        M13 = o13*(o11*dudx + o21*dvdx + o31*dwdx) + o23*(o11*dudy + o21*dvdy + o31*dwdy) + o33*(-o11*velW/coordR + o31*velU/coordR);
        M21 = o11*(o12*dudx + o22*dvdx + o32*dwdx) + o21*(o12*dudy + o22*dvdy + o32*dwdy) + o31*(-o12*velW/coordR + o32*velU/coordR);
        M22 = o12*(o12*dudx + o22*dvdx + o32*dwdx) + o22*(o12*dudy + o22*dvdy + o32*dwdy) + o32*(-o12*velW/coordR + o32*velU/coordR);
        M23 = o13*(o12*dudx + o22*dvdx + o32*dwdx) + o23*(o12*dudy + o22*dvdy + o32*dwdy) + o33*(-o12*velW/coordR + o32*velU/coordR);
        M31 = o11*(o13*dudx + o23*dvdx + o33*dwdx) + o21*(o13*dudy + o23*dvdy + o33*dwdy) + o31*(-o13*velW/coordR + o33*velU/coordR);
        M32 = o12*(o13*dudx + o23*dvdx + o33*dwdx) + o22*(o13*dudy + o23*dvdy + o33*dwdy) + o32*(-o13*velW/coordR + o33*velU/coordR);
        M33 = o13*(o13*dudx + o23*dvdx + o33*dwdx) + o23*(o13*dudy + o23*dvdy + o33*dwdy) + o33*(-o13*velW/coordR + o33*velU/coordR);

//        M11 = o11*(o11*dudx + o21*dudy - o31*velW/coordR) + o21*(o11*dvdx + o21*dvdy) + o31*(o11*dwdx + o21*dwdy + o31*velU/coordR);
//        M12 = o12*(o11*dudx + o21*dudy - o31*velW/coordR) + o22*(o11*dvdx + o21*dvdy) + o33*(o11*dwdx + o21*dwdy + o31*velU/coordR);
//        M13 = o13*(o11*dudx + o21*dudy - o31*velW/coordR) + o23*(o11*dvdx + o21*dvdy) + o33*(o11*dwdx + o21*dwdy + o31*velU/coordR);
//        M21 = o11*(o12*dudx + o22*dudy - o32*velW/coordR) + o21*(o12*dvdx + o22*dvdy) + o31*(o12*dwdx + o22*dwdy + o32*velU/coordR);
//        M22 = o12*(o12*dudx + o22*dudy - o32*velW/coordR) + o22*(o12*dvdx + o22*dvdy) + o32*(o12*dwdx + o22*dwdy + o32*velU/coordR);
//        M23 = o13*(o12*dudx + o22*dudy - o32*velW/coordR) + o23*(o12*dvdx + o22*dvdy) + o33*(o12*dwdx + o22*dwdy + o32*velU/coordR);
//        M31 = o11*(o13*dudx + o23*dudy - o33*velW/coordR) + o21*(o13*dvdx + o23*dvdy) + o31*(o13*dwdx + o23*dwdy + o33*velU/coordR);
//        M32 = o12*(o13*dudx + o23*dudy - o33*velW/coordR) + o22*(o13*dvdx + o23*dvdy) + o32*(o13*dwdx + o23*dwdy + o33*velU/coordR);
//        M33 = o13*(o13*dudx + o23*dudy - o33*velW/coordR) + o23*(o13*dvdx + o23*dvdy) + o33*(o13*dwdx + o23*dwdy + o33*velU/coordR);
    }

    B11 = o11*o11*M11 + o12*o12*M22;
    B12 = o11*o21*M11 + o22*o12*M22;
    B22 = o21*o21*M11 + o22*o22*M22;
    if( M.tipoCoord==AXI_CILINDRICO ) {
        B11 += o13*o13*M33;
        B12 += o13*o23*M33;
        B22 += o23*o23*M33;
        B13 = o11*o31*M11 + o12*o32*M22 + o13*o33*M33;
        B23 = o21*o31*M11 + o22*o32*M22 + o23*o33*M33;
        B33 = o31*o31*M11 + o32*o32*M22 + o33*o33*M33;
    }

    if( M.tipoCoord==CARTESIANO ) {
        S_Tio12 = (M12*Lambda2[i][j] + M21*Lambda1[i][j])/(Lambda2[i][j] - Lambda1[i][j] + 1e-20);
        S12 = -o12*o21*S_Tio12 + o11*o22*S_Tio12;
    }
    else {
        S_Tio12 = (M12*Lambda2[i][j] + M21*Lambda1[i][j])/(Lambda2[i][j] - Lambda1[i][j] + 1e-20);
        S_Tio13 = (M13*Lambda3[i][j] + M31*Lambda1[i][j])/(Lambda3[i][j] - Lambda1[i][j] + 1e-20);
        S_Tio23 = (M23*Lambda3[i][j] + M32*Lambda2[i][j])/(Lambda3[i][j] - Lambda2[i][j] + 1e-20);

        S12 = o21*(- o12*S_Tio12 - o13*S_Tio13) + o22*(o11*S_Tio12 - o13*S_Tio23) + o23*(o11*S_Tio13 + o12*S_Tio23);
        S13 = o31*(- o12*S_Tio12 - o13*S_Tio13) + o32*(o11*S_Tio12 - o13*S_Tio23) + o33*(o11*S_Tio13 + o12*S_Tio23);
        S23 = o31*(- o22*S_Tio12 - o23*S_Tio13) + o32*(o21*S_Tio12 - o23*S_Tio23) + o33*(o21*S_Tio13 + o22*S_Tio23);
    }


    static double **parteO_Lambda_O = NULL;
    if( parteO_Lambda_O==NULL ) {
        double valorInicial = 0.0;
        parteO_Lambda_O = (double **)AlocaMatriz(4, 4, sizeof(double), &valorInicial);
    }

    LOG_Termo_Wi(o11, o12, o13, o21, o22, o23, o31, o32, o33,
                 Lambda1[i][j], Lambda2[i][j], Lambda3[i][j],
                 M.alpha, (M.tipoCoord==AXI_CILINDRICO) ? 3 : 2, parteO_Lambda_O);

    double parteSigma, parteWi;

    parteSigma = 2.0*( S12*PsiVelho_xy[i][j] + S13*PsiVelho_xt[i][j] + B11 );
    parteWi = (1.0/M.We)*parteO_Lambda_O[1][1];
    PsiNovo_xx[i][j] = PsiVelho_xx[i][j] + M.dt*( - convecPsiXX + parteSigma + parteWi );

    parteSigma = S12*PsiVelho_yy[i][j] + S13*PsiVelho_yt[i][j] - S12*PsiVelho_xx[i][j] + S23*PsiVelho_xt[i][j] + 2.0*B12;
    parteWi = (1.0/M.We)*parteO_Lambda_O[1][2];
    PsiNovo_xy[i][j] = PsiVelho_xy[i][j] + M.dt*( - convecPsiXY + parteSigma + parteWi );


    parteSigma = 2.0*( - S12*PsiVelho_xy[i][j] + S23*PsiVelho_yt[i][j] + B22 );
    parteWi = (1.0/M.We)*parteO_Lambda_O[2][2];
    PsiNovo_yy[i][j] = PsiVelho_yy[i][j] + M.dt*( - convecPsiYY + parteSigma + parteWi );

    if( M.tipoCoord==AXI_CILINDRICO ) {
        parteSigma = S12*PsiVelho_yt[i][j] + S13*PsiVelho_tt[i][j] - S13*PsiVelho_xx[i][j] - S23*PsiVelho_xy[i][j] + 2.0*B13;
        parteWi = (1.0/M.We)*parteO_Lambda_O[1][3];
        PsiNovo_xt[i][j] = PsiVelho_xt[i][j] + M.dt*( - convecPsiXT + parteSigma + parteWi );

        parteSigma = - S12*PsiVelho_xt[i][j] + S23*PsiVelho_tt[i][j] - S13*PsiVelho_xy[i][j] - S23*PsiVelho_yy[i][j] + 2.0*B23;
        parteWi = (1.0/M.We)*parteO_Lambda_O[2][3];
        PsiNovo_yt[i][j] = PsiVelho_yt[i][j] + M.dt*( - convecPsiYT + parteSigma + parteWi );

        parteSigma = 2.0*( - S13*PsiVelho_xt[i][j] - S23*PsiVelho_yt[i][j] + B33 );
        parteWi = (1.0/M.We)*parteO_Lambda_O[3][3];
        PsiNovo_tt[i][j] = PsiVelho_tt[i][j] + M.dt*( - convecPsiTT + parteSigma + parteWi );
    }


    return;
}

/// === Esta funcao calcula o termo abaixo da equacao constitutiva LOG-giesekus
/// === Termo: O*LambdaInv*O' - I + alfa*( 2I - O*LambdaInv*O' - O*Lambda*O' )
void LOG_Termo_Wi(double o11, double o12, double o13, double o21, double o22, double o23, double o31, double o32, double o33,
                  double lambda1, double lambda2, double lambda3, double alfa, int dimensao, double **Result)
{
    int i, j, k;
    double coef;
    static double O[4][4], Lambda[4];

    O[1][1] = o11;
    O[1][2] = o12;
    O[1][3] = o13;
    O[2][1] = o21;
    O[2][2] = o22;
    O[2][3] = o23;
    O[3][1] = o31;
    O[3][2] = o32;
    O[3][3] = o33;
    Lambda[1] = lambda1;
    Lambda[2] = lambda2;
    Lambda[3] = lambda3;


    for( i=1; i<=dimensao; i++ ) {
        for( j=i; j<=dimensao; j++ ) {

            Result[i][j] = (i==j) ? -1.0 + 2.0*alfa : 0.0;
            for(k=1; k<=dimensao; k++) {
                coef = O[i][k]*O[j][k];
                Result[i][j] += coef/Lambda[k] - alfa*coef*(Lambda[k] + 1.0/Lambda[k]);

//                double temp1 = coef/Lambda[k] - alfa*coef*(Lambda[k] + 1.0/Lambda[k]);
//                if( temp1!=0 )
//                    PrintDebug("RRRRRESULT: %lf %lf\n", coef, temp1);
//                if( Result[i][j]!=0 )
//                    PrintDebug("RRRRRESULT: %lf\n", Result[i][j]);
            }


        }
    }

}

double EVPT_CalculaSigmaDev(MALHA M, int i, int j, 
                            double **U, double **V, double **P, 
                            double **Txx, double **Txy, double **Tyy, double **Ttt)
{
    
    

    // Sigma eh o total stress tensor. Soh preciso das componentes xx e yy, pois vou fazer o trace
    static double sigma[3][3];


    // sigma[0][0] = - M.Re*P[i][j] + 2.0*M.eta_inf*dudx + Txx[i][j];
    // sigma[1][1] = - M.Re*P[i][j] + 2.0*M.eta_inf*dvdy + Tyy[i][j];
    // sigma[0][1] = sigma[1][0] = M.eta_inf*(dudy + dvdx) + Txy[i][j];
    // double sigma_trace = sigma[0][0] + sigma[1][1];
    
    // sigma[0][0] = sigma[0][0] - sigma_trace/3.0;
    // sigma[1][1] = sigma[1][1] - sigma_trace/3.0;
    
    // // Elevando a matriz sigma ao quadrado
    // sigma[0][0] = sigma[0][0]*sigma[0][0] + sigma[0][1]*sigma[0][1];
    // sigma[1][1] = sigma[1][1]*sigma[1][1] + sigma[0][1]*sigma[0][1];

    // sigma_trace = sigma[0][0] + sigma[1][1];
    // return sqrt( 0.5*sigma_trace );



    /// ==== Cassio version - EVPT - Freeflow
    // double dudx = DerivadaDUDX(U, i, j, M);
    // double dudy = DerivadaDUDY(U, i, j, (M.passoTemporal+1)*M.dt, M);
    // double dvdx = DerivadaDVDX(V, i, j, M);
    // double dvdy = DerivadaDVDY(V, i, j, M);
    // sigma[0][0] = 2.0*M.eta_inf*dudx + Txx[i][j];
    // sigma[1][1] = 2.0*M.eta_inf*dvdy + Tyy[i][j];
    // sigma[0][1] = sigma[1][0] = M.eta_inf*(dudy + dvdx) + Txy[i][j];
    // return sqrt(0.5*(sigma[0][0] * sigma[0][0] + 2 * sigma[0][1] * sigma[0][1] + sigma[1][1] * sigma[1][1]));


    /// ==== Saramito version - EVP model
    double trace = Txx[i][j] + Tyy[i][j] + Ttt[i][j];
    double one_third = 1.0/3.0;
    sigma[0][0] = Txx[i][j] - one_third*trace;
    sigma[1][1] = Tyy[i][j] - one_third*trace;
    sigma[2][2] = Ttt[i][j] - one_third*trace;
    sigma[0][1] = sigma[1][0] = Txy[i][j];
    sigma[0][2] = sigma[2][0] = 0.0;
    
    return sqrt( sigma[0][0]*sigma[0][0] + 2.0*sigma[0][1]*sigma[0][1] + sigma[1][1]*sigma[1][1] + sigma[2][2]*sigma[2][2] );
}

double DerivadaDUDX(double **U, int i, int j, MALHA M)
{
    return (U[i+1][j] - U[i][j])/M.dx[i];
}

double DerivadaDUDY(double **U, int i, int j, double t, MALHA M)
{
    static double uCimaDir = 0.0, uCimaEsq = 0.0, uBaixoDir = 0.0, uBaixoEsq = 0.0, uCentroEsq = 0.0, uCentroDir = 0.0;
    static double h1, h2;


    uCentroEsq = U[i][j];
    uCentroDir = U[i+1][j];

    ///Encontrando o uCimaEsq
    if( CIMA_U(i, j) ) {
        if( M.pontosV[i][j+1].tipoBoundary == NOSLIP )
            uCimaEsq = - U[i][j] + 2.0*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[i][j+1].tipoBoundary == SLIP )
            uCimaEsq = - U[i][j] + 2.0*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[i][j+1].tipoBoundary == INFLOW )
            uCimaEsq = - U[i][j];
        else if( M.pontosV[i][j+1].tipoBoundary == NEUMANN )
            uCimaEsq = U[i][j];
        else if( M.pontosV[i][j+1].tipoBoundary == SIMETRIA )
            uCimaEsq = U[i][j];
        else DeuProblema("Problema DerivadaDUDY Cima esquerda %d %d\n", i, j);
    }
    else
        uCimaEsq = U[i][j+1];

    ///Encontrando uCimaDir
    if( CIMA_U(i+1, j) ) {
        if( M.pontosV[i][j+1].tipoBoundary == NOSLIP )
            uCimaDir = - U[i+1][j] + 2.0*FuncaoValorNoSlipHorizontal( M, i+1, j, M.x[i+1], t, U );
        else  if( M.pontosV[i][j+1].tipoBoundary == SLIP )
            uCimaDir = - U[i+1][j] + 2.0*FuncaoValorNoSlipHorizontal( M, i+1, j, M.x[i+1], t, U );
        else if( M.pontosV[i][j+1].tipoBoundary == INFLOW )
            uCimaDir = - U[i+1][j];
        else if( M.pontosV[i][j+1].tipoBoundary == NEUMANN )
            uCimaDir = U[i+1][j];
        else if( M.pontosV[i][j+1].tipoBoundary == SIMETRIA )
            uCimaDir = U[i+1][j];
        else DeuProblema("Problema DerivadaDUDY Cima direita %d %d\n", i, j);
    }
    else
        uCimaDir = U[i+1][j+1];

    ///Encontrando uBaixoEsq
    if( BAIXO_U(i, j) ) {
        if( M.pontosV[i][j].tipoBoundary == NOSLIP )
            uBaixoEsq = - U[i][j] + 2.0*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else  if( M.pontosV[i][j].tipoBoundary == SLIP )
            uBaixoEsq = - U[i][j] + 2.0*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[i][j].tipoBoundary == INFLOW )
            uBaixoEsq = - U[i][j];
        else if( M.pontosV[i][j].tipoBoundary == NEUMANN )
            uBaixoEsq = U[i][j];
        else if( M.pontosV[i][j].tipoBoundary == SIMETRIA )
            uBaixoEsq = U[i][j];
        else DeuProblema("Problema DerivadaDUDY baixo esquerda %d %d\n", i, j);
    }
    else {
        if( j==0 )
            DeuProblema("AAAA %d %d\n", i, j);
        uBaixoEsq = U[i][j-1];
    }

    if( BAIXO_U(i+1, j) ) {
        if( M.pontosV[i][j].tipoBoundary == NOSLIP )
            uBaixoDir = - U[i+1][j] + 2.0*FuncaoValorNoSlipHorizontal( M, i+1, j, M.x[i+1], t, U );
        else if( M.pontosV[i][j].tipoBoundary == SLIP )
            uBaixoDir = - U[i+1][j] + 2.0*FuncaoValorNoSlipHorizontal( M, i+1, j, M.x[i+1], t, U );
        else if( M.pontosV[i][j].tipoBoundary == INFLOW )
            uBaixoDir = - U[i+1][j];
        else if( M.pontosV[i][j].tipoBoundary == NEUMANN )
            uBaixoDir = U[i+1][j];
        else if( M.pontosV[i][j].tipoBoundary == SIMETRIA )
            uBaixoDir = U[i+1][j];
        else DeuProblema("Problema DerivadaDUDY baixo direita %d %d\n", i, j);
    }
    else
        uBaixoDir = U[i+1][j-1];


    h1 = (j==0) ? M.dy[j] : 0.5*( M.dy[j-1] + M.dy[j] );
    h2 = (j==M.Ny-1) ? M.dy[j] : 0.5*( M.dy[j] + M.dy[j+1] );

    return ( (-h2/(h1*(h1+h2)))*0.5*( uBaixoEsq + uBaixoDir ) ) +
                 ( ((h2-h1)/(h1*h2))*0.5*( uCentroEsq + uCentroDir ) ) +
                 ( (h1/(h2*(h1+h2)))*0.5*( uCimaEsq + uCimaDir ) );
}

double DerivadaDVDX(double **V, int i, int j, MALHA M)
{
    static double vDirCima = 0.0, vEsqCima = 0.0, vDirBaixo = 0.0, vEsqBaixo = 0.0, vCentroCima = 0.0, vCentroBaixo = 0.0;
    static double h1 = 0.0, h2 = 0.0;

    vCentroBaixo = V[i][j];
    vCentroCima = V[i][j+1];


    if( ESQUERDA_V(i, j) ) {
        if( M.pontosU[i][j].tipoBoundary == NEUMANN )
            vEsqBaixo = V[i][j];
        else if( M.pontosU[i][j].tipoBoundary == NOSLIP )
            vEsqBaixo = - V[i][j] + 2.0*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i][j].tipoBoundary == SLIP )
            vEsqBaixo = - V[i][j] + 2.0*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i][j].tipoBoundary == INFLOW )
            vEsqBaixo = - V[i][j];
        else if( M.pontosU[i][j].tipoBoundary == EIXO_AXISSIMETRIA )
            vEsqBaixo = V[i][j];
        else if( M.pontosU[i][j].tipoBoundary == KNOWN_PRESSURE )
            vEsqBaixo = 0.0;
        else DeuProblema("Problema Esquerda DerivadaDVDX: %d %d\n", i, j);
    }
    else
        vEsqBaixo = V[i-1][j];

    if( ESQUERDA_V(i, j+1) ) {
        if( M.pontosU[i][j].tipoBoundary == NEUMANN )
            vEsqCima = V[i][j+1];
        else if( M.pontosU[i][j].tipoBoundary == NOSLIP )
            vEsqCima = - V[i][j+1] + 2.0*FuncaoValorNoSlipVertical( M, i, j+1, M.y[j], 0.0, V );
        else if( M.pontosU[i][j].tipoBoundary == SLIP )
            vEsqCima = - V[i][j+1] + 2.0*FuncaoValorNoSlipVertical( M, i, j+1, M.y[j], 0.0, V );
        else if( M.pontosU[i][j].tipoBoundary == INFLOW )
            vEsqCima = - V[i][j+1];
        else if( M.pontosU[i][j].tipoBoundary == EIXO_AXISSIMETRIA )
            vEsqCima = V[i][j+1];
        else if( M.pontosU[i][j].tipoBoundary == KNOWN_PRESSURE )
            vEsqCima = 0.0;
        else DeuProblema("Problema Esquerda DerivadaDVDX: %d %d\n", i, j);
    }
    else
        vEsqCima = V[i-1][j+1];

    ///Encontrando vDirBaixo
    if( DIREITA_V(i, j) ) {
        if( M.pontosU[i+1][j].tipoBoundary==NEUMANN )
            vDirBaixo = V[i][j];
        else if( M.pontosU[i+1][j].tipoBoundary==NOSLIP )
            vDirBaixo = - V[i][j] + 2.0*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][j].tipoBoundary==SLIP )
            vDirBaixo = - V[i][j] + 2.0*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][j].tipoBoundary==INFLOW )
            vDirBaixo = - V[i][j];
        else if( M.pontosU[i+1][j].tipoBoundary==KNOWN_PRESSURE )
            vDirBaixo = 0.0;
        else DeuProblema("Problema Direita DerivadaDVDX: %d %d\n", i, j);
    }
    else
        vDirBaixo = V[i+1][j];

    if( DIREITA_V(i, j+1) ) {
        if( M.pontosU[i+1][j].tipoBoundary==NEUMANN )
            vDirCima = V[i][j+1];
        else if( M.pontosU[i+1][j].tipoBoundary==NOSLIP )
            vDirCima = - V[i][j+1] + 2.0*FuncaoValorNoSlipVertical( M, i, j+1, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][j].tipoBoundary==SLIP )
            vDirCima = - V[i][j+1] + 2.0*FuncaoValorNoSlipVertical( M, i, j+1, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][j].tipoBoundary==INFLOW )
            vDirCima = - V[i][j+1];
        else if( M.pontosU[i+1][j].tipoBoundary==KNOWN_PRESSURE )
            vDirCima = 0.0;
        else DeuProblema("Problema Direita DerivadaDVDX: %d %d\n", i, j);
    }
    else
        vDirCima = V[i+1][j+1];

    h1 = (i==0) ? M.dx[i] : 0.5*( M.dx[i-1] + M.dx[i] );
    h2 = (i==M.Nx-1) ? M.dx[i] : 0.5*( M.dx[i] + M.dx[i+1] );

    return ( (-h2/(h1*(h1+h2)))*0.5*( vEsqBaixo + vEsqCima ) ) +
                 ( ((h2-h1)/(h1*h2))*0.5*( vCentroBaixo + vCentroCima ) ) +
                 ( (h1/(h2*(h1+h2)))*0.5*( vDirBaixo + vDirCima ) );
}

double DerivadaDVDY(double **V, int i, int j, MALHA M)
{
    return (V[i][j+1] - V[i][j])/M.dy[j];
}

double DerivadaDWDX(double **W, int i, int j, MALHA M)
{
    double wEsq = 0.0, wDir = 0.0, wCentro = 0.0;
    double valorDirichlet;

    wCentro = W[i][j];

    //Encontrando o W(i-1, j)
    if( i!=0 && M.pontosW[i-1][j].tipo==FULL ) {
        wEsq = W[i-1][j];
    }
    else {
        if( M.pontosU[i][j].movimento==MOVIMENTO_THETA )
            valorDirichlet = M.pontosU[i][j].valorDirichlet;
        else
            valorDirichlet = 0.0;

        if( M.pontosU[i][j].tipoBoundary == NOSLIP )
            wEsq = - W[i][j] + 2.0*valorDirichlet;
        else if( M.pontosU[i][j].tipoBoundary == INFLOW )
            wEsq = - W[i][j];
        else if( M.pontosU[i][j].tipoBoundary == EIXO_AXISSIMETRIA )
            wEsq = W[i][j];
        else {
            DeuProblema("problema Esquerda - DWDY... %d %d\n", i, j);
        }
    }

    //Encontarndo o W(i+1, j)
    if( i!=M.Nx-1 && M.pontosW[i+1][j].tipo==FULL ) {
        wDir = W[i+1][j];
    }
    else {

        if( M.pontosU[i+1][j].movimento==MOVIMENTO_THETA )
            valorDirichlet = M.pontosU[i+1][j].valorDirichlet;
        else
            valorDirichlet = 0.0;

        if( M.pontosU[i+1][j].tipoBoundary == NOSLIP )
            wDir = - W[i][j] + 2.0*valorDirichlet;
        else if( M.pontosU[i+1][j].tipoBoundary == INFLOW )
            wDir = - W[i][j];
        else {
            DeuProblema("problema Direita - DWDX... %d %d\n", i, j);
        }
    }

    double h1 = (i==0) ? M.dx[i] : 0.5*( M.dx[i-1] + M.dx[i] );
    double h2 = (i==M.Nx-1) ? M.dx[i] : 0.5*( M.dx[i] + M.dx[i+1] );

    return ( (-h2/(h1*(h1+h2)))*wEsq ) +
                 ( ((h2-h1)/(h1*h2))*wCentro ) +
                 ( (h1/(h2*(h1+h2)))*wDir );
}

double DerivadaDWDY(double **W, int i, int j, MALHA M)
{
    double wBaixo = 0.0, wCima = 0.0, wCentro = 0.0;

    wCentro = W[i][j];

    //Encontrando o W(i, j-1)
    if( j!=0 && M.pontosW[i][j-1].tipo==FULL ) {
        wBaixo = W[i][j-1];
    }
    else {
        if( M.pontosV[i][j].tipoBoundary == NOSLIP )
            wBaixo = - W[i][j];
        else if( M.pontosV[i][j].tipoBoundary == INFLOW )
            wBaixo = - W[i][j];
        else if( M.pontosV[i][j].tipoBoundary == NEUMANN )
            wBaixo = W[i][j];
        else {
            DeuProblema("problema Baixo - DWDY... %d %d\n", i, j);
        }
    }

    //Encontrando o W(i, j+1)
    if( j!=M.Ny-1 && M.pontosW[i][j+1].tipo==FULL ) {
        wCima = W[i][j+1];
    }
    else {
        if( M.pontosV[i][j+1].tipoBoundary == NOSLIP )
            wCima = - W[i][j];
        else if( M.pontosV[i][j+1].tipoBoundary == INFLOW )
            wCima = - W[i][j];
        else if( M.pontosV[i][j+1].tipoBoundary == NEUMANN )
            wCima = W[i][j];
        else {
            DeuProblema("problema Cima - DWDY\n", i, j);
        }
    }

    double h1 = (j==0) ? M.dy[j] : 0.5*( M.dy[j-1] + M.dy[j] );
    double h2 = (j==M.Ny-1) ? M.dy[j] : 0.5*( M.dy[j] + M.dy[j+1] );

    return ( (-h2/(h1*(h1+h2)))*wBaixo ) +
                 ( ((h2-h1)/(h1*h2))*wCentro ) +
                 ( (h1/(h2*(h1+h2)))*wCima );
}

double DerivadaDUDT(double **UVelho, double **UNovo, int i, int j, MALHA M)
{
    return 0.5*(UNovo[i+1][j] + UNovo[i][j] - UVelho[i+1][j] - UVelho[i][j])/M.dt;
}

double DerivadaDVDT(double **VVelho, double **VNovo, int i, int j, MALHA M)
{
    return 0.5*(VNovo[i][j+1] + VNovo[i][j] - VVelho[i][j+1] - VVelho[i][j])/M.dt;
}

double EVPT_ElasticModulus(MALHA M, double Lambda)
{
    return exp( M.evpt_m*( 1.0/Lambda - 1.0 ) );
}

double EVPT_Viscosity(MALHA M, double Lambda)
{
    return (pow(M.eta_0/M.eta_inf, Lambda) - 1.0)*M.eta_inf;
}

int ProximoDaParede_Hibrido(int i, int j, MALHA M)
{
    int grauHibrido;
    //int grauDiagonal;

    grauHibrido = M.grauHibrido;
    //grauDiagonal = 1;

    ///BLOCOS CSF - teste cassio
//    int alturaBloco = 50;
//    if( i<70 ) {
//        if( (j<alturaBloco) || ((259-j)<alturaBloco) )
//            return 1;
//    }



//    int iQuina1, iQuina2, jQuina1, jQuina2;
//
//    int raioQuina = 40;
//
//    iQuina1 = 69;
//    jQuina1 = 90;
//    iQuina2 = 69;
//    jQuina2 = 169;
//
//    if( (abs(i-iQuina1)<=raioQuina) && (abs(j-jQuina1)<=raioQuina) )
//        return 0;
//    if( (abs(i-iQuina2)<=raioQuina) && (abs(j-jQuina2)<=raioQuina) )
//        return 0;

//
//    if( ((260-j)<=grauHibrido) && (i>=50) ) {
//        return 1;
//    }
//    if( (j<grauHibrido) && (i>=50) ) {
//        return 1;
//    }
//
//    if( (j>=210) && ((70-i)<=grauHibrido) ) {
//        return 1;
//    }
//    if( (j<50) && ((70-i)<=grauHibrido) ) {
//        return 1;
//    }

    //Casos obvios, pois estao nos extremos do dominio geral
    if( (i<grauHibrido) || ((M.Nx-i)<=grauHibrido) || (j<grauHibrido) || ((M.Ny-j)<=grauHibrido) )
        return 1;

    //Casos nao-obvios, pois estao proximos a uma parede no meio do dominio (tipo o problema da contracao)
    for( ; grauHibrido>0; grauHibrido-- ) {
        if( M.pontosP[i-grauHibrido][j].tipo==EMPTY )
            return 1;
        if( M.pontosP[i+grauHibrido][j].tipo==EMPTY )
            return 1;
        if( M.pontosP[i][j-grauHibrido].tipo==EMPTY )
            return 1;
        if( M.pontosP[i][j+grauHibrido].tipo==EMPTY )
            return 1;

//        //Aplicando tambem na quina
//        if( M.pontosP[i-grauDiagonal][j-grauDiagonal].tipo==EMPTY )
//            return 1;
//        if( M.pontosP[i-grauDiagonal][j+grauDiagonal].tipo==EMPTY )
//            return 1;
//        if( M.pontosP[i+grauDiagonal][j-grauDiagonal].tipo==EMPTY )
//            return 1;
//        if( M.pontosP[i+grauDiagonal][j+grauDiagonal].tipo==EMPTY )
//            return 1;
    }

    return 0;
}


double CalculaDetMinA(MALHA M, double **Txx, double **Txy, double **Tyy)
{
    int i, j;
    double menorDet = 1e+10, det;
    double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
    static double A[3][3];

    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            A[1][1] = 1.0 + Txx[i][j]/xi;
            A[1][2] = Txy[i][j]/xi;
            A[2][1] = Txy[i][j]/xi;
            A[2][2] = 1.0 + Tyy[i][j]/xi;

            det = A[1][1]*A[2][2] - A[1][2]*A[2][1];

            if( menorDet>det )
                menorDet = det;
        }
    }

    return menorDet;
}

void CalculaResiduosEqsConstitutivas(MALHA *M, double **UVelho, double **UNovo , double **VVelho, double **VNovo,
                          double **TxxVelho, double **TxxNovo, double **TxyVelho, double **TxyNovo, double **TyyVelho, double **TyyNovo,
                          double **LambdaVelho, double **LambdaNovo, double **MuVelho, double **MuNovo, double **NuVelho, double **NuNovo,
                          double **LambdaAux, double **NuAux,
                          double **PsiVelho_xx, double **PsiNovo_xx, double **PsiVelho_xy, double **PsiNovo_xy, double **PsiVelho_yy, double **PsiNovo_yy,
                          double **ResTxx, double **ResTxy, double **ResTyy,
                          double **ResLambda, double **ResMu, double **ResNu,
                          double **ResPsiXX, double **ResPsiXY, double **ResPsiYY)
{
    int i, j;
    double dudx, dudy, dvdx, dvdy, dudt, dvdt;
    double convecTxx, convecTxy, convecTyy, convecLambda, convecMu, convecNu, convecPsiXX, convecPsiXY, convecPsiYY;
    double parteTempo, parteConv, parteWi, parte2, parteDivW, parteDelT;
    double velU, velV, modV;

    //double tol = M->tolNSF;
    double tol = 0.0;

    if( M->formulacaoVisco==VISCO_NSF || M->formulacaoVisco==VISCO_HIBRIDO ) {
        int linha;
        for( linha=0; linha<M->qtdIncognitasP; linha++ ) {
            i = M->indicesP[linha].i;
            j = M->indicesP[linha].j;

            velU = 0.5*( UNovo[i+1][j] + UNovo[i][j] );
            velV = 0.5*( VNovo[i][j+1] + VNovo[i][j] );
            modV = velU*velU + velV*velV;

            LambdaAux[i][j] = LambdaNovo[i][j]/( modV + tol );
            NuAux[i][j] = modV*NuNovo[i][j];
//
//            LambdaAux[i][j] = LambdaVelho[i][j]/( modV + tol);
//            NuAux[i][j] = modV*NuVelho[i][j];
        }
    }

    for( i=0; i<M->Nx; i++) {
        for( j=0; j<M->Ny; j++ ) {

            //Vou calcular apenas em pontos FULL
            if( M->pontosP[i][j].tipo!=FULL )
                continue;

            dudx = DerivadaDUDX(UNovo, i, j, *M);
            dudy = DerivadaDUDY(UNovo, i, j, 0.0, *M);
            dvdx = DerivadaDVDX(VNovo, i, j, *M);
            dvdy = DerivadaDVDY(VNovo, i, j, *M);
            dudt = DerivadaDUDT(UVelho, UNovo, i, j, *M);
            dvdt = DerivadaDVDT(VVelho, VNovo, i, j, *M);
            convecTxx = Conv_Tensor(TxxNovo, UNovo, VNovo, i, j, "txx", *M);
            convecTxy = Conv_Tensor(TxyNovo, UNovo, VNovo, i, j, "txy", *M);
            convecTyy = Conv_Tensor(TyyNovo, UNovo, VNovo, i, j, "tyy", *M);

            parteTempo = (TxxNovo[i][j] - TxxVelho[i][j])/(M->dt);
            parteConv = convecTxx;
            parte2 = - 2.0*( dudx*TxxNovo[i][j] + dudy*TxyNovo[i][j] );
            parteWi = - (2.0/M->We)*(1.0 - M->beta)*dudx;
            ResTxx[i][j] = ((1.0/M->We)*TxxNovo[i][j]) + parteTempo + parteConv + parte2 + parteWi;

            parteTempo = (TxyNovo[i][j] - TxyVelho[i][j])/(M->dt);
            parteConv = convecTxy;
            parte2 = - ( dvdx*TxxNovo[i][j] + dudy*TyyNovo[i][j] );
            parteWi = - (1.0/M->We)*(1.0 - M->beta)*(dudy + dvdx);
            ResTxy[i][j] = ((1.0/M->We)*TxyNovo[i][j]) + parteTempo + parteConv + parte2 + parteWi;

            parteTempo = (TyyNovo[i][j] - TyyVelho[i][j])/(M->dt);
            parteConv = convecTyy;
            parte2 = - 2.0*( dvdx*TxyNovo[i][j] + dvdy*TyyNovo[i][j] );
            parteWi = - (2.0/M->We)*(1.0 - M->beta)*dvdy;
            ResTyy[i][j] = ((1.0/M->We)*TyyNovo[i][j]) + parteTempo + parteConv + parte2 + parteWi;




            velU = 0.5*( UNovo[i+1][j] + UNovo[i][j] );
            velV = 0.5*( VNovo[i][j+1] + VNovo[i][j] );
            modV = velU*velU + velV*velV;
            convecLambda = Conv_Tensor(LambdaAux, UNovo, VNovo, i, j, "lambdaAux", *M);
            convecMu = Conv_Tensor(MuNovo, UNovo, VNovo, i, j, "mu", *M);
            convecNu = Conv_Tensor(NuAux, UNovo, VNovo, i, j, "nuAux", *M);

            parteDelT = (LambdaNovo[i][j] - LambdaVelho[i][j])/M->dt;
            parteTempo = (2.0*MuNovo[i][j]/(modV + tol))*(velV*dudt - velU*dvdt);
            parteConv = modV*convecLambda;
            parteDivW = 2.0*MuNovo[i][j]*( (velV*velV - velU*velU)*( dvdx + dudy ) + 4.0*velU*velV*dudx )/( modV + tol );
            parteWi = (LambdaNovo[i][j] - 1.0)/M->We;
            ResLambda[i][j] = parteDelT + parteTempo + parteConv + parteDivW + parteWi;

            parteDelT = (MuNovo[i][j] - MuVelho[i][j])/M->dt;
            parteTempo = ((LambdaNovo[i][j] - NuNovo[i][j])/(modV + tol ))*(velU*dvdt - velV*dudt);
            parteConv = convecMu;
            parteDivW = NuNovo[i][j]*( (velV*velV - velU*velU)*( dvdx + dudy ) + 4.0*velU*velV*dudx )/( modV + tol );
            parteWi = (MuNovo[i][j]/M->We); //Parte Geisekus
            ResMu[i][j] = parteDelT + parteTempo + parteConv + parteDivW + parteWi;

            parteDelT = (NuNovo[i][j] - NuVelho[i][j])/M->dt;
            parteTempo = (2.0*MuNovo[i][j]/(modV + tol))*(velU*dvdt - velV*dudt);
            parteConv = convecNu/( modV + tol );
            parteDivW = 0.0;
            parteWi = ((NuNovo[i][j] - 1.0)/M->We);
            ResNu[i][j] = parteDelT + parteTempo + parteConv + parteDivW + parteWi;






            convecPsiXX = Conv_Tensor(PsiNovo_xx, UNovo, VNovo, i, j, "PsiXX", *M);
            convecPsiXY = Conv_Tensor(PsiNovo_xy, UNovo, VNovo, i, j, "PsiXY", *M);
            convecPsiYY = Conv_Tensor(PsiNovo_yy, UNovo, VNovo, i, j, "PsiYY", *M);

            double xi = (1.0 - M->beta)/(/*M.Re**/M->We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            A[1][1] = 1.0 + TxxNovo[i][j]/xi;
            A[1][2] = TxyNovo[i][j]/xi;
            A[2][1] = TxyNovo[i][j]/xi;
            A[2][2] = 1.0 + TyyNovo[i][j]/xi;

            CalculaAutovalores_Jacobi(A, 2, auto_valores, aux1, aux2, O, &nrot);

            double o11 = O[1][1];
            double o12 = O[1][2];
            double o21 = O[2][1];
            double o22 = O[2][2];

            double M11 = o11*o11*dudx + o11*o21*dvdx + o11*o21*dudy + o21*o21*dvdy;
            double M12 = o11*o12*dudx + o12*o21*dvdx + o11*o22*dudy + o21*o22*dvdy;
            double M21 = o11*o12*dudx + o11*o22*dvdx + o12*o21*dudy + o21*o22*dvdy;
            double M22 = o12*o12*dudx + o12*o22*dvdx + o12*o22*dudy + o22*o22*dvdy;

            double B11 = o11*o11*M11 + o12*o12*M22;
            double B12 = o11*o21*M11 + o22*o12*M22;
            double B22 = o21*o21*M11 + o22*o22*M22;

            double S_Tio = (M12*auto_valores[2] + M21*auto_valores[1])/(auto_valores[2] - auto_valores[1] + 1e-20);
            double S12 = -o12*o21*S_Tio + o11*o22*S_Tio;

            double parteSigma, parteWi, parteTempo;

            parteTempo = (PsiNovo_xx[i][j] - PsiVelho_xx[i][j])/M->dt;
            parteSigma = 2.0*( PsiNovo_xy[i][j]*S12 + B11 );
            parteWi = (1.0/M->We)*( (o11*o11/auto_valores[1]) + (o12*o12/auto_valores[2]) - 1.0 );
            ResPsiXX[i][j] = parteTempo + convecPsiXX - parteSigma - parteWi;

            parteTempo = (PsiNovo_xy[i][j] - PsiVelho_xy[i][j])/M->dt;
            parteSigma = S12*( PsiNovo_yy[i][j] - PsiNovo_xx[i][j] ) + 2.0*B12;
            parteWi = (1.0/M->We)*( (o11*o21/auto_valores[1]) + (o12*o22/auto_valores[2]) );
            ResPsiXY[i][j] = parteTempo + convecPsiXY - parteSigma - parteWi;

            parteTempo = (PsiNovo_yy[i][j] - PsiVelho_yy[i][j])/M->dt;
            parteSigma = 2.0*( - PsiNovo_xy[i][j]*S12 + B22 );
            parteWi = (1.0/M->We)*( (o21*o21/auto_valores[1]) + (o22*o22/auto_valores[2]) - 1.0 );
            ResPsiYY[i][j] = parteTempo + convecPsiYY - parteSigma - parteWi;

            ///Apenas testando pra ser se da zero usando valores em n
//            velU = 0.5*( UNovo[i+1][j] + UNovo[i][j] );
//            velV = 0.5*( VNovo[i][j+1] + VNovo[i][j] );
//            modV = velU*velU + velV*velV;
//            convecLambda = Conv_Tensor(LambdaAux, UNovo, VNovo, i, j, "lambdaAux", *M);
//            convecMu = Conv_Tensor(MuVelho, UNovo, VNovo, i, j, "mu", *M);
//            convecNu = Conv_Tensor(NuAux, UNovo, VNovo, i, j, "nuAux", *M);
//
//            parteDelT = (LambdaNovo[i][j] - LambdaVelho[i][j])/M->dt;
//            parteTempo = (2.0*MuVelho[i][j]/(modV + tol))*(velV*dudt - velU*dvdt);
//            parteConv = modV*convecLambda;
//            parteDivW = 2.0*MuVelho[i][j]*( (velV*velV - velU*velU)*( dvdx + dudy ) + 4.0*velU*velV*dudx )/( modV + tol );
//            parteWi = (LambdaVelho[i][j] - 1.0)/M->We;
//            ResLambda[i][j] = parteDelT + parteTempo + parteConv + parteDivW + parteWi;
//
//            parteDelT = (MuNovo[i][j] - MuVelho[i][j])/M->dt;
//            parteTempo = ((LambdaVelho[i][j] - NuVelho[i][j])/(modV + tol))*(velU*dvdt - velV*dudt);
//            parteConv = convecMu;
//            parteDivW = NuVelho[i][j]*( (velV*velV - velU*velU)*( dvdx + dudy ) + 4.0*velU*velV*dudx )/( modV + tol );
//            parteWi = (MuVelho[i][j]/M->We); //Parte Geisekus
//            ResMu[i][j] = parteDelT + parteTempo + parteConv + parteDivW + parteWi;
//
//            parteDelT = (NuNovo[i][j] - NuVelho[i][j])/M->dt;
//            parteTempo = (2.0*MuVelho[i][j]/(modV + tol ))*(velU*dvdt - velV*dudt);
//            parteConv = convecNu/( modV + tol );
//            parteDivW = 0.0;
//            parteWi = ((NuVelho[i][j] - 1.0)/M->We);
//            ResNu[i][j] = parteDelT + parteTempo + parteConv + parteDivW + parteWi;
        }
    }
}

void CalculaTermoDominanteEqsConstituvas(int i, int j,
                                         double **Txx_Velho, double **Txy_Velho, double **Tyy_Velho, double **Txt_Velho, double **Tyt_Velho, double **Ttt_Velho,
                                         double **Txx_Novo, double **Txy_Novo, double **Tyy_Novo, double **Txt_Novo, double **Tyt_Novo, double **Ttt_Novo,
                                         double **U, double **V, double **W,
                                         double **TermoDominanteXX, double **TermoDominanteXY, double **TermoDominanteYY,
                                         int n, MALHA M)
{
    double dudx, dudy, dvdx, dvdy;//, dwdx, dwdy;
    double convecTxx, convecTxy, convecTyy;//, convecTxt, convecTyt, convecTtt;
    double xi;

    // Ignora se for beta==1
    if( M.beta==1.0 )
        return;

    xi = (1.0 - M.beta)/(/*M.Re**/M.We);

    dudx = DerivadaDUDX(U, i, j, M);
    dudy = DerivadaDUDY(U, i, j, (n+1)*M.dt, M);
    dvdx = DerivadaDVDX(V, i, j, M);
    dvdy = DerivadaDVDY(V, i, j, M);
//    dwdx = (M.tipoCoord==AXI_CILINDRICO) ? DerivadaDWDX(W, i, j, M) : 0.0;
//    dwdy = (M.tipoCoord==AXI_CILINDRICO) ? DerivadaDWDY(W, i, j, M) : 0.0;
    convecTxx = Conv_Tensor(Txx_Velho, U, V, i, j, "txx", M);
    convecTxy = Conv_Tensor(Txy_Velho, U, V, i, j, "txy", M);
    convecTyy = Conv_Tensor(Tyy_Velho, U, V, i, j, "tyy", M);
//    convecTxt = (M.tipoCoord==AXI_CILINDRICO) ? Conv_Tensor(Txt_Velho, U, V, i, j, "txt", M) : 0.0;
//    convecTyt = (M.tipoCoord==AXI_CILINDRICO) ? Conv_Tensor(Tyt_Velho, U, V, i, j, "tyt", M) : 0.0;
//    convecTtt = (M.tipoCoord==AXI_CILINDRICO) ? Conv_Tensor(Ttt_Velho, U, V, i, j, "ttt", M) : 0.0;

    /// FABIANO: Verificando qual o termo dominante
    /// === Verificando o termo dominante pra imprimir no VTK
    double termo1, termo2, termo3;
    // Primeiramente componente xx
    termo1 = fabs(Txx_Novo[i][j])/M.We;
    termo2 = fabs( ((Txx_Novo[i][j] - Txx_Velho[i][j])/M.dt) + convecTxx - 2.0*dudx*Txx_Novo[i][j] - 2.0*dudy*Txy_Novo[i][j] );
    termo3 = fabs( 2.0*xi*dudx );
    if( (termo1>termo2) && (termo1>termo3) )
        TermoDominanteXX[i][j] = 1.0;
    else if( termo2>termo3 )
        TermoDominanteXX[i][j] = 2.0;
    else
        TermoDominanteXX[i][j] = 3.0;

    // Agora componente xy
    termo1 = fabs(Txy_Novo[i][j])/M.We;
    termo2 = fabs( ((Txy_Novo[i][j] - Txy_Velho[i][j])/M.dt) + convecTxy - dvdx*Txx_Novo[i][j] - dudy*Tyy_Novo[i][j] );
    termo3 = fabs( xi*(dudy + dvdx) );
    if( (termo1>termo2) && (termo1>termo3) )
        TermoDominanteXY[i][j] = 1.0;
    else if( termo2>termo3 )
        TermoDominanteXY[i][j] = 2.0;
    else
        TermoDominanteXY[i][j] = 3.0;

    // Agora componente YY
    termo1 = fabs(Tyy_Novo[i][j])/M.We;
    termo2 = fabs( ((Tyy_Novo[i][j] - Tyy_Velho[i][j])/M.dt) + convecTyy - 2.0*dvdx*Txy_Novo[i][j] - 2.0*dvdy*Tyy_Novo[i][j] );
    termo3 = fabs( 2.0*xi*dvdy );
    if( (termo1>termo2) && (termo1>termo3) )
        TermoDominanteXY[i][j] = 1.0;
    else if( termo2>termo3 )
        TermoDominanteXY[i][j] = 2.0;
    else
        TermoDominanteXY[i][j] = 3.0;

    return;
}


