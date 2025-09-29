#include "Visualizacao.h"


void ImprimeArquivosUVPGnuplot(double **U, double **V, double **W, double **P, MALHA Malha, int n)
{
    FILE *arqU, *arqV, *arqW, *arqP;
    int i, j;
    char nomeArq[100];

    //Criando os arquivos
    sprintf(nomeArq, "ArquivosPlot/U/VelocU%d.txt", n);
    arqU = fopen(nomeArq, "wt");
    sprintf(nomeArq, "ArquivosPlot/V/VelocV%d.txt", n);
    arqV = fopen(nomeArq, "wt");
    sprintf(nomeArq, "ArquivosPlot/W/VelocW%d.txt", n);
    arqW = fopen(nomeArq, "wt");
    sprintf(nomeArq, "ArquivosPlot/P/PressaoP%d.txt", n);
    arqP = fopen(nomeArq, "wt");

    //Plotando U
    for(j=0; j<Malha.Ny; j++) {
        for(i=0; i<=Malha.Nx; i++) {
            /*if( i==0 || j==0 || i==Malha.Nx || j==Malha.Ny-1 )
                continue;
            if( Malha.pontosU[i][j].tipo==SURFACE || Malha.pontosU[i+1][j].tipo==SURFACE || Malha.pontosU[i-1][j].tipo==SURFACE
               || Malha.pontosU[i][j-1].tipo==SURFACE || Malha.pontosU[i][j+1].tipo==SURFACE )
                fprintf(arqU, "%lf %lf %e\n", Malha.y[j] + 0.5*Malha.dy[j], Malha.x[i], 0.0);
            else*/
//            if( Malha.pontosU[i][j].tipo==FULL )
//                fprintf(arqU, "%lf %lf %e\n", Malha.y[j] + 0.5*Malha.dy[j], Malha.x[i], 0.0);
//            else
                fprintf(arqU, "%lf %lf %e\n", Malha.y[j] + 0.5*Malha.dy[j], Malha.x[i], U[i][j]);
                //fprintf(arqU, "%lf %lf %d\n", Malha.y[j] + 0.5*Malha.dy[j], Malha.x[i], Malha.pontosU[i][j].tipo);
        }
        fprintf(arqU, "\n");
    }

    //Plotando V
    for(j=0; j<=Malha.Ny; j++) {
        for(i=0; i<Malha.Nx; i++) {
//            if( Malha.pontosV[i][j].tipo==FULL )
//                fprintf(arqV, "%lf %lf %lf\n", Malha.y[j], Malha.x[i] + 0.5*Malha.dx[i], 0.0);
//            else
                fprintf(arqV, "%lf %lf %lf\n", Malha.y[j], Malha.x[i] + 0.5*Malha.dx[i], V[i][j]);
                //fprintf(arqV, "%lf %lf %d\n", Malha.y[j], Malha.x[i] + 0.5*Malha.dx[i], Malha.pontosV[i][j].tipo);
        }
        fprintf(arqV, "\n");
    }

    if( Malha.tipoCoord==AXI_CILINDRICO ) {
        for(j=0; j<Malha.Ny; j++) {
            for(i=0; i<Malha.Nx; i++)
                fprintf(arqW, "%lf %lf %lf\n", Malha.y[j] + 0.5*Malha.dy[j], Malha.x[i] + 0.5*Malha.dx[i], W[i][j]);
            fprintf(arqW, "\n");
        }
    }

    //Plotando P
    for(j=0; j<Malha.Ny; j++) {
        for(i=0; i<Malha.Nx; i++)
            fprintf(arqP, "%lf %lf %lf\n", Malha.y[j] + 0.5*Malha.dy[j], Malha.x[i] + 0.5*Malha.dx[i], P[i][j]);
            //fprintf(arqP, "%lf %lf %d\n", Malha.y[j] + 0.5*Malha.dy[j], Malha.x[i] + 0.5*Malha.dx[i], Malha.pontosP[i][j].tipo);
        fprintf(arqP, "\n");
    }


    fclose(arqU);
    fclose(arqV);
    fclose(arqW);
    fclose(arqP);
    return;
}

void ImprimeArquivoVTK(double **U, double **V, double **W, double **P,
                       double **Txx, double **Txy, double **Tyy, double **Txt, double **Tyt, double **Ttt,
                       double **Lambda, double **Mu, double **Nu,
                       double **EVPT_Structure, double **Norm_tau_dev,
                       double **MuNewtGen,
                       double **TermoDominanteXX, double **TermoDominanteXY, double **TermoDominanteYY,
                       double **ResidualFaceX, double **ResidualFaceY, double **ResidualCell,
                       MALHA Malha, int n)
{
    FILE *arq;
    int i, j;
    char nomeArq[500];
    double uBaixo, uCima, vEsq, vDir;


    //Criando os arquivos
    sprintf(nomeArq, "ArquivosPlot/VTK/%s/Dados-N%d.vtk", Malha.nomeArquivo, n);
    arq = fopen(nomeArq, "wt");
    if( !arq )
        DeuProblema("\n\n ImprimeArquivoVTK: Problema no fopen do arquivo. \n\n");

    //Cabecalho
    fprintf(arq, "# vtk DataFile Version 2.0\nStructured Grid Example\nASCII\n");
    fprintf(arq, "DATASET RECTILINEAR_GRID\n");
    fprintf(arq, "DIMENSIONS %d %d 1\n", Malha.Nx+1, Malha.Ny+1);
    fprintf(arq, "X_COORDINATES %d float\n", Malha.Nx+1);
    for( i=0; i<=Malha.Nx; i++ )
        fprintf(arq, "%lf ", Malha.x[i]);
    fprintf(arq, "\nY_COORDINATES %d float\n", Malha.Ny+1);
    for( i=0; i<=Malha.Ny; i++ )
        fprintf(arq, "%lf ", Malha.y[i]);
    fprintf(arq, "\nZ_COORDINATES %d float\n", 1);
    fprintf(arq, "0.0\n\n");

    //Vetores velocidade (estou imprimindo valores de U e V nos vertices das celulas fazendo interpolacao)
    fprintf(arq, "POINT_DATA %d\n", (Malha.Nx+1)*(Malha.Ny+1));
    fprintf(arq, "VECTORS velocidade float\n");
    for( j=0; j<=Malha.Ny; j++ ){
        for( i=0; i<=Malha.Nx; i++ ) {

            vEsq = (i==0) ? V[i][j] : V[i-1][j];
            vDir = (i==Malha.Nx) ? V[i-1][j] : V[i][j];
            uBaixo = (j==0) ? U[i][j] : U[i][j-1];
            uCima = (j==Malha.Ny) ? U[i][j-1] : U[i][j];

            double u, v;
            u = fabs( (0.5*(uBaixo + uCima)) ) < 1e-10 ? 0.0 : 0.5*(uBaixo + uCima);
            v = fabs( (0.5*(vEsq + vDir)) ) < 1e-10 ? 0.0 :  0.5*(vEsq + vDir);
            //fprintf(arq, "%e %e 0.0000\n", 0.5*(uBaixo + uCima), 0.5*(vEsq + vDir) );
            fprintf(arq, "%e %e 0.0000\n", u, v );
        }
    }
    fprintf(arq, "\n");


    /* Velocity scalars */
    fprintf(arq, "SCALARS velocidade-u float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Malha.Ny; j++) {
        for (i=0; i<=Malha.Nx; i++) {
            uBaixo = (j==0) ? U[i][j] : U[i][j-1];
            uCima = (j==Malha.Ny) ? U[i][j-1] : U[i][j];
            double u;
            u = fabs( (0.5*(uBaixo + uCima)) ) < 1e-10 ? 0.0 : 0.5*(uBaixo + uCima);
            //fprintf(arq, "%e\n", 0.5*(uBaixo + uCima));
            fprintf(arq, "%e\n",u);
        }
    }
    fprintf(arq, "\n");

    fprintf(arq, "SCALARS velocidade-v float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Malha.Ny; j++) {
        for(i=0; i<=Malha.Nx; i++) {
            vEsq = (i==0) ? V[i][j] : V[i-1][j];
            vDir = (i==Malha.Nx) ? V[i-1][j] : V[i][j];

            double v;
            v = fabs( (0.5*(vEsq + vDir)) ) < 1e-10 ? 0.0 :  0.5*(vEsq + vDir);
            fprintf(arq, "%e\n", v);
            //fprintf(arq, "%e\n", 0.5*(vEsq + vDir));
        }
    }

    // fprintf(arq, "SCALARS residual-u float 1\n");
    // fprintf(arq, "LOOKUP_TABLE default\n");
    // for(j=0; j<=Malha.Ny; j++) {
    //     for (i=0; i<=Malha.Nx; i++) {
    //         uBaixo = (j==0) ? ResidualFaceX[i][j] : ResidualFaceX[i][j-1];
    //         uCima = (j==Malha.Ny) ? ResidualFaceX[i][j-1] : ResidualFaceX[i][j];
    //         double u;
    //         u = fabs( (0.5*(uBaixo + uCima)) ) < 1e-10 ? 0.0 : 0.5*(uBaixo + uCima);
    //         //fprintf(arq, "%e\n", 0.5*(uBaixo + uCima));
    //         fprintf(arq, "%e\n",u);
    //     }
    // }
    // fprintf(arq, "\n");

    // fprintf(arq, "SCALARS residual-v float 1\n");
    // fprintf(arq, "LOOKUP_TABLE default\n");
    // for(j=0; j<=Malha.Ny; j++) {
    //     for(i=0; i<=Malha.Nx; i++) {
    //         vEsq = (i==0) ? ResidualFaceY[i][j] : ResidualFaceY[i-1][j];
    //         vDir = (i==Malha.Nx) ? ResidualFaceY[i-1][j] : ResidualFaceY[i][j];

    //         double v;
    //         v = fabs( (0.5*(vEsq + vDir)) ) < 1e-10 ? 0.0 :  0.5*(vEsq + vDir);
    //         fprintf(arq, "%e\n", v);
    //         //fprintf(arq, "%e\n", 0.5*(vEsq + vDir));
    //     }
    // }

    fprintf(arq, "CELL_DATA %d\n", (Malha.Nx)*(Malha.Ny));
    
    
    fprintf(arq, "SCALARS pressao float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<Malha.Ny; j++)
        for(i=0; i<Malha.Nx; i++)
            fprintf(arq, "%e\n", P[i][j]);
    fprintf(arq, "\n");

    fprintf(arq, "SCALARS curvatura float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<Malha.Ny; j++)
        for(i=0; i<Malha.Nx; i++)
            fprintf(arq, "%e\n", Malha.curvaturaCelulas[i][j]);
    fprintf(arq, "\n");

    // fprintf(arq, "SCALARS ResidualVisc float 1\n");
    // fprintf(arq, "LOOKUP_TABLE default\n");
    // for(j=0; j<Malha.Ny; j++)
    //     for(i=0; i<Malha.Nx; i++)
    //         fprintf(arq, "%e\n", ResidualCell[i][j]);
    // fprintf(arq, "\n");

     if( Malha.tipoCoord==AXI_CILINDRICO ) {
         fprintf(arq, "SCALARS velocidade-w float 1\n");
         fprintf(arq, "LOOKUP_TABLE default\n");
         for(j=0; j<Malha.Ny; j++)
             for(i=0; i<Malha.Nx; i++)
                 fprintf(arq, "%e\n", (Malha.tipoCoord==AXI_CILINDRICO) ? W[i][j] : 0.0);
         fprintf(arq, "\n");
     }

     fprintf(arq, "SCALARS Txx float 1\n");
     fprintf(arq, "LOOKUP_TABLE default\n");
     for(j=0; j<Malha.Ny; j++)
         for(i=0; i<Malha.Nx; i++)
             fprintf(arq, "%e\n", Txx[i][j]);
     fprintf(arq, "\n");

     fprintf(arq, "SCALARS Txy float 1\n");
     fprintf(arq, "LOOKUP_TABLE default\n");
     for(j=0; j<Malha.Ny; j++)
         for(i=0; i<Malha.Nx; i++)
             fprintf(arq, "%e\n", Txy[i][j]);
     fprintf(arq, "\n");

     fprintf(arq, "SCALARS Tyy float 1\n");
     fprintf(arq, "LOOKUP_TABLE default\n");
     for(j=0; j<Malha.Ny; j++)
         for(i=0; i<Malha.Nx; i++)
             fprintf(arq, "%e\n", Tyy[i][j]);
     fprintf(arq, "\n");

    // fprintf(arq, "SCALARS NormTauDev float 1\n");
    // fprintf(arq, "LOOKUP_TABLE default\n");
    // for(j=0; j<Malha.Ny; j++)
    //     for(i=0; i<Malha.Nx; i++)
    //         fprintf(arq, "%e\n", Norm_tau_dev[i][j]);
    // fprintf(arq, "\n");

    // fprintf(arq, "SCALARS EffectiveVisc float 1\n");
    // fprintf(arq, "LOOKUP_TABLE default\n");
    // for(j=0; j<Malha.Ny; j++)
    //     for(i=0; i<Malha.Nx; i++)
    //         fprintf(arq, "%e\n", MuNewtGen[i][j]);
    // fprintf(arq, "\n");

    // if( Malha.tipoCoord==AXI_CILINDRICO ) {
    //     fprintf(arq, "SCALARS Txt float 1\n");
    //     fprintf(arq, "LOOKUP_TABLE default\n");
    //     for(j=0; j<Malha.Ny; j++)
    //         for(i=0; i<Malha.Nx; i++)
    //             fprintf(arq, "%e\n", (Malha.tipoCoord==AXI_CILINDRICO) ? Txt[i][j] : 0.0);
    //     fprintf(arq, "\n");

    //     fprintf(arq, "SCALARS Tyt float 1\n");
    //     fprintf(arq, "LOOKUP_TABLE default\n");
    //     for(j=0; j<Malha.Ny; j++)
    //         for(i=0; i<Malha.Nx; i++)
    //             fprintf(arq, "%e\n", (Malha.tipoCoord==AXI_CILINDRICO) ? Tyt[i][j] : 0.0);
    //     fprintf(arq, "\n");

    //     fprintf(arq, "SCALARS Ttt float 1\n");
    //     fprintf(arq, "LOOKUP_TABLE default\n");
    //     for(j=0; j<Malha.Ny; j++)
    //         for(i=0; i<Malha.Nx; i++)
    //             fprintf(arq, "%e\n", (Malha.tipoCoord==AXI_CILINDRICO) ? Ttt[i][j] : 0.0);
    //     fprintf(arq, "\n");
    // }

//    fprintf(arq, "SCALARS Lambda float 1\n");
//    fprintf(arq, "LOOKUP_TABLE default\n");
//    for(j=0; j<Malha.Ny; j++)
//        for(i=0; i<Malha.Nx; i++)
//            fprintf(arq, "%e\n", Lambda[i][j]);
//    fprintf(arq, "\n");
//
//    fprintf(arq, "SCALARS Mu float 1\n");
//    fprintf(arq, "LOOKUP_TABLE default\n");
//    for(j=0; j<Malha.Ny; j++)
//        for(i=0; i<Malha.Nx; i++)
//            fprintf(arq, "%e\n", Mu[i][j]);
//    fprintf(arq, "\n");
//
//    fprintf(arq, "SCALARS Nu float 1\n");
//    fprintf(arq, "LOOKUP_TABLE default\n");
//    for(j=0; j<Malha.Ny; j++)
//        for(i=0; i<Malha.Nx; i++)
//            fprintf(arq, "%e\n", Nu[i][j]);
//    fprintf(arq, "\n");
//
//    fprintf(arq, "SCALARS TermoDominanteXX float 1\n");
//    fprintf(arq, "LOOKUP_TABLE default\n");
//    for(j=0; j<Malha.Ny; j++)
//        for(i=0; i<Malha.Nx; i++)
//            fprintf(arq, "%e\n", TermoDominanteXX[i][j]);
//    fprintf(arq, "\n");
//
//    fprintf(arq, "SCALARS TermoDominanteXY float 1\n");
//    fprintf(arq, "LOOKUP_TABLE default\n");
//    for(j=0; j<Malha.Ny; j++)
//        for(i=0; i<Malha.Nx; i++)
//            fprintf(arq, "%e\n", TermoDominanteXY[i][j]);
//    fprintf(arq, "\n");
//
//    fprintf(arq, "SCALARS TermoDominanteYY float 1\n");
//    fprintf(arq, "LOOKUP_TABLE default\n");
//    for(j=0; j<Malha.Ny; j++)
//        for(i=0; i<Malha.Nx; i++)
//            fprintf(arq, "%e\n", TermoDominanteYY[i][j]);
//    fprintf(arq, "\n");


    // if( Malha.tipo_modelo!=MODELO_VE ) {
    //     fprintf(arq, "SCALARS StructureEVPT float 1\n");
    //     fprintf(arq, "LOOKUP_TABLE default\n");
    //     for(j=0; j<Malha.Ny; j++)
    //         for(i=0; i<Malha.Nx; i++)
    //             fprintf(arq, "%e\n", EVPT_Structure[i][j]);
    //     fprintf(arq, "\n");
    // }

    fclose(arq);
    //VecRestoreArray(xGlobal, &sol);
    return;
}

void ImprimeInterfaceVTK(MALHA Malha, int n)
{
    FILE *arq;
    PONTO *ponto;
    CURVA *curva;
    int i, *qtdPontos, qtdTotalPontos, qtdCurvas, indicePonto;
    char nomeArq[500];// nomeLocal[600], temp[100];

    if( Malha.interface.curvaInicial==NULL )
        return;

     //Criando o nome do arquivo
//    strcpy(nomeLocal, Malha.nomeArquivo);
//    sprintf(temp, "%d", n);
//    while( SubstituiString(nomeLocal, "[Passo]", temp) );
//    sprintf(temp, "%g", Malha.dt);
//    while( SubstituiString(nomeLocal, "[dt]", temp) );
//    sprintf(temp, "%g", Malha.beta);
//    while( SubstituiString(nomeLocal, "[beta]", temp) );
//    sprintf(temp, "%g", Malha.We);
//    while( SubstituiString(nomeLocal, "[We]", temp) );
//    sprintf(temp, "%1.0e", Malha.tolNSF);
//    while( SubstituiString(nomeLocal, "[tolNSF]", temp) );


    //Criando os arquivos
    //Criando os arquivos
    sprintf(nomeArq, "ArquivosPlot/VTK/%s/Interface-N%d.vtk", Malha.nomeArquivo, n);
    arq = fopen(nomeArq, "wt");
    if( !arq )
        DeuProblema("\n\n ImprimeInterfaceVTK: Problema no fopen do arquivo \n\n");


    /// Primeiramente contando quantos pontos e curvas tem
    qtdCurvas = 0;
    for( curva=Malha.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva )
        qtdCurvas++;

    qtdPontos = (int *)malloc( qtdCurvas*sizeof(int) );
    i = 0;
    qtdTotalPontos = 0;
    for( curva=Malha.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        qtdPontos[i] = 0;
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;
            qtdPontos[i]++;
            qtdTotalPontos++;
        }
        i++;
    }




//    Cabecalho
    PetscFPrintf(PETSC_COMM_SELF, arq, "# vtk DataFile Version 2.0\nStructured Grid Example\nASCII\n");
    PetscFPrintf(PETSC_COMM_SELF, arq, "DATASET POLYDATA\n");
    PetscFPrintf(PETSC_COMM_SELF, arq, "POINTS %d float\n", qtdTotalPontos);

    for( curva=Malha.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;
            fprintf(arq, "%e %e %e\n", ponto->x, ponto->y, 0.0);
        }
    }

    i = 0;
    indicePonto = 0;
    PetscFPrintf(PETSC_COMM_SELF, arq, "POLYGONS %d %d\n", qtdCurvas, qtdTotalPontos+qtdCurvas);
    for( curva=Malha.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        PetscFPrintf(PETSC_COMM_SELF, arq, "%d", qtdPontos[i]);
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;
            PetscFPrintf(PETSC_COMM_SELF, arq, " %d", indicePonto++);
        }
        PetscFPrintf(PETSC_COMM_SELF, arq, "\n");
        i++;
    }

    free(qtdPontos);
    fclose(arq);
    return;
}

void ImprimeArquivoBinario(double **U, double **V, MALHA Malha, int n)
{
    FILE *arq;
    int i, j;
    char nomeArq[900];
    static double **U_grid = NULL, **V_grid = NULL;

    if( U_grid==NULL ) {
        double valorInicial = 0.0;
        U_grid = (double **)AlocaMatriz(Malha.Nx+1, Malha.Ny+1, sizeof(double), &valorInicial);
        V_grid = (double **)AlocaMatriz(Malha.Nx+1, Malha.Ny+1, sizeof(double), &valorInicial);
    }

    // === Colocando U e V no grid
    for( i=0; i<=Malha.Nx; i++ ) {
        for( j=0; j<Malha.Ny; j++ ) {

            // Forcando um noslip na parede inferior
            if( j==0 )
                U_grid[i][j] = V_grid[i][j] = 0.0;
            else {
                double vEsq = (i==0) ? V[i][j] : V[i-1][j];
                double vDir = (i==Malha.Nx) ? V[i-1][j] : V[i][j];
                double uBaixo = (j==0) ? U[i][j] : U[i][j-1];
                double uCima = (j==Malha.Ny) ? U[i][j-1] : U[i][j];
                U_grid[i][j] = 0.5*( uBaixo + uCima );
                V_grid[i][j] = 0.5*( vEsq + vDir );
            }
            
        }
    }



    //Criando os arquivos
    sprintf(nomeArq, "ArquivosPlot/VTK/%s/velBinary-N%d.bin", Malha.nomeArquivo, n);
    arq = fopen(nomeArq, "w");
    if( !arq )
        DeuProblema("\n\n ImprimeArquivoVTK: Problema no fopen do arquivo. \n\n");
    
    for( i=0; i<=Malha.Nx; i++ )
        fwrite(U_grid[i], sizeof(double), Malha.Ny+1, arq);
    for( i=0; i<=Malha.Nx; i++ )
        fwrite(V_grid[i], sizeof(double), Malha.Ny+1, arq);
    fclose(arq);
    return;
}

void ImprimeCorte(double **U, double **V, double **W, int Indice, const char *Direcao, const char *Variavel, const char *NomeArquivo, MALHA M)
{
    FILE *arq;
    int i;

    arq = fopen(NomeArquivo, "wt");
//    fprintf(arq, "matriz=[");

    if( !strcmp(Direcao, "vertical") ) {

        if( !strcmp(Variavel, "U") ) {
            for( i=0; i<M.Ny; i++ )
                fprintf(arq, "%lf %lf;\n", 0.5*(M.y[i] + M.y[i+1]), U[Indice][i]);
        }
        else if( !strcmp(Variavel, "V") ) {
            for( i=0; i<=M.Ny; i++ )
                fprintf(arq, "%lf %lf;\n", M.y[i], V[Indice][i]);
        }
        else if( !strcmp(Variavel, "W") ) {
            for( i=0; i<M.Ny; i++ )
                fprintf(arq, "%lf %lf;\n", 0.5*(M.y[i] + M.y[i+1]), W[Indice][i]);
        }

    }
    else if( !strcmp(Direcao, "horizontal") ) {

        if( !strcmp(Variavel, "U") ) {
            for( i=0; i<=M.Nx; i++ )
                fprintf(arq, "%lf %lf;\n", M.x[i], U[i][Indice]);
        }
        else if( !strcmp(Variavel, "V") ) {
            for( i=0; i<M.Nx; i++ )
                fprintf(arq, "%lf %lf;\n", 0.5*(M.x[i] + M.x[i+1]), V[i][Indice]);
        }
        else if( !strcmp(Variavel, "W") ) {
            for( i=0; i<M.Nx; i++ )
                fprintf(arq, "%lf %lf;\n", 0.5*(M.x[i] + M.x[i+1]), W[i][Indice]);
        }

    }

//    fprintf(arq, "];");
    fclose(arq);

    return;
}

void PlotU(MALHA Malha, int n)
{
//    FILE *arq;
//    char nomeArq[200], comando[200];
//
//    sprintf(nomeArq, "ArquivosPlot/U/ScriptU%d.txt", n);
//    arq = fopen(nomeArq, "wt");
//
//    sprintf(comando, "gnuplot ArquivosPlot/U/ScriptU%d.txt", n);
//
//    fprintf(arq, "unset key\n");
//    fprintf(arq, "set tic scale 1\n");
//    fprintf(arq, "set palette defined ( 0 \"blue\", 1 \"cyan\", 2 \"green\", 3 \"yellow\", 4 \"red\" )#\n");
//    fprintf(arq, "set xlabel \"x\"\n");
//    fprintf(arq, "set ylabel \"y\"\n");
//    fprintf(arq, "#set cbrange [-0.5:1.5]\n");
//    fprintf(arq, "#set cblabel \"Escala\"\n");
//    fprintf(arq, "#unset cbtics\n");
//    fprintf(arq, "set xrange [%lf:%lf]\n", Malha.xMin, Malha.xMax);
//    fprintf(arq, "set yrange [-%lf:%lf]\n", Malha.yMin, Malha.yMax);
//    //fprintf(arq, "set size ratio -1\n");
//    fprintf(arq, "set view map\n");
//    fprintf(arq, "set term pdf\n");
//    fprintf(arq, "set output \"ArquivosPlot/U/PlotU%d.pdf\"\n", n);
//    fprintf(arq, "plot 'ArquivosPlot/U/VelocU%d.txt' using 2:1:3 with image\n", n);
//    fprintf(arq, "unset output\n");
//
//    fclose(arq);
//    system(comando);
}

void PlotV(MALHA Malha, int n)
{
//    FILE *arq;
//    char nomeArq[200], comando[200];
//
//    sprintf(nomeArq, "ArquivosPlot/V/ScriptV%d.txt", n);
//    arq = fopen(nomeArq, "wt");
//
//    sprintf(comando, "gnuplot ArquivosPlot/V/ScriptV%d.txt", n);
//
//    fprintf(arq, "unset key\n");
//    fprintf(arq, "set tic scale 1\n");
//    fprintf(arq, "set palette defined ( 0 \"blue\", 1 \"cyan\", 2 \"green\", 3 \"yellow\", 4 \"red\" )#\n");
//    fprintf(arq, "set xlabel \"x\"\n");
//    fprintf(arq, "set ylabel \"y\"\n");
//    fprintf(arq, "#set cbrange [-0.5:1.5]\n");
//    fprintf(arq, "#set cblabel \"Escala\"\n");
//    fprintf(arq, "#unset cbtics\n");
//    fprintf(arq, "set xrange [%lf:%lf]\n", Malha.xMin, Malha.xMax);
//    fprintf(arq, "set yrange [-%lf:%lf]\n", Malha.yMin, Malha.yMax);
//    //fprintf(arq, "set size ratio -1\n");
//    fprintf(arq, "set view map\n");
//    fprintf(arq, "set term pdf\n");
//    fprintf(arq, "set output \"ArquivosPlot/V/PlotV%d.pdf\"\n", n);
//    fprintf(arq, "plot 'ArquivosPlot/V/VelocV%d.txt' using 2:1:3 with image\n", n);
//    fprintf(arq, "unset output\n");
//
//    fclose(arq);
//    system(comando);
}

void PlotP(MALHA Malha, int n)
{
//    FILE *arq;
//    char nomeArq[200], comando[200];
//
//    sprintf(nomeArq, "ArquivosPlot/P/ScriptP%d.txt", n);
//    arq = fopen(nomeArq, "wt");
//
//    sprintf(comando, "gnuplot ArquivosPlot/P/ScriptP%d.txt", n);
//
//    fprintf(arq, "unset key\n");
//    fprintf(arq, "set tic scale 1\n");
//    fprintf(arq, "set palette defined ( 0 \"blue\", 1 \"cyan\", 2 \"green\", 3 \"yellow\", 4 \"red\" )#\n");
//    fprintf(arq, "set xlabel \"x\"\n");
//    fprintf(arq, "set ylabel \"y\"\n");
//    fprintf(arq, "#set cbrange [-0.5:1.5]\n");
//    fprintf(arq, "#set cblabel \"Escala\"\n");
//    fprintf(arq, "#unset cbtics\n");
//    fprintf(arq, "set xrange [%lf:%lf]\n", Malha.xMin, Malha.xMax);
//    fprintf(arq, "set yrange [-%lf:%lf]\n", Malha.yMin, Malha.yMax);
//    //fprintf(arq, "set size ratio -1\n");
//    fprintf(arq, "set view map\n");
//    fprintf(arq, "set term pdf\n");
//    fprintf(arq, "set output \"ArquivosPlot/P/P%d.pdf\"\n", n);
//    fprintf(arq, "plot 'ArquivosPlot/P/PressaoP%d.txt' using 2:1:3 with image\n", n);
//    fprintf(arq, "unset output\n");
//
//    fclose(arq);
//    system(comando);
}

void PlotW(MALHA Malha, int n)
{
//    FILE *arq;
//    char nomeArq[200], comando[200];
//
//    if( Malha.tipoCoord!=AXI_CILINDRICO )
//        return;
//
//    sprintf(nomeArq, "ArquivosPlot/W/ScriptW%d.txt", n);
//    arq = fopen(nomeArq, "wt");
//
//    sprintf(comando, "gnuplot ArquivosPlot/W/ScriptW%d.txt", n);
//
//    fprintf(arq, "unset key\n");
//    fprintf(arq, "set tic scale 1\n");
//    fprintf(arq, "set palette defined ( 0 \"blue\", 1 \"cyan\", 2 \"green\", 3 \"yellow\", 4 \"red\" )#\n");
//    fprintf(arq, "set xlabel \"x\"\n");
//    fprintf(arq, "set ylabel \"y\"\n");
//    fprintf(arq, "#set cbrange [-0.5:1.5]\n");
//    fprintf(arq, "#set cblabel \"Escala\"\n");
//    fprintf(arq, "#unset cbtics\n");
//    fprintf(arq, "set xrange [%lf:%lf]\n", Malha.xMin, Malha.xMax);
//    fprintf(arq, "set yrange [-%lf:%lf]\n", Malha.yMin, Malha.yMax);
//    fprintf(arq, "set view map\n");
//    fprintf(arq, "set term pdf\n");
//    fprintf(arq, "set output \"ArquivosPlot/W/W%d.pdf\"\n", n);
//    fprintf(arq, "plot 'ArquivosPlot/W/VelocW%d.txt' using 2:1:3 with image\n", n);
//    fprintf(arq, "unset output\n");
//
//    fclose(arq);
//    system(comando);
}



///Superficie livre
void ImprimePontosCurva(CURVA Curva, FILE *Arquivo)
{
    PONTO *ponto;

    //Percorre cada ponto desta curva
    LOOP_CURVA(&Curva) {
        ponto = loop_stuff.ponto;
        fprintf(Arquivo, "[%G %G]\n", ponto->x, ponto->y);
    }

    return;
}

void ImprimeCurvasEPontos(INTERFACE Interface)
{
    CURVA *curva;
    PONTO *ponto;
    int i = 0;

    for( curva=Interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        printf("\n NOVA CURVA %d\n", i++);
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;
            printf("[%.15lf %.15lf]\n", ponto->x, ponto->y);
        }
    }


    return;
}

void PlotInterfaceGnuplot(MALHA Malha, char DesenhaMalha, char DesenhaNormais, int n)
{
//    CURVA *curva;
//    PONTO *ponto;
//    FILE *arq;
//    char NomeArquivo[300];
//    double x, y;
//    int i = 0, j;
//
//    //Verifica se o pipe gnuplot ja foi criado, se nao foi cria ele
//    //if( pipeGnuplot==NULL )
//    //    InicializaPipeGnuplot();
//
//    //Percorre cada curva da interface
//    for( curva=Malha.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//        sprintf(NomeArquivo, "ArquivosPlot/Interface/PontosGnuplot%d.txt", i++);
//        arq = fopen(NomeArquivo, "wt");
//
//        //Percorre cada ponto desta curva e imprime no arquivo
//        for(ponto = curva->nodeInicial; ponto!=NULL; ponto=ponto->prox ) {
//            //Descomenta essa linha pra plotar vetor normal
//            //fprintf(arq, "set arrow 10 from %lf, %lf to %lf, %lf\n", ponto->x, ponto->y, ponto->vetNormalX, ponto->vetNormalY);
//            //printf("%e %e %e %e\n", ponto->x, ponto->y, ponto->vetNormalX, ponto->vetNormalY);
//
//            //VetorNormalMinimosQuadrados(Malha.interface.curvaInicial, ponto->x, ponto->y, Malha.minDxDy, &(ponto->vetNormalX), &(ponto->vetNormalY));
//
//            //if( i%20==0 )
//            fprintf(arq, "%.15e %.15e %.15e %.15e\n", ponto->x, ponto->y, 0.05*ponto->vetNormalX, 0.05*ponto->vetNormalY);
//        }
//
//        fclose(arq);
//    }
//
//    //Fazendo o script pra plotar
//    arq = fopen("ArquivosPlot/Interface/ScriptGnuplot.txt", "wt");
//
//    //Apenas se quiser imprimir no PDF
//    fprintf(arq, "set terminal pdf\n");
//    fprintf(arq, "set output 'ArquivosPlot/Interface/Interface.pdf'\n");
//
//
//    fprintf(arq, "set xrange[%lf:%lf]\n", Malha.xMin, Malha.xMax);
//    fprintf(arq, "set yrange[%lf:%lf]\n", Malha.yMin, Malha.yMax);
//    fprintf(arq, "set size ratio -1\n");
//
//    ///Desenhando a malha (se quiser)
//    if( DesenhaMalha ) {
//        x = Malha.xMin;
//        for( i=0; i<=Malha.Nx; i++ ) {
//            fprintf(arq, "set arrow from %lf, %lf to %lf, %lf nohead\n", x, Malha.yMin, x, Malha.yMax);
//            if( i!=Malha.Nx )
//                x += Malha.dx[i];
//        }
//
//        y = Malha.yMin;
//        for( j=0; j<=Malha.Ny; j++ ) {
//            fprintf(arq, "set arrow from %lf, %lf to %lf, %lf nohead\n", Malha.xMin, y, Malha.xMax, y);
//            if( j!=Malha.Ny )
//                y += Malha.dy[j];
//        }
//    }
//    ///Fim da malha
//
//    i = 0;
//    for( curva=Malha.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//        if( i==0 ) {
//            //sprintf(NomeArquivo, "plot \"PontosGnuplot%d.txt\" using 1:2 with lines notitle lc 0\n", i);
//            //InsereComandoGnuplot(NomeArquivo);
//            fprintf(arq, "plot \"ArquivosPlot/Interface/PontosGnuplot%d.txt\" using 1:2 with lines notitle\n", i);
//            if( DesenhaNormais )
//                fprintf(arq, "replot \"ArquivosPlot/Interface/PontosGnuplot%d.txt\" using 1:2:3:4 with vectors notitle\n", i);
//        }
//        else {
//            //sprintf(NomeArquivo, "replot \"PontosGnuplot%d.txt\" using 1:2 with lines notitle lc 0\n", i);
//            //InsereComandoGnuplot(NomeArquivo);
//            fprintf(arq, "replot \"ArquivosPlot/Interface/PontosGnuplot%d.txt\" using 1:2 with lines notitle\n", i);
//            if( DesenhaNormais )
//                fprintf(arq, "replot \"ArquivosPlot/Interface/PontosGnuplot%d.txt\" using 1:2:3:4 with vectors notitle\n", i);
//        }
//        i++;
//    }
//
//    fprintf(arq, "unset output\n");
//    fprintf(arq, "set output 'ArquivosPlot/Interface/Interface%d.pdf'\n", n);
//    fprintf(arq, "replot\n");
//    fprintf(arq, "unset output\n");
//    fclose(arq);
//
//
//
//    //Chama o gnuplot pra plotar
//    system("gnuplot ArquivosPlot/Interface/ScriptGnuplot.txt");

    return;
}

void ImprimeArquivoColetaVTK(float **flowType, float **magnitudeVelocidade, float **energiaDissipativaCelulas, MALHA Malha, int n) {
    FILE *arq;
    int i, j;
    char nomeArq[500];

    //Criando os arquivos
    sprintf(nomeArq, "ArquivosPlot/VTK/%s/Coleta-N%d.vtk", Malha.nomeArquivo, n);
    arq = fopen(nomeArq, "wt");
    if( !arq )
        DeuProblema("\n\n ImprimeArquivoVTK: Problema no fopen do arquivo. \n\n");

    //Cabecalho
    fprintf(arq, "# vtk DataFile Version 2.0\nStructured Grid Example\nASCII\n");
    fprintf(arq, "DATASET RECTILINEAR_GRID\n");
    fprintf(arq, "DIMENSIONS %d %d 1\n", Malha.Nx+1, Malha.Ny+1);
    fprintf(arq, "X_COORDINATES %d float\n", Malha.Nx+1);
    for( i=0; i<=Malha.Nx; i++ )
        fprintf(arq, "%lf ", Malha.x[i]);
    fprintf(arq, "\nY_COORDINATES %d float\n", Malha.Ny+1);
    for( i=0; i<=Malha.Ny; i++ )
        fprintf(arq, "%lf ", Malha.y[i]);
    fprintf(arq, "\nZ_COORDINATES %d float\n", 1);
    fprintf(arq, "0.0\n\n");  
    
    fprintf(arq, "SCALARS flowType float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Malha.Ny; j++)
        for(i=0; i<=Malha.Nx; i++)
            fprintf(arq, "%e\n", flowType[i][j]);
    fprintf(arq, "\n");

    fprintf(arq, "SCALARS magnitude float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Malha.Ny; j++)
        for(i=0; i<=Malha.Nx; i++)
            fprintf(arq, "%e\n", magnitudeVelocidade[i][j]);
    fprintf(arq, "\n");

    fprintf(arq, "SCALARS dissipativa float 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for(j=0; j<=Malha.Ny; j++)
        for(i=0; i<=Malha.Nx; i++)
            fprintf(arq, "%e\n", energiaDissipativaCelulas[i][j]);
    fprintf(arq, "\n");

    fclose(arq);

    return;
}


