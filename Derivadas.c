#include "Derivadas.h"

extern int velocidadeNewtoniana;


//*   Operadores     */
//*   Operadores     */
double Lx(double **U, int i, int j, double t, MALHA M)
{
    double uEsq = 0.0, uDir = 0.0, uBaixo = 0.0, uCima = 0.0, uCentro = 0.0;
    int iV;

    return 0.0;

   if( i==M.Nx )
        iV = i-1;
    else
        iV = i;

    //Verificando o U(i, j-1)
    if( BAIXO_U(iV, j) ) {
        if( M.pontosV[iV][j].tipoBoundary == NOSLIP )
            uBaixo = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j].tipoBoundary == SLIP )
            uBaixo = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j].tipoBoundary == INFLOW )
            uBaixo = - U[i][j];
        else if( M.pontosV[iV][j].tipoBoundary == NEUMANN )
            uBaixo = U[i][j];
        else if( M.pontosV[iV][j].tipoBoundary == SIMETRIA )
            uBaixo = U[i][j];
        else
            DeuProblema("Lx: problema baixo %d %d\n", i, j);
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        uBaixo = U[i][j-1];

    //Encontrando o U(i, j+1)
    if( CIMA_U(i, j) ) {
        if( (M.pontosV[iV][j+1].tipoBoundary == NOSLIP) )
            uCima = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( (M.pontosV[iV][j+1].tipoBoundary == SLIP) )
            uCima = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j+1].tipoBoundary == INFLOW )
            uCima = - U[i][j];
        else if( M.pontosV[iV][j+1].tipoBoundary == NEUMANN )
            uCima = U[i][j];
        else if( M.pontosV[iV][j+1].tipoBoundary == SIMETRIA )
            uCima = U[i][j];
        else
            DeuProblema("Lx: problema cima\n");
    }
    else
        uCima = U[i][j+1];

    //Encontrando o U(i-1, j)
    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
        DeuProblema("Lx: Nao era pra entrar aqui");
        if( M.pontosU[i][j].tipoBoundary==NEUMANN ) {
            uEsq = U[i][j]; //Ordem 1
        }
        else
            DeuProblema("Lx: problema esquerda\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else {
        uEsq = U[i-1][j];
    }


    //Encontrando o U(i+1, j)
    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
        DeuProblema("Lx: Nao era pra entrar aqui");
        if( M.pontosU[i][j].tipoBoundary==NEUMANN ) {
            uDir = U[i][j]; //Ordem 1
        }
        else
            DeuProblema("Lx: problema direita\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else {
        uDir = U[i+1][j];
    }

    if( (i==0) || (i==M.Nx) )
        DeuProblema("\n\n Lx: NAO ERA PRA TER ENTRADO AQUI!! \n\n");


    double h1, h2;

    h1 = M.dx[i-1];
    h2 = M.dx[i];
    double parteX = (2.0/(h1*(h1+h2)))*uEsq - (2.0/(h1*h2))*uCentro + (2.0/(h2*(h1+h2)))*uDir;

    h1 = (j==0) ? M.dy[j] : 0.5*(M.dy[j-1] + M.dy[j]);
    h2 = (j==M.Ny-1) ? M.dy[j] : 0.5*(M.dy[j+1] + M.dy[j]);

    double parteY = (2.0/(h1*(h1+h2)))*uBaixo - (2.0/(h1*h2))*uCentro + (2.0/(h2*(h1+h2)))*uCima;

    return parteX + parteY;
}

double Ly(double **V, int i, int j, MALHA M)
{
    static double vEsq = 0.0, vDir = 0.0, vBaixo = 0.0, vCima = 0.0;
    int jU;

    return 0.0;

    if( j==M.Ny )
        jU = j-1;
    else
        jU = j;

    //Encontrando o V(i-1, j)
    if( ESQUERDA_V(i, j) ) {
        if( M.pontosU[i][jU].tipoBoundary == NOSLIP )
            vEsq = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i][jU].tipoBoundary == SLIP )
            vEsq = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i][jU].tipoBoundary == INFLOW )
            vEsq = - V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == NEUMANN )
            vEsq = V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == SIMETRIA )
            vEsq = V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == EIXO_AXISSIMETRIA )
            vEsq = V[i][j];
        else
            DeuProblema("Ly: problema esquerda\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        vEsq = V[i-1][j];

    //Encontrando o V(i+1, j)
    if( DIREITA_V(i, j) ) {
        if( M.pontosU[i+1][jU].tipoBoundary==NEUMANN )
            vDir = V[i][j];
        else if( M.pontosU[i+1][jU].tipoBoundary==SIMETRIA )
            vDir = V[i][j];
        else if( M.pontosU[i+1][jU].tipoBoundary==NOSLIP )
            vDir = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][jU].tipoBoundary==SLIP )
            vDir = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][jU].tipoBoundary==INFLOW )
            vDir = - V[i][j];
        else
            DeuProblema("Ly: problema direita\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        vDir = V[i+1][j];

    //Encontrando o V(i, j+1)
    if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].extremo==CIMA) {
        printf("\nDO JEITO QUE TA AGORA NAO ERA PRA ENTRAR AQUI\n");
        if( M.pontosV[i][j].tipoBoundary==NEUMANN )
            vCima = V[i][j-1]; //Neumann ordem 2
        else
            DeuProblema("Ly: problema cima\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        vCima = V[i][j+1];

    //Encontrando o V(i, j-1)
    if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].extremo==BAIXO ) {
        printf("\nDO JEITO QUE TA AGORA NAO ERA PRA ENTRAR AQUI\n");
        if( M.pontosV[i][j].tipoBoundary==NEUMANN )
            vBaixo = V[i][j-1]; //Neumann ordem 2
        else
            DeuProblema("Ly: problema baixo\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        vBaixo = V[i][j+1];

    /// ========= Precisa mudar isto aqui pra malha nao-uniforme
    return ((vDir - 2*V[i][j] + vEsq)/(M.dx[0]*M.dx[0])) +  (vCima - 2*V[i][j] + vBaixo)/(M.dy[0]*M.dy[0]) ;
}

double Cx(double **U, double **V, double **W, int i, int j, double t, MALHA M)
{
    static double uBaixo = 0.0, uCima = 0.0, uEsq = 0.0, uDir = 0.0;
    static double wDir = 0.0, wEsq = 0.0, wCentro = 0.0;
    static double vBaixoEsq = 0.0, vBaixoDir = 0.0, vCimaEsq = 0.0, vCimaDir = 0.0;
    double rImenosMeio;
    double h1, h2;
    char boundaryEsq, boundaryDir, boundaryBaixo, boundaryCima;
    int iV;

//    return 0.0;

    if( i==M.Nx )
        iV = i-1;
    else
        iV = i;

    //Verificando o U(i, j-1)
    if( BAIXO_U(iV, j) ) {
        if( M.pontosV[iV][j].tipoBoundary == NOSLIP )
            uBaixo = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j].tipoBoundary == SLIP )
            uBaixo = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j].tipoBoundary == INFLOW )
            uBaixo = - U[i][j];
        else if( M.pontosV[iV][j].tipoBoundary == NEUMANN )
            uBaixo = U[i][j];
        else if( M.pontosV[iV][j].tipoBoundary == SIMETRIA )
            uBaixo = U[i][j];
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("Cx: problema baixo %d %d %d\n", i, j, M.pontosV[iV][j].tipoBoundary);
        }
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        uBaixo = U[i][j-1];

    //Encontrando o U(i, j+1)
    if( CIMA_U(i, j) ) {
        if( (M.pontosV[iV][j+1].tipoBoundary == NOSLIP) )
            uCima = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( (M.pontosV[iV][j+1].tipoBoundary == SLIP) )
            uCima = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j+1].tipoBoundary == INFLOW )
            uCima = - U[i][j];
        else if( M.pontosV[iV][j+1].tipoBoundary == NEUMANN )
            uCima = U[i][j];
        else if( M.pontosV[iV][j+1].tipoBoundary == SIMETRIA )
            uCima = U[i][j];
        else if( M.pontosU[i][j+1].tipoBoundary==PERIODICIDADE )
            uCima = 0.0; //period mudar isso?
        else
            DeuProblema("Cx: problema cima %d %d\n", i, j);
    }
    else
        uCima = U[i][j+1];

    //Encontrando o U(i-1, j)
    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
//        if( M.pontosU[i][j].tipoBoundary==NEUMANN ) {
//            uEsq = U[i][j]; //Ordem 1
//            vBaixoEsq = V[i][j]; //Neumann
//            vCimaEsq = V[i][j+1]; //Neumann
//            wEsq = (M.tipoCoord==AXI_CILINDRICO) ? W[i][j] : 0.0;
//        }
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE ) {
            uEsq = U[M.Nx-1][j];
            vBaixoEsq = V[M.Nx-1][j];
            vCimaEsq = V[M.Nx-1][j+1];
            wEsq = (M.tipoCoord==AXI_CILINDRICO) ? W[M.Nx-1][j] : 0.0;
        }
        else
            DeuProblema("Cx: problema esquerda\n\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else {
        uEsq = U[i-1][j];
        vBaixoEsq = V[i-1][j];
        vCimaEsq = V[i-1][j+1];
        wEsq = (M.tipoCoord==AXI_CILINDRICO) ? W[i-1][j] : 0.0;
    }


    //Encontrando o U(i+1, j)
    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
//        if( M.pontosU[i][j].tipoBoundary==NEUMANN ) {
//            uDir = U[i][j]; //Ordem 1
//            vBaixoDir = V[i-1][j]; //Neumann
//            vCimaDir = V[i-1][j+1]; //Neumann
//            wEsq = (M.tipoCoord==AXI_CILINDRICO) ? W[i-1][j] : 0.0;
//        }
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE ) {
            uDir = U[1][j];
            vBaixoDir = V[0][j];
            vCimaDir = V[0][j+1];
            wDir = (M.tipoCoord==AXI_CILINDRICO) ? W[0][j] : 0.0;
        }
        else
            DeuProblema("Cx: problema direita\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else {
        uDir = U[i+1][j];
        vBaixoDir = V[i][j];
        vCimaDir = V[i][j+1];
        wDir = (M.tipoCoord==AXI_CILINDRICO) ? W[i][j] : 0.0;
    }

    if( M.tipoCoord == AXI_CILINDRICO )
        rImenosMeio = M.r[i] - 0.5*M.dx[i];
    else
        rImenosMeio = 1.0;

    if( (i==0) || (i==M.Nx) ) {
        if( M.pontosU[i][j].tipoBoundary!=PERIODICIDADE )
            DeuProblema("\n\n Cx: NAO ERA PRA TER ENTRADO AQUI!! \n\n");
    }


    boundaryDir = boundaryEsq = boundaryBaixo = boundaryCima = 0;
    if( (i<2) || M.pontosU[i][j].extremo==ESQUERDA || M.pontosU[i-1][j].extremo==ESQUERDA
        || M.pontosU[i-1][j].tipo==SURFACE || M.pontosU[i-2][j].tipo==SURFACE)
        boundaryEsq = 1;
    if( (i>M.Nx-2) || M.pontosU[i][j].extremo==DIREITA || M.pontosU[i+1][j].extremo==DIREITA
       || M.pontosU[i+1][j].tipo==SURFACE || M.pontosU[i+2][j].tipo==SURFACE )
        boundaryDir = 1;
    if( (j<2) || M.pontosV[iV][j].tipo==BOUNDARY || M.pontosV[iV][j-1].tipo==BOUNDARY
       || M.pontosU[i][j-1].tipo==SURFACE || M.pontosU[i][j-2].tipo==SURFACE )
        boundaryBaixo = 1;
    if( (j>M.Ny-3) || M.pontosV[iV][j+1].tipo==BOUNDARY || M.pontosV[iV][j+2].tipo==BOUNDARY
       || M.pontosU[i][j+1].tipo==SURFACE || M.pontosU[i][j+2].tipo==SURFACE )
        boundaryCima = 1;

    /* CUBISTA */
    double fluxoUU, fluxoUV;
    //Aproximando a derivada UU
    double vel, phiR, phiU, phiD, parte1, parte2;
    double coordR, coordU, coordD, coordF;
    vel = (i==M.Nx) ? M.r[i-1]*0.5*(uDir + U[i][j]) : M.r[i]*0.5*(uDir + U[i][j]);
    coordF = (i==M.Nx) ? M.x[i] + 0.5*M.dx[i-1] : M.x[i] + 0.5*M.dx[i];
    if( vel>=0 ) {
        phiR = uEsq;
        phiU = U[i][j];
        phiD = uDir;

        coordR = (i==0) ? M.x[0] - M.dx[0] : M.x[i-1];
        coordU = M.x[i];
        coordD = (i==M.Nx) ? M.x[i] + M.dx[i-1] : M.x[i+1];
    }
    else {
        phiR = boundaryDir ? 0.0 : U[i+2][j];
        phiU = uDir;
        phiD = U[i][j];

        coordR = boundaryDir ? 0.0 : M.x[i+2];
        coordU = (i==M.Nx) ? M.x[i] + M.dx[i-1] : M.x[i+1];
        coordD = M.x[i];
    }
    parte1 = boundaryDir ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    vel = (i==0) ? (M.r[i] - M.dx[0])*0.5*(uEsq + U[i][j]) :  M.r[i-1]*0.5*(uEsq + U[i][j]);
    coordF = (i==0) ? M.x[i] - 0.5*M.dx[i] : M.x[i-1] + 0.5*M.dx[i-1];
    if( vel>=0 ) {
        phiR = boundaryEsq ? 0.0 : U[i-2][j];
        phiU = uEsq;
        phiD = U[i][j];

        coordR = boundaryEsq ? 0.0 : M.x[i-2];
        coordU = (i==0) ? M.x[i] - M.dx[i] : M.x[i-1];
        coordD = M.x[i];

    }
    else {
        phiR = uDir;
        phiU = U[i][j];
        phiD = uEsq;

        coordR = (i==M.Nx) ? M.x[i] + M.dx[i-1] : M.x[i+1];
        coordU = M.x[i];
        coordD = (i==0) ? M.x[i] - M.dx[i] : M.x[i-1];
    }
    parte2 = boundaryEsq ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    h1 = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
    h2 = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];
    fluxoUU = (1.0/rImenosMeio)*( (-parte2*h2/(h1*(h1+h2))) + (rImenosMeio*QUAD(U[i][j])*(h2-h1)/(h1*h2)) + (parte1*h1/(h2*(h1+h2))) );







    //Aproximando a derivada UV
    double dxEsq = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
    double dxDir = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];
    vel = Interpolacao(dxEsq, dxDir, vCimaEsq, vCimaDir);
    coordF = M.y[j+1];
    if( vel>=0 ) {
        phiR = uBaixo;
        phiU = U[i][j];
        phiD = uCima;

        coordR = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
        coordU = M.y[j] + 0.5*M.dy[j];
        coordD = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
    }
    else {
        phiR = boundaryCima ? 0.0 : U[i][j+2];
        phiU = uCima;
        phiD = U[i][j];

        coordR = boundaryCima ? 0.0 : M.y[j+2] + 0.5*M.dy[j+2];
        coordU = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
        coordD = M.y[j] + 0.5*M.dy[j];
    }
    parte1 = boundaryCima ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    vel = Interpolacao(dxEsq, dxDir, vBaixoEsq, vBaixoDir);
    coordF = M.y[j];
    if( vel>=0 ) {
        phiR = boundaryBaixo ? 0.0 : U[i][j-2];
        phiU = uBaixo;
        phiD = U[i][j];

        coordR = boundaryBaixo ? 0.0 : M.y[j-2] + 0.5*M.dy[j-2];
        coordU = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
        coordD = M.y[j] + 0.5*M.dy[j];
    }
    else {
        phiR = uCima;
        phiU = U[i][j];
        phiD = uBaixo;

        coordR = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
        coordU = M.y[j] + 0.5*M.dy[j];
        coordD = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
    }
    parte2 = boundaryBaixo ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    fluxoUV = (parte1 - parte2)/M.dy[j];

    wCentro = Interpolacao(dxEsq, dxDir, wEsq, wDir);
    return fluxoUU + fluxoUV - (QUAD(wCentro)/rImenosMeio);
}

double Cy(double **U, double **V, int i, int j, MALHA M)
{
    static double vEsq = 0.0, vDir = 0.0, vBaixo = 0.0, vCima = 0.0;
    static double uBaixoEsq = 0.0, uCimaEsq = 0.0, uBaixoDir = 0.0, uCimaDir = 0.0;
    double rImaisMeio, rImenosMeio;
    double h1, h2;
    char boundaryEsq, boundaryDir, boundaryBaixo, boundaryCima;
    int jU;

//    return 0.0;

    if( j==M.Ny )
        jU = j-1;
    else
        jU = j;

    //Encontrando o V(i-1, j)
    if( ESQUERDA_V(i, j) ) {
        if( M.pontosU[i][jU].tipoBoundary == NOSLIP )
            vEsq = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i][jU].tipoBoundary == SLIP )
            vEsq = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i][jU].tipoBoundary == INFLOW )
            vEsq = - V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == NEUMANN )
            vEsq = V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == SIMETRIA )
            vEsq = V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == PERIODICIDADE )
            vEsq = V[M.Nx-1][j];
        else if( M.pontosU[i][jU].tipoBoundary == EIXO_AXISSIMETRIA )
            vEsq = V[i][j];
        else
            DeuProblema("Derivadas Cy: Problema esquerda\n\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        vEsq = V[i-1][j];

    //Encontrando o V(i+1, j)
    if( DIREITA_V(i, j) ) {
        if( M.pontosU[i+1][jU].tipoBoundary==NEUMANN )
            vDir = V[i][j];
        else if( M.pontosU[i+1][jU].tipoBoundary==SIMETRIA )
            vDir = V[i][j];
        else if( M.pontosU[i+1][jU].tipoBoundary==NOSLIP )
            vDir = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][jU].tipoBoundary==SLIP )
            vDir = V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][jU].tipoBoundary==INFLOW )
            vDir = - V[i][j];
        else if( M.pontosU[i+1][jU].tipoBoundary==PERIODICIDADE )
            vDir = V[0][j];
        else
            DeuProblema("Derivadas Cy: Problema direita\n\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        vDir = V[i+1][j];

    //Encontrando o V(i, j-1)
    if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].extremo==BAIXO )
        DeuProblema("Cy: Nao era pra ter entrado aqui %d %d!!! \n\n", i, j);
    else {
        vBaixo = V[i][j-1]; //Neumann ordem 2
        uBaixoEsq = U[i][j-1];
        uBaixoDir = U[i+1][j-1];
    }

    //Encontrando o V(i, j+1)
    if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].extremo==CIMA )
        DeuProblema("Cy: Nao era pra ter entrado aqui %d %d!!! \n\n", i, j);
    else {
        vCima = V[i][j+1]; //Neumann ordem 2
        uCimaEsq = U[i][j];
        uCimaDir = U[i+1][j];
    }

    if( M.tipoCoord==AXI_CILINDRICO ) {
        rImaisMeio = M.r[i] + 0.5*M.dx[i];
        rImenosMeio = M.r[i] - 0.5*M.dx[i];
    }
    else
        rImaisMeio = rImenosMeio = 1.0;

 /*   if( (j<2) || (j>M.Ny-2) || (i<2) || (i>M.Nx-3) ||
       M.pontosV[i][j].tipo==BOUNDARY || M.pontosV[i][j-1].tipo==BOUNDARY || M.pontosV[i][j+1].tipo==BOUNDARY ||
       M.pontosU[i][j].tipo==BOUNDARY || M.pontosU[i-1][j].tipo==BOUNDARY || M.pontosU[i+1][j].tipo==BOUNDARY )
        return (1.0/(4*M.dy[0]))*( QUAD(vCima+V[i][j])-QUAD(V[i][j]+vBaixo) )
            + (1.0/(4*M.r[i]*M.dx[0]))*rImaisMeio*(uCimaDir+uBaixoDir)*(vDir+V[i][j]) - (1.0/(4*M.r[i]*M.dx[0]))*rImenosMeio*(uCimaEsq+uBaixoEsq)*(V[i][j]+vEsq);
*/










    boundaryDir = boundaryEsq = boundaryBaixo = boundaryCima = 0;
    if( (i<2) || M.pontosU[i][j].tipo==BOUNDARY || M.pontosU[i-1][j].tipo==BOUNDARY
       || M.pontosV[i-1][j].tipo==SURFACE || M.pontosV[i-2][j].tipo==SURFACE )
        boundaryEsq = 1;
    if( (i>M.Nx-3) || M.pontosU[i+1][j].tipo==BOUNDARY || M.pontosU[i+2][j].tipo==BOUNDARY
       || M.pontosV[i+1][j].tipo==SURFACE || M.pontosV[i+2][j].tipo==SURFACE )
        boundaryDir = 1;
    if( (j<2) || M.pontosV[i][j].extremo==BAIXO || M.pontosV[i][j-1].extremo==BAIXO
       || M.pontosV[i][j-1].tipo==SURFACE || M.pontosV[i][j-2].tipo==SURFACE )
        boundaryBaixo = 1;
    if( (j>M.Ny-2) || M.pontosV[i][j].extremo==CIMA || M.pontosV[i][j+1].extremo==CIMA
       || M.pontosV[i][j+1].tipo==SURFACE || M.pontosV[i][j+2].tipo==SURFACE )
        boundaryCima = 1;

    if( (j==0) || (j==M.Ny) ) {
        DeuProblema("\n\n Cy: NAO ERA PRA TER ACONTECIDO ISSO Cy %d %d\n\n", i, j);
    }

    /* CALCULANDO OS PONTOS R, U, D  DO ESQUEMA CUBISTA (OU OUTRO ESQUEMA, TANTO FAZ)*/
    double fluxoVV, fluxoUV;
    //Aproximando a derivada UV
    double vel, phiR, phiU, phiD, parte1, parte2;
    double coordF, coordR, coordU, coordD;
    vel = rImaisMeio*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], uBaixoDir, uCimaDir);
    coordF = M.x[i+1];
    if( vel>=0 ) {
        phiR = vEsq;
        phiU = V[i][j];
        phiD = vDir;

        coordR = (i==0) ? M.x[i] - 0.5*M.dx[i] : M.x[i-1] + 0.5*M.dx[i-1];
        coordU = M.x[i] + 0.5*M.dx[i];
        coordD = (i==M.Nx-1) ? M.x[i+1] + 0.5*M.dx[i] : M.x[i+1] + 0.5*M.dx[i+1];
    }
    else {
        phiR = boundaryDir ? 0.0 : V[i+2][j];
        phiU = vDir;
        phiD = V[i][j];



        coordR = boundaryDir ? 0.0 : M.x[i+2] + 0.5*M.dx[i+2];
        coordU = (i==M.Nx-1) ? M.x[i+1] + 0.5*M.dx[i] : M.x[i+1] + 0.5*M.dx[i+1];
        coordD = M.x[i] + 0.5*M.dx[i];
    }
    parte1 = boundaryDir ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);



    vel = rImenosMeio*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], uBaixoEsq, uCimaEsq);
    coordF = M.x[i];
    if( vel>=0 ) {
        phiR = boundaryEsq ? 0.0 : V[i-2][j];
        phiU = vEsq;
        phiD = V[i][j];

        coordR = boundaryEsq ? 0.0 : M.x[i-2] + 0.5*M.dx[i-2];
        coordU = (i==0) ? M.x[i] - 0.5*M.dx[i] : M.x[i-1] + 0.5*M.dx[i-1];
        coordD = M.x[i] + 0.5*M.dx[i];
    }
    else {
        phiR = vDir;
        phiU = V[i][j];
        phiD = vEsq;

        coordR = (i==M.Nx-1) ? M.x[i+1] + 0.5*M.dx[i] : M.x[i+1] + 0.5*M.dx[i+1];
        coordU = M.x[i] + 0.5*M.dx[i];
        coordD = (i==0) ? M.x[i] - 0.5*M.dx[i] : M.x[i-1] + 0.5*M.dx[i-1];
    }
    parte2 = boundaryEsq ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    fluxoUV = (1.0/M.r[i])*(parte1 - parte2)/M.dx[i];

    //Aproximando a derivada VV
    vel = 0.5*(V[i][j] + vCima);
    coordF = M.y[j] + 0.5*M.dy[j];
    if( vel>=0 ) {
        phiR = vBaixo;
        phiU = V[i][j];
        phiD = vCima;

        coordR = boundaryBaixo ? 0.0 : M.y[j-1];
        coordU = M.y[j];
        coordD = (j==M.Ny) ? M.y[j] + M.dy[j-1] : M.y[j+1];
    }
    else {
        phiR = boundaryCima ? 0.0 : V[i][j+2];
        phiU = vCima;
        phiD = V[i][j];

        coordR = boundaryCima ? 0.0 : M.y[j+2];
        coordU = (j==M.Ny) ? M.y[j] + M.dy[j-1] : M.y[j+1];
        coordD = M.y[j];
    }
    parte1 = boundaryCima ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);


    vel = 0.5*(V[i][j] + vBaixo);
    coordF = M.y[j-1] + 0.5*M.dy[j-1];
    if( vel>=0 ) {
        phiR = boundaryBaixo ? 0.0 : V[i][j-2];
        phiU = vBaixo;
        phiD = V[i][j];

        coordR = boundaryBaixo ? 0.0 : M.y[j-2];
        coordU = (j==0) ? M.y[j] - M.dy[j] : M.y[j-1];
        coordD = M.y[j];
    }
    else {
        phiR = vCima;
        phiU = V[i][j];
        phiD = vBaixo;

        coordR = (j==M.Ny) ? M.y[j] + M.dy[j-1] : M.y[j+1];
        coordU = M.y[j];
        coordD = (j==0) ? M.y[j] - M.dy[j] : M.y[j-1];
    }
    parte2 = boundaryBaixo ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    h1 = (j==0) ? 0.5*M.dy[j] : 0.5*M.dy[j-1];
    h2 = 0.5*M.dy[j];
    fluxoVV = (-parte2*h2/(h1*(h1+h2))) + (QUAD(V[i][j])*(h2-h1)/(h1*h2)) + (parte1*h1/(h2*(h1+h2)));


    return fluxoVV + fluxoUV;
}

double Ct(double **U, double **V, double **W, int i, int j, MALHA M)
{
    double wBaixo, wCima, wEsq, wDir;
    double coordF, coordR, coordU, coordD;
    double valorDirichlet;

    wBaixo = wCima = wEsq = wDir = 1e+10; //inicializando com LIXO. Esse valor tem que ser substituido abaixo

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
        else
            printf("\n\n Problema Esquerda - Convectivo Theta\n\n");
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
        else
            printf("\n\n Problema Direita - Convectivo Theta\n\n");
    }

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
        else
            printf("problema Baixo - Convectivo Theta");
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
        else
            printf("problema Cima - Convectivo Theta");
    }

//    if( (j<2) || (j>M.Ny-3) || (i<2) || (i>M.Nx-3) ||
//       M.pontosW[i-1][j].tipo!=FULL || M.pontosW[i-2][j].tipo!=FULL || M.pontosW[i+1][j].tipo!=FULL || M.pontosW[i+2][j].tipo!=FULL ||
//       M.pontosW[i][j-1].tipo!=FULL || M.pontosW[i][j-2].tipo!=FULL || M.pontosW[i][j+1].tipo!=FULL || M.pontosW[i][j+2].tipo!=FULL )
//        return (1.0/(2.0*M.dx[i]*M.r[i]))*( (M.r[i]+0.5*M.dx[i])*U[i+1][j]*(wDir + W[i][j]) + (M.r[i]-0.5*M.dx[i])*U[i][j]*(W[i][j] + wEsq) )
//            + 0.5*(1.0/M.dy[j])*( (wCima + W[i][j])*V[i][j+1] - (W[i][j] + wBaixo)*V[i][j] )
//            + ((U[i+1][j] + U[i][j])*W[i][j])/(2.0*M.r[i]);



    int boundaryDir, boundaryEsq, boundaryBaixo, boundaryCima = 0;
    boundaryDir = boundaryEsq = boundaryBaixo = boundaryCima = 0;
    if( (i<2) || M.celulas[i-1][j]==EMPTY || M.celulas[i-2][j]==EMPTY )
        boundaryEsq = 1;
    if( (i>M.Nx-3) || M.celulas[i+1][j]==EMPTY || M.celulas[i+2][j]==EMPTY )
        boundaryDir = 1;
    if( (j<2) || M.celulas[i][j-1]==EMPTY || M.celulas[i][j-2]==EMPTY )
        boundaryBaixo = 1;
    if( (j>M.Ny-3) || M.celulas[i][j+1]==EMPTY || M.celulas[i][j+2]==EMPTY )
        boundaryCima = 1;

    double rImaisMeio = M.r[i] + 0.5*M.dx[i];
    double rImenosMeio = M.r[i] - 0.5*M.dx[i];

    double lixo = 1e+10;

    //USANDO CUBISTA
    double fluxoUW, fluxoVW;
    double vel, phiR, phiU, phiD, parte1, parte2;

    //Aproximando a derivada UW
    vel = rImaisMeio*rImaisMeio*U[i+1][j];
    coordF = M.x[i+1];
    if( vel>=0 ) {
        phiR = wEsq; //W[i-1][j];
        phiU = W[i][j];
        phiD = wDir; //W[i+1][j];

        coordR = (i==0) ? M.r[i] - M.dx[i] : M.r[i-1];
        coordU = M.r[i];
        coordD = (i==M.Nx-1) ? M.r[i] + M.dx[i] : M.r[i+1];
    }
    else {
        phiR = boundaryDir ? lixo : W[i+2][j]; //Nao eh pra usar phiR pra nada se for boundaryDir
        phiU = wDir; //W[i+1][j];
        phiD = W[i][j];

        coordR = boundaryDir ? lixo : M.r[i+2];
        coordU = (i==M.Nx-1) ? M.r[i] + M.dx[i] : M.r[i+1];
        coordD = M.r[i];
    }
    parte1 = boundaryDir ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    vel = rImenosMeio*rImenosMeio*U[i][j];
    coordF = M.x[i];
    if( vel>=0 ) {
        phiR = boundaryEsq ? lixo : W[i-2][j];
        phiU = wEsq; //W[i-1][j];
        phiD = W[i][j];

        coordR = boundaryEsq ? lixo : M.r[i-2];
        coordU = (i==0) ? M.r[i] - M.dx[i] : M.r[i-1];
        coordD = M.r[i];
    }
    else {
        phiR = wDir; //W[i+1][j];
        phiU = W[i][j];
        phiD = wEsq; //W[i-1][j];

        coordR = (i==M.Nx-1) ? M.r[i] + M.dx[i] : M.r[i+1];
        coordU = M.r[i];
        coordD = (i==0) ? M.r[i] - M.dx[i] : M.r[i-1];
    }
    parte2 = boundaryEsq ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    fluxoUW = (1.0/(M.r[i]*M.r[i]))*(parte1 - parte2)/M.dx[i];


    //Aproximando a derivada VW
    vel = V[i][j+1];
    coordF = M.y[j+1];
    if( vel>=0 ) {
        phiR = wBaixo; //W[i][j-1];
        phiU = W[i][j];
        phiD = wCima; //W[i][j+1];

        coordR = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
        coordU = M.y[j] + 0.5*M.dy[j];
        coordD = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
    }
    else {
        phiR = boundaryCima ? lixo : W[i][j+2];
        phiU = wCima; //W[i][j+1];
        phiD = W[i][j];

        coordR = boundaryCima ? lixo : M.y[j+2] + 0.5*M.dy[j+2];
        coordU = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
        coordD = M.y[j] + 0.5*M.dy[j];
    }
    parte1 = boundaryCima ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    vel = V[i][j];
    coordF = M.y[j];
    if( vel>=0 ) {
        phiR = boundaryBaixo ? lixo : W[i][j-2];
        phiU = wBaixo; //W[i][j-1];
        phiD = W[i][j];

        coordR = boundaryBaixo ? lixo : M.y[j-2] + 0.5*M.dy[j-2];
        coordU = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
        coordD = M.y[j] + 0.5*M.dy[j];
    }
    else {
        phiR = wCima; //W[i][j+1];
        phiU = W[i][j];
        phiD = wBaixo; //W[i][j-1];

        coordR = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
        coordU = M.y[j] + 0.5*M.dy[j];
        coordD = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
    }
    parte2 = boundaryBaixo ? vel*FuncaoConvCentral(phiU, phiD, coordF, coordU, coordD) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);

    fluxoVW = (parte1 - parte2)/M.dy[j];

    return fluxoUW + fluxoVW;
}

double GradPy(double **P, int i, int j, MALHA M)
{
    double pBaixo, pCima, pCentro;
    double h1, h2;

    h1 = (j==0) ? 0.5*M.dy[j] : 0.5*M.dy[j-1];
    h2 = (j==M.Ny) ? 0.5*M.dy[j-1] : 0.5*M.dy[j];

    pCima = (j==M.Ny) ? 0.0 : P[i][j];
    pBaixo = (j==0) ? 0.0 : P[i][j-1];
    pCentro = Interpolacao(h1, h2, pBaixo, pCima);

    if( (j==0) || (j==M.Ny) )
        DeuProblema("GradPy: problema no indice (%d %d)\n", i, j);

    return (-pBaixo*h2/(h1*(h1+h2))) + (pCentro*(h2-h1)/(h1*h2)) + (pCima*h1/(h2*(h1+h2)));
}

double GradPx(double **P, int i, int j, MALHA M)
{
    double pEsq, pDir, pCentro;
    double h1, h2;

    h1 = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
    h2 = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];

    /// Periodicidade
//    pDir = (i==M.Nx) ? P[0][j]: P[i][j];
//    pEsq = (i==0) ? P[M.Nx-1][j] : P[i-1][j];

    /// Eixo de simetria
//    pEsq = (i==0) ? P[i][j] : P[i-1][j];

    /// FULL normalmente
    pDir = P[i][j];
    pEsq = P[i-1][j];

    pCentro = Interpolacao(h1, h2, pEsq, pDir);

    if( (i==0) || (i==M.Nx) )
        DeuProblema("GradPx: problema no indice (%d %d)\n", i, j);

    return (-pEsq*h2/(h1*(h1+h2))) + (pCentro*(h2-h1)/(h1*h2)) + (pDir*h1/(h2*(h1+h2)));
}

double DivT_u(double **Txx, double **Txy, double **Ttt, int i, int j, MALHA M)
{
    double valor = 0.0, h1, h2;
    static double txxDir = 0.0, txxEsq = 0.0, txxCentro = 0.0;
    static double txyDirCima = 0.0, txyEsqCima = 0.0, txyDirCentro = 0.0, txyEsqCentro = 0.0, txyDirBaixo = 0.0, txyEsqBaixo = 0.0;
    double rImenos1, rI, rImenosMeio;
    int iV;

    iV = (i==M.Nx) ? i-1 : i;

    if( velocidadeNewtoniana )
        return 0.0;

    if( M.tipoCoord==AXI_CILINDRICO ) {
        rImenos1 = M.r[i-1];
        rI = M.r[i];
        rImenosMeio = M.x[i];
    }
    else
        rImenos1 = rI = rImenosMeio = 1.0;

    /// Fazendo o termo da derivada em x
    h1 = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
    h2 = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];


    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
            txxDir = rI*Txx[0][j];
        else
            DeuProblema("DivTu: Problema condicao direita.\n\n");
    }
    else
        txxDir = rI*Txx[i][j];

    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
            txxEsq = rImenos1*Txx[M.Nx-1][j];
        else
            DeuProblema("DivTu: Problema condicao esquerda.\n\n");
    }
    else
        txxEsq = rImenos1*Txx[i-1][j];

    txxCentro = rImenosMeio*Interpolacao(h1, h2, txxEsq, txxDir);

    valor = ( 1.0/rImenosMeio )*( (-h2/(h1*(h1+h2)))*txxEsq + ((h2-h1)/(h1*h2))*txxCentro + (h1/(h2*(h1+h2)))*txxDir );


    /// Fazendo o termo da derivada em y
    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
            txyDirCentro = Txy[0][j];
        else
            DeuProblema("DivTu: Problema condicao esquerda.\n\n");
    }
    else
        txyDirCentro = Txy[i][j];

    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
            txyEsqCentro = Txy[M.Nx-1][j];
        else
            DeuProblema("DivTu: Problema condicao esquerda.\n\n");
    }
    else
        txyEsqCentro = Txy[i-1][j];

    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].tipoBoundary!=PERIODICIDADE ) {
        DeuProblema("Deu problema DivT_u\n");
    }
    else if( BAIXO_U(iV, j) ) {

        if( (i!=M.Nx) && M.pontosP[i][j+1].tipo!=EMPTY )
            txyDirCima = Txy[i][j+1];
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
            if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
                txyDirCima = Txy[0][j+1];
            else
                DeuProblema("DivTu: problema TxyEsqCima\n\n");
        }
        else if( M.pontosP[i-1][j].tipo!=EMPTY ) { //Joga pra baixo e esquerda
            txyDirCima = Txy[i-1][j];
//            DesenhaMalhaVTK(M, 0);
//            ImprimeInterfaceVTK(M, 10000);
//            DeuProblema("DIVT_u 1 %d %d\n", i, j);
        }
        else
            DeuProblema("1: Problema Tensor DivT_u %d %d\n", i, j);

        if( (i!=0) && M.pontosP[i-1][j+1].tipo!=EMPTY )
            txyEsqCima = Txy[i-1][j+1];
        else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
            if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
                txyEsqCima = Txy[M.Nx-1][j+1];
            else
                DeuProblema("DivTu: problema TxyEsqCima\n\n");
        }
        else if( M.pontosP[i-1][j].tipo!=EMPTY && M.pontosP[i][j+1].tipo!=EMPTY ) { //Extrapolando pra direita e pra baixo
            txyEsqCima = 0.5*(Txy[i-1][j] + Txy[i][j+1]);
//            DesenhaMalhaVTK(M, 0);
//            ImprimeInterfaceVTK(M, 100000);
//            DeuProblema("DIVT_u 2 %d %d\n", i, j);
        }
        else if( M.pontosP[i-1][j].tipo!=EMPTY ) {//Extrapola apenas pra baixo
            txyEsqCima = Txy[i-1][j];
//            DeuProblema("DIVT_u 3 %d %d\n", i, j);
        }
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("2: Problema Tensor DivT_u %d %d\n", i, j);
        }

        double dxEsq = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        double dxDir = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];

        valor += ( Interpolacao(dxEsq, dxDir, txyEsqCima, txyDirCima) -
                 Interpolacao(dxEsq, dxDir, txyEsqCentro, txyDirCentro) )/( 0.5*(M.dy[j] + M.dy[j+1]) );
    }
    else if( CIMA_U(iV, j) ) {

        if( M.pontosP[i][j-1].tipo!=EMPTY )
            txyDirBaixo = Txy[i][j-1];
        else if( M.pontosP[i][j].tipo!=EMPTY ) {//Extrapolando pra cima
            txyDirBaixo = Txy[i][j];
//            DesenhaMalhaVTK(M, 0);
//            ImprimeInterfaceVTK(M, 100000000);
//            DeuProblema("DIVT_u 4 %d %d\n", i, j);
        }
        else {
//            ImprimeInterfaceVTK(M, 100000);
            DesenhaMalhaVTK(M, 0);
            DeuProblema("3: Problema Tensor DivT_u %d %d\n", i, j);
        }

        if( (i!=0) && M.pontosP[i-1][j-1].tipo!=EMPTY )
            txyEsqBaixo = Txy[i-1][j-1];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
            txyEsqBaixo = Txy[M.Nx-1][j-1];
        else if( M.pontosP[i-1][j].tipo!=EMPTY ) {//Extrapolando pra cima
            txyEsqBaixo = Txy[i-1][j];
//            DesenhaMalhaVTK(M, 0);
//            ImprimeInterfaceVTK(M, 10000000);
//            DeuProblema("DIVT_u 5 %d %d\n", i, j);
        }
        else
            DeuProblema("4: Problema Tensor DivT_u\n");

        double dxEsq = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        double dxDir = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];

        valor += ( Interpolacao(dxEsq, dxDir, txyEsqCentro, txyDirCentro) -
                 Interpolacao(dxEsq, dxDir, txyEsqBaixo, txyDirBaixo) )/( 0.5*(M.dy[j-1] + M.dy[j]) );

    }
    else {


        if( (i!=M.Nx) && M.pontosP[i][j-1].tipo!=EMPTY )
            txyDirBaixo = Txy[i][j-1];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j].extremo==DIREITA )
            txyDirBaixo = Txy[0][j-1];
        else if( M.pontosP[i][j].tipo!=EMPTY && M.pontosP[i-1][j-1].tipo!=EMPTY ) //Extrapola pra cima e esquerda
            txyDirBaixo = 0.5*(Txy[i][j] + Txy[i-1][j-1]);
        else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapola pra cima
            txyDirBaixo = Txy[i][j];
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("5: Problema Tensor DivT_u %d %d\n", i, j);
        }



        if( (i!=0) && M.pontosP[i-1][j-1].tipo!=EMPTY )
            txyEsqBaixo = Txy[i-1][j-1];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j].extremo==ESQUERDA )
            txyEsqBaixo = Txy[M.Nx-1][j-1];
        else if( M.pontosP[i-1][j].tipo!=EMPTY ) //Extrapola pra cima
            txyEsqBaixo = Txy[i-1][j];
        else
            DeuProblema("6: Problema Tensor DivT_u %d %d\n", i, j);


        if( (i!=M.Nx) && M.pontosP[i][j+1].tipo!=EMPTY )
            txyDirCima = Txy[i][j+1];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j].extremo==DIREITA )
            txyDirCima = Txy[0][j+1];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j].extremo==ESQUERDA )
            txyDirCima = Txy[M.Nx-1][j+1];
        else if( (i!=0) && M.pontosP[i][j].tipo!=EMPTY && M.pontosP[i-1][j+1].tipo!=EMPTY ) //Extrapola pra esquerda e pra baixo
            txyDirCima = 0.5*(Txy[i][j] + Txy[i-1][j+1]);
        else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapola pra baixo
            txyDirCima = Txy[i][j];
        else
            DeuProblema("7: Problema Tensor DivT_u %d %d\n", i, j);



        if( (i!=0) && M.pontosP[i-1][j+1].tipo!=EMPTY )
            txyEsqCima = Txy[i-1][j+1];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j].extremo==ESQUERDA )
            txyEsqCima = Txy[M.Nx-1][j-1];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE && M.pontosU[i][j].extremo==DIREITA )
            txyEsqCima = Txy[0][j-1];
        else if( M.pontosP[i-1][j].tipo!=EMPTY ) //Extrapola pra baixo
            txyEsqCima = Txy[i-1][j];
        else
            DeuProblema("8: Problema Tensor DivT_u %d %d\n", i, j);

        h1 = 0.5*( M.dy[j-1] + M.dy[j] );
        h2 = 0.5*( M.dy[j] + M.dy[j+1] );

        double dxEsq = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        double dxDir = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];

        valor += (-h2/(h1*(h1+h2)))*Interpolacao(dxEsq, dxDir, txyEsqBaixo, txyDirBaixo) +
                 ((h2-h1)/(h1*h2))*Interpolacao(dxEsq, dxDir, txyEsqCentro, txyDirCentro) +
                 (h1/(h2*(h1+h2)))*Interpolacao(dxEsq, dxDir, txyEsqCima, txyDirCima);
    }



    /// Fazendo o termo sem derivada nenhuma (apenas no axissimetrico)
    if( M.tipoCoord == AXI_CILINDRICO )
        valor += - (1.0/rImenosMeio)*Interpolacao(0.5*M.dx[i-1], 0.5*M.dx[i], Ttt[i-1][j], Ttt[i][j]);

    return valor;
}

double DivT_v(double **Txy, double **Tyy, int i, int j, MALHA M)
{
    double valor = 0.0, h1, h2;
    static double tyyBaixo = 0.0, tyyCima = 0.0, tyyCentro = 0.0;
    static double txyEsqBaixo = 0.0, txyEsqCima = 0.0, txyCentroBaixo = 0.0, txyCentroCima = 0.0, txyDirBaixo = 0.0, txyDirCima = 0.0;

    if( velocidadeNewtoniana )
        return 0.0;

    static double rImenos1, rI, rImais1;
    if( M.tipoCoord == AXI_CILINDRICO ) {
        rI = M.r[i];
        rImenos1 = (i==0) ? M.r[i] - M.dx[i] : M.r[i-1];
        rImais1 = (i==M.Nx-1) ? M.r[i] + M.dx[i] : M.r[i+1];
    }
    else
        rI = rImenos1 = rImais1 = 1.0;


    h1 = 0.5*M.dy[j-1];
    h2 = 0.5*M.dy[j];
    tyyBaixo = Tyy[i][j-1];
    tyyCima = Tyy[i][j];
    tyyCentro = Interpolacao(h1, h2, tyyBaixo, tyyCima);

    valor = (-h2/(h1*(h1+h2)))*tyyBaixo + ((h2-h1)/(h1*h2))*tyyCentro + (h1/(h2*(h1+h2)))*tyyCima;

    if( M.pontosP[i][j-1].tipo==EMPTY )
        DeuProblema("Problema DivTv centroBaixo %d %d\n", i, j);
    else
        txyCentroBaixo = Txy[i][j-1];

    txyCentroCima = Txy[i][j];

    if( M.pontosV[i][j].tipo==BOUNDARY ) {
        DeuProblema("Deu problema DivT_v");
    }
    else if( M.pontosU[i][j].tipo==BOUNDARY || M.pontosU[i][j-1].tipo==BOUNDARY ) {

        if( M.pontosP[i+1][j-1].tipo!=EMPTY )
            txyDirBaixo = Txy[i+1][j-1];
        else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Extrapolando pra esquerda
            txyDirBaixo = Txy[i][j-1];
        else
            DeuProblema("1: Problema Tensor DivT_v %d %d\n", i, j);



        if( M.pontosP[i+1][j].tipo!=EMPTY )
            txyDirCima = Txy[i+1][j];
        else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapolando pra esquerda
            txyDirCima = Txy[i][j];
        else
            DeuProblema("2: Problema Tensor DivT_v %d %d\n", i, j);

        valor += (1.0/rI)*( rImais1*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], txyDirBaixo, txyDirCima) -
                            rI*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], txyCentroBaixo, txyCentroCima) )/( 0.5*(M.dx[i] + M.dx[i+1]) );
    }
    else if( M.pontosU[i+1][j].tipo==BOUNDARY || M.pontosU[i+1][j-1].tipo==BOUNDARY ) {

        if( M.pontosP[i-1][j-1].tipo!=EMPTY )
            txyEsqBaixo = Txy[i-1][j-1];
        else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Joga pra direita
            txyEsqBaixo = Txy[i][j-1];
        else
            DeuProblema("3: Problema Tensor DivT_v %d %d\n", i, j);

        if( M.pontosP[i-1][j].tipo!=EMPTY )
            txyEsqCima = Txy[i-1][j];
        else if( M.pontosP[i][j].tipo!=EMPTY )
            txyEsqCima = Txy[i][j];
        else
            DeuProblema("4: Problema Tensor DivT_v %d %d\n", i, j);

        valor += (1.0/rI)*( rI*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], txyCentroBaixo, txyCentroCima) -
                            rImenos1*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], txyEsqBaixo, txyEsqCima) )/( 0.5*(M.dx[i-1] + M.dx[i]) );
    }
    else {


        if( M.pontosP[i-1][j-1].tipo!=EMPTY )
            txyEsqBaixo = Txy[i-1][j-1];
        else if( M.pontosP[i][j-1].tipo!=EMPTY && M.pontosP[i-1][j].tipo!=EMPTY ) //Extrapola pra direita e cima
            txyEsqBaixo = 0.5*(Txy[i][j-1] + Txy[i-1][j]);
        else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Extrapola pra direita
            txyEsqBaixo = Txy[i][j-1];
        else
            DeuProblema("5: Problema Tensor DivT_v %d %d\n", i, j);



        if( M.pontosP[i-1][j].tipo!=EMPTY )
            txyEsqCima = Txy[i-1][j];
        else if( M.pontosP[i][j].tipo!=EMPTY && M.pontosP[i-1][j-1].tipo!=EMPTY ) //Extrapola pra direita e baixo
            txyEsqCima = 0.5*(Txy[i][j] + Txy[i-1][j-1]);
        else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapola pra direita
            txyEsqCima = Txy[i][j];
        else
            DeuProblema("6: Problema Tensor DivT_v %d %d\n", i, j);




        if( M.pontosP[i+1][j-1].tipo!=EMPTY )
            txyDirBaixo = Txy[i+1][j-1];
        else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Exrapola pra esquerda
            txyDirBaixo = Txy[i][j-1];
        else
            DeuProblema("7: Problema Tensor DivT_v %d %d\n", i, j);


        if( M.pontosP[i+1][j].tipo!=EMPTY )
            txyDirCima = Txy[i+1][j];
        else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapola pra esquerda
            txyDirCima = Txy[i][j];
        else
            DeuProblema("8: Problema Tensor DivT_v %d %d\n", i, j);

        h1 = 0.5*( M.dx[i-1] + M.dx[i] );
        h2 = 0.5*( M.dx[i] + M.dx[i+1] );

        valor += (1.0/rI)*( (-h2/(h1*(h1+h2)))*rImenos1*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], txyEsqBaixo, txyEsqCima) +
                 ((h2-h1)/(h1*h2))*rI*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], txyCentroBaixo, txyCentroCima) +
                 (h1/(h2*(h1+h2)))*rImais1*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], txyDirBaixo, txyDirCima) );
    }

    return valor;
}

double exact_velocity_channel(MALHA M, double y)
{
    double P_grad = 100.0;

    if( y<M.Bi/(M.Re*P_grad) )
		return 0.5*M.Re*P_grad*(1.0 - y*y) - M.Bi*(y + 1.0);
	else if( y<M.Bi/(M.Re*P_grad) )
		return 0.5*M.Re*P_grad - M.Bi + 0.5*M.Bi*M.Bi/(M.Re*P_grad);
    // else
	return 0.5*M.Re*P_grad*(1.0 - y*y) + M.Bi*(y - 1.0);
}

double exact_viscosity_channel(MALHA M, double y, int print)
{
    double dy = 0.0001;
    double u_top = exact_velocity_channel(M, y + dy);
    double u_bottom = exact_velocity_channel(M, y - dy);
    double dudy = (u_top - u_bottom)/(2*dy);
    double strain_rate = fabs(dudy)/sqrt(2.0);

    double visc = 1.0 + M.Bi*(1.0 - exp(-strain_rate/1e-04))/(2.0*strain_rate + 1e-20);
    return visc;
}

// The equivalent of the "divergent" part of the stress tensor term for generalized newtonian models
double DivT_u_NewtGen(MALHA M, double **U, double **V, double **ViscosityFunction, int i, int j)
{
    double uBaixo = 0.0, uCima = 0.0, uCentro = 0.0, uDir = 0.0, uEsq = 0.0;
    double vBaixoDir = 0.0, vBaixoEsq = 0.0, vCimaEsq = 0.0, vCimaDir = 0.0;
    double t = (M.passoTemporal+1)*M.dt;

    if( (i==0) || (M.pontosU[i-1][j].tipo==EMPTY) )
        DeuProblema("DivT_u_NewtGen: problema no indice i...\n");
    if( (i==M.Nx) || (M.pontosU[i+1][j].tipo==EMPTY) )
        DeuProblema("DivT_u_NewtGen: problema no indice i...\n");

    int iV = (i==M.Nx) ? i-1 : i;

    //Verificando o U(i, j-1)
    if( BAIXO_U(iV, j) ) {
        if( M.pontosV[iV][j].tipoBoundary == NOSLIP )
            uBaixo = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j].tipoBoundary == SLIP )
            uBaixo = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j].tipoBoundary == INFLOW )
            uBaixo = - U[i][j];
        else if( M.pontosV[iV][j].tipoBoundary == NEUMANN )
            uBaixo = U[i][j];
        else if( M.pontosV[iV][j].tipoBoundary == SIMETRIA )
            uBaixo = U[i][j];
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("DivT_u_NewtGen: problema baixo %d %d %d\n", i, j, M.pontosV[iV][j].tipoBoundary);
        }
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        uBaixo = U[i][j-1];

    //Encontrando o U(i, j+1)
    if( CIMA_U(i, j) ) {
        if( (M.pontosV[iV][j+1].tipoBoundary == NOSLIP) )
            uCima = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( (M.pontosV[iV][j+1].tipoBoundary == SLIP) )
            uCima = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], t, U );
        else if( M.pontosV[iV][j+1].tipoBoundary == INFLOW )
            uCima = - U[i][j];
        else if( M.pontosV[iV][j+1].tipoBoundary == NEUMANN )
            uCima = U[i][j];
        else if( M.pontosV[iV][j+1].tipoBoundary == SIMETRIA )
            uCima = U[i][j];
        else if( M.pontosU[i][j+1].tipoBoundary==PERIODICIDADE )
            uCima = 0.0; //period mudar isso?
        else
            DeuProblema("DivT_u_NewtGen: problema cima %d %d\n", i, j);
    }
    else
        uCima = U[i][j+1];

    //Encontrando o U(i-1, j)
    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE ) {
            uEsq = U[M.Nx-1][j];
            vBaixoEsq = V[M.Nx-1][j];
            vCimaEsq = V[M.Nx-1][j+1];
        }
        else
            DeuProblema("Cx: problema esquerda\n\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else {
        uEsq = U[i-1][j];
        vBaixoEsq = V[i-1][j];
        vCimaEsq = V[i-1][j+1];
    }


    //Encontrando o U(i+1, j)
    if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
        if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE ) {
            uDir = U[1][j];
            vBaixoDir = V[0][j];
            vCimaDir = V[0][j+1];
        }
        else
            DeuProblema("Cx: problema direita\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else {
        uDir = U[i+1][j];
        vBaixoDir = V[i][j];
        vCimaDir = V[i][j+1];
    }

    /// Derivatives: dudx and dudy
    double hx1 = (i==0) ? M.dx[i] : M.dx[i-1];
    double hx2 = (i==M.Nx) ? M.dx[i-1] : M.dx[i];
    double hy1 = (j==0) ? M.dy[j] : 0.5*(M.dy[j-1] + M.dy[j]);
    double hy2 = (j==M.Ny-1) ? M.dy[j] : 0.5*(M.dy[j] + M.dy[j+1]);
    double dudx = -uEsq*hx2/(hx1*(hx1+hx2)) + uCentro*(hx2-hx1)/(hx1*hx2) + uDir*hx1/(hx2*(hx1+hx2));
    double dudy = -uBaixo*hy2/(hy1*(hy1+hy2)) + uCentro*(hy2-hy1)/(hy1*hy2) + uCima*hy1/(hy2*(hy1+hy2));
    
    /// Derivative: dvdx
    hx1 = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
    hx2 = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];
    double vDir = 0.5*(vBaixoDir + vCimaDir);
    double vEsq = 0.5*(vBaixoEsq + vCimaEsq);
    double vCentro = Interpolacao(hx1, hx2, vEsq, vDir);
    double dvdx = -vEsq*hx2/(hx1*(hx1+hx2)) + vCentro*(hx2-hx1)/(hx1*hx2) + vDir*hx1/(hx2*(hx1+hx2));
    

    /// Derivative: dEta_dx
    double h1 = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
    double h2 = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];
    double etaDir = ViscosityFunction[i][j];
    double etaEsq = ViscosityFunction[i-1][j];
    double etaCentro = Interpolacao(h1, h2, etaEsq, etaDir);
    double dEta_dx = (-h2/(h1*(h1+h2)))*etaEsq + ((h2-h1)/(h1*h2))*etaCentro + (h1/(h2*(h1+h2)))*etaDir;


    /// Derivative: dEta_dy
    /// This one is complicated, because we change the derivative direction if we are near a top/bottom boundary
    /// Near top boundary: i do backwards finite differences (going downwards to avoid the wall at the top)
    /// Near bottom boundary: i do forward finite differences (going upwards to avoid the wall at the bottom)
    /// Near everywhere else: centered finite differences, no problem!
    if( M.celulas[i-1][j]==EMPTY || M.celulas[i][j]==EMPTY ) 
        DeuProblema("CELULA INVALIDA\n");
    
    double dEta_dy;
    double etaDirCentro = ViscosityFunction[i][j];
    double etaEsqCentro = ViscosityFunction[i-1][j];
    double etaEsqBaixo = 0.0, etaDirBaixo = 0.0, etaEsqCima = 0.0, etaDirCima = 0.0;

    /// Case 1: near a boundary at the bottom
    if( BAIXO_U(iV, j) ) {

        if( (i!=M.Nx) && M.pontosP[i][j+1].tipo!=EMPTY )
            etaDirCima = ViscosityFunction[i][j+1];
        // else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==DIREITA ) {
        //     if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
        //         etaDirCima = ViscosityFunction[0][j+1];
        //     else
        //         DeuProblema("DivTu: problema TxyEsqCima\n\n");
        // }
        // else if( M.pontosP[i-1][j].tipo!=EMPTY ) { //Joga pra baixo e esquerda
        //     etaDirCima = ViscosityFunction[i-1][j];
        // }
        else
            DeuProblema("1: Problema Tensor DivT_u %d %d\n", i, j);

        if( (i!=0) && M.pontosP[i-1][j+1].tipo!=EMPTY )
            etaEsqCima = ViscosityFunction[i-1][j+1];
        // else if( M.pontosU[i][j].tipo==BOUNDARY && M.pontosU[i][j].extremo==ESQUERDA ) {
        //     if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
        //         etaEsqCima = ViscosityFunction[M.Nx-1][j+1];
        //     else
        //         DeuProblema("DivTu: problema TxyEsqCima\n\n");
        // }
        // else if( M.pontosP[i-1][j].tipo!=EMPTY && M.pontosP[i][j+1].tipo!=EMPTY ) { //Extrapolando pra direita e pra baixo
        //     etaEsqCima = 0.5*(ViscosityFunction[i-1][j] + ViscosityFunction[i][j+1]);
        // }
        // else if( M.pontosP[i-1][j].tipo!=EMPTY ) {//Extrapola apenas pra baixo
        //     etaEsqCima = ViscosityFunction[i-1][j];
        // }
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("2: Problema Tensor DivT_u %d %d\n", i, j);
        }

        double dxEsq = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        double dxDir = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];

        dEta_dy = ( Interpolacao(dxEsq, dxDir, etaEsqCima, etaDirCima) -
                 Interpolacao(dxEsq, dxDir, etaEsqCentro, etaDirCentro) )/( 0.5*(M.dy[j] + M.dy[j+1]) );
    }
    else if( CIMA_U(iV, j) ) {

        if( M.pontosP[i][j-1].tipo!=EMPTY )
            etaDirBaixo = ViscosityFunction[i][j-1];
//         else if( M.pontosP[i][j].tipo!=EMPTY ) {//Extrapolando pra cima
//             etaDirBaixo = ViscosityFunction[i][j];
// //            DesenhaMalhaVTK(M, 0);
// //            ImprimeInterfaceVTK(M, 100000000);
// //            DeuProblema("DIVT_u 4 %d %d\n", i, j);
//         }
        else {
//            ImprimeInterfaceVTK(M, 100000);
            DesenhaMalhaVTK(M, 0);
            DeuProblema("3: Problema Tensor DivT_u %d %d\n", i, j);
        }

        if( (i!=0) && M.pontosP[i-1][j-1].tipo!=EMPTY )
            etaEsqBaixo = ViscosityFunction[i-1][j-1];
        // else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
        //     etaEsqBaixo = ViscosityFunction[M.Nx-1][j-1];
        // else if( M.pontosP[i-1][j].tipo!=EMPTY ) {//Extrapolando pra cima
        //     etaEsqBaixo = ViscosityFunction[i-1][j];
//            DesenhaMalhaVTK(M, 0);
//            ImprimeInterfaceVTK(M, 10000000);
//            DeuProblema("DIVT_u 5 %d %d\n", i, j);
        // }
        else
            DeuProblema("4: Problema Tensor DivT_u\n");

        double dxEsq = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        double dxDir = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];

        dEta_dy = ( Interpolacao(dxEsq, dxDir, etaEsqCentro, etaDirCentro) -
                 Interpolacao(dxEsq, dxDir, etaEsqBaixo, etaDirBaixo) )/( 0.5*(M.dy[j-1] + M.dy[j]) );
    }
    else {


        if( (i!=M.Nx) && M.pontosP[i][j-1].tipo!=EMPTY )
            etaDirBaixo = ViscosityFunction[i][j-1];
        // else if( M.pontosP[i][j].tipo!=EMPTY && M.pontosP[i-1][j-1].tipo!=EMPTY ) //Extrapola pra cima e esquerda
        //     etaDirBaixo = 0.5*(ViscosityFunction[i][j] + ViscosityFunction[i-1][j-1]);
        // else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapola pra cima
        //     etaDirBaixo = ViscosityFunction[i][j];
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("5: Problema Tensor DivT_u %d %d\n", i, j);
        }



        if( (i!=0) && M.pontosP[i-1][j-1].tipo!=EMPTY )
            etaEsqBaixo = ViscosityFunction[i-1][j-1];
        // else if( M.pontosP[i-1][j].tipo!=EMPTY ) //Extrapola pra cima
        //     etaEsqBaixo = ViscosityFunction[i-1][j];
        else
            DeuProblema("6: Problema Tensor DivT_u %d %d\n", i, j);


        if( (i!=M.Nx) && M.pontosP[i][j+1].tipo!=EMPTY )
            etaDirCima = ViscosityFunction[i][j+1];
        // else if( (i!=0) && M.pontosP[i][j].tipo!=EMPTY && M.pontosP[i-1][j+1].tipo!=EMPTY ) { //Extrapola pra esquerda e pra baixo
        //     etaDirCima = 0.5*(ViscosityFunction[i][j] + ViscosityFunction[i-1][j+1]);
        //     DeuProblema("test 1\n");
        // }
        // else if( M.pontosP[i][j].tipo!=EMPTY ) { //Extrapola pra baixo
        //     etaDirCima = ViscosityFunction[i][j];
        //     DeuProblema("test 2\n");
        // }
        else
            DeuProblema("7: Problema Tensor DivT_u %d %d\n", i, j);



        if( (i!=0) && M.pontosP[i-1][j+1].tipo!=EMPTY )
            etaEsqCima = ViscosityFunction[i-1][j+1];
        // else if( M.pontosP[i-1][j].tipo!=EMPTY ) //Extrapola pra baixo
        //     etaEsqCima = ViscosityFunction[i-1][j];
        else
            DeuProblema("8: Problema Tensor DivT_u %d %d\n", i, j);

       

        h1 = 0.5*( M.dy[j-1] + M.dy[j] );
        h2 = 0.5*( M.dy[j] + M.dy[j+1] );

        double dxEsq = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        double dxDir = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];


        dEta_dy = (-h2/(h1*(h1+h2)))*Interpolacao(dxEsq, dxDir, etaEsqBaixo, etaDirBaixo) +
                 ((h2-h1)/(h1*h2))*Interpolacao(dxEsq, dxDir, etaEsqCentro, etaDirCentro) +
                 (h1/(h2*(h1+h2)))*Interpolacao(dxEsq, dxDir, etaEsqCima, etaDirCima);

        // if( M.celulas[i-1][j]==SURFACE )
        //     PrintDebug("LEFT %d %d: %lf.....\n", i, j, dEta_dy);
        // if( M.celulas[i][j]==SURFACE )
        //     PrintDebug("RIGHT %d %d: %lf.....\n", i, j, dEta_dy);
    }

    // int print = ( M.celulas[i][j]==SURFACE ) ? 1 : 0;
    // // Just a quick test using the exact solution...
    // double y = 0.5*( M.y[j] + M.y[j+1] );
    // double dy = 0.0001;
    // double eta_top = exact_viscosity_channel(M, y + dy, print*2);
    // double eta_bottom = exact_viscosity_channel(M, y - dy, print*3);
    // double exact_dEta_dy = (eta_top - eta_bottom)/(2*dy);


    return 2.0*dEta_dx*dudx + dEta_dy*( dudy + dvdx );
}



double DivT_v_NewtGen(MALHA M, double **U, double **V, double **ViscosityFunction, int i, int j)
{
    double vEsq = 0.0, vDir = 0.0, vBaixo = 0.0, vCima = 0.0, vCentro;
    double uBaixoEsq = 0.0, uBaixoDir = 0.0, uCimaEsq = 0.0, uCimaDir = 0.0;
    
    int jU = (j==M.Ny) ? j-1 : j;
    
    // Encontrando o V(i-1, j)
    if( ESQUERDA_V(i, j) ) {
        if( M.pontosU[i][jU].tipoBoundary == NOSLIP )
            vEsq = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i][jU].tipoBoundary == SLIP )
            vEsq = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i][jU].tipoBoundary == INFLOW )
            vEsq = - V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == NEUMANN )
            vEsq = V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == SIMETRIA )
            vEsq = V[i][j];
        else if( M.pontosU[i][jU].tipoBoundary == PERIODICIDADE )
            vEsq = V[M.Nx-1][j];
        else if( M.pontosU[i][jU].tipoBoundary == EIXO_AXISSIMETRIA )
            vEsq = V[i][j];
        else
            DeuProblema("Derivadas Cy: Problema esquerda\n\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        vEsq = V[i-1][j];

    //Encontrando o V(i+1, j)
    if( DIREITA_V(i, j) ) {
        if( M.pontosU[i+1][jU].tipoBoundary==NEUMANN )
            vDir = V[i][j];
        else if( M.pontosU[i+1][jU].tipoBoundary==SIMETRIA )
            vDir = V[i][j];
        else if( M.pontosU[i+1][jU].tipoBoundary==NOSLIP )
            vDir = - V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][jU].tipoBoundary==SLIP )
            vDir = V[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, V );
        else if( M.pontosU[i+1][jU].tipoBoundary==INFLOW )
            vDir = - V[i][j];
        else if( M.pontosU[i+1][jU].tipoBoundary==PERIODICIDADE )
            vDir = V[0][j];
        else
            DeuProblema("Derivadas Cy: Problema direita\n\n");
        //Fazer os outros casos de tipoBoundary (Neumann, noslip, inflow)
    }
    else
        vDir = V[i+1][j];

    //Encontrando o V(i, j-1)
    if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].extremo==BAIXO )
        DeuProblema("Cy: Nao era pra ter entrado aqui %d %d!!! \n\n", i, j);
    else {
        vBaixo = V[i][j-1]; //Neumann ordem 2
        uBaixoEsq = U[i][j-1];
        uBaixoDir = U[i+1][j-1];
    }

    //Encontrando o V(i, j+1)
    if( M.pontosV[i][j].tipo==BOUNDARY && M.pontosV[i][j].extremo==CIMA )
        DeuProblema("Cy: Nao era pra ter entrado aqui %d %d!!! \n\n", i, j);
    else {
        vCima = V[i][j+1]; //Neumann ordem 2
        uCimaEsq = U[i][j];
        uCimaDir = U[i+1][j];
    }

    vCentro = V[i][j];

    /// Derivatives: dvdx and dvdy
    double hy1 = (j==0) ? M.dy[j] : M.dy[j-1];
    double hy2 = (j==M.Ny) ? M.dy[j-1] : M.dy[j];
    double hx1 = (i==0) ? M.dx[i] : 0.5*(M.dx[i-1] + M.dx[i]);
    double hx2 = (i==M.Nx-1) ? M.dx[i] : 0.5*(M.dx[i] + M.dx[i+1]);
    double dvdx = -vEsq*hx2/(hx1*(hx1+hx2)) + vCentro*(hx2-hx1)/(hx1*hx2) + vDir*hx1/(hx2*(hx1+hx2));
    double dvdy = -vBaixo*hy2/(hy1*(hy1+hy2)) + vCentro*(hy2-hy1)/(hy1*hy2) + vCima*hy1/(hy2*(hy1+hy2));
    
    /// Derivative: dudy
    hy1 = (j==0) ? 0.5*M.dy[j] : 0.5*M.dy[j-1];
    hy2 = (j==M.Ny) ? 0.5*M.dy[j-1] : 0.5*M.dy[j];
    double uBaixo = 0.5*(uBaixoDir + uBaixoEsq);
    double uCima = 0.5*(uCimaDir + uCimaEsq);
    double uCentro = Interpolacao(hy1, hy2, uBaixo, uCima);
    double dudy = -uBaixo*hy2/(hy1*(hy1+hy2)) + uCentro*(hy2-hy1)/(hy1*hy2) + uCima*hy1/(hy2*(hy1+hy2));
    

    /// Derivative: dEta_dy
    double h1 = (j==0) ? 0.5*M.dy[j] : 0.5*M.dy[j-1];
    double h2 = (j==M.Ny) ? 0.5*M.dy[j-1] : 0.5*M.dy[j];
    double etaCima = ViscosityFunction[i][j];
    double etaBaixo = ViscosityFunction[i][j-1];
    double etaCentro = Interpolacao(h1, h2, etaBaixo, etaCima);
    double dEta_dy = (-h2/(h1*(h1+h2)))*etaBaixo + ((h2-h1)/(h1*h2))*etaCentro + (h1/(h2*(h1+h2)))*etaCima;




    /// Derivative: dEta_dx
    /// This one is complicated, because we change the derivative direction if we are near a right/left boundary
    /// Near right boundary: i do backwards finite differences (going left  to avoid the wall at the right)
    /// Near left boundary: i do forward finite differences (going right to avoid the wall at the left)
    /// Near everywhere else: centered finite differences, no problem!
    double dEta_dx;
    double etaCentroBaixo = 0.0, etaCentroCima = 0.0, etaDirBaixo = 0.0, etaDirCima = 0.0, etaEsqBaixo = 0.0, etaEsqCima = 0.0;
    if( M.pontosP[i][j-1].tipo==EMPTY )
        DeuProblema("Problema DivTv centroBaixo %d %d\n", i, j);
    else
        etaCentroBaixo = ViscosityFunction[i][j-1];

    etaCentroCima = ViscosityFunction[i][j];

    if( M.pontosV[i][j].tipo==BOUNDARY ) {
        DeuProblema("Deu problema DivT_v");
        dEta_dx = 0.0;
    }
    else if( M.pontosU[i][j].tipo==BOUNDARY || M.pontosU[i][j-1].tipo==BOUNDARY ) {

        if( M.pontosP[i+1][j-1].tipo!=EMPTY )
            etaDirBaixo = ViscosityFunction[i+1][j-1];
        else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Extrapolando pra esquerda
            etaDirBaixo = ViscosityFunction[i][j-1];
        else
            DeuProblema("1: Problema Tensor DivT_v %d %d\n", i, j);



        if( M.pontosP[i+1][j].tipo!=EMPTY )
            etaDirCima = ViscosityFunction[i+1][j];
        else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapolando pra esquerda
            etaDirCima = ViscosityFunction[i][j];
        else
            DeuProblema("2: Problema Tensor DivT_v %d %d\n", i, j);

        dEta_dx = ( Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], etaDirBaixo, etaDirCima) -
                        Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], etaCentroBaixo, etaCentroCima) )/( 0.5*(M.dx[i] + M.dx[i+1]) );
    }
    else if( M.pontosU[i+1][j].tipo==BOUNDARY || M.pontosU[i+1][j-1].tipo==BOUNDARY ) {

        if( M.pontosP[i-1][j-1].tipo!=EMPTY )
            etaEsqBaixo = ViscosityFunction[i-1][j-1];
        else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Joga pra direita
            etaEsqBaixo = ViscosityFunction[i][j-1];
        else
            DeuProblema("3: Problema Tensor DivT_v %d %d\n", i, j);

        if( M.pontosP[i-1][j].tipo!=EMPTY )
            etaEsqCima = ViscosityFunction[i-1][j];
        else if( M.pontosP[i][j].tipo!=EMPTY )
            etaEsqCima = ViscosityFunction[i][j];
        else
            DeuProblema("4: Problema Tensor DivT_v %d %d\n", i, j);

        dEta_dx = ( Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], etaCentroBaixo, etaCentroCima) -
                            Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], etaEsqBaixo, etaEsqCima) )/( 0.5*(M.dx[i-1] + M.dx[i]) );
    }
    else {


        if( M.pontosP[i-1][j-1].tipo!=EMPTY )
            etaEsqBaixo = ViscosityFunction[i-1][j-1];
        else if( M.pontosP[i][j-1].tipo!=EMPTY && M.pontosP[i-1][j].tipo!=EMPTY ) //Extrapola pra direita e cima
            etaEsqBaixo = 0.5*(ViscosityFunction[i][j-1] + ViscosityFunction[i-1][j]);
        else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Extrapola pra direita
            etaEsqBaixo = ViscosityFunction[i][j-1];
        else
            DeuProblema("5: Problema Tensor DivT_v %d %d\n", i, j);



        if( M.pontosP[i-1][j].tipo!=EMPTY )
            etaEsqCima = ViscosityFunction[i-1][j];
        else if( M.pontosP[i][j].tipo!=EMPTY && M.pontosP[i-1][j-1].tipo!=EMPTY ) //Extrapola pra direita e baixo
            etaEsqCima = 0.5*(ViscosityFunction[i][j] + ViscosityFunction[i-1][j-1]);
        else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapola pra direita
            etaEsqCima = ViscosityFunction[i][j];
        else
            DeuProblema("6: Problema Tensor DivT_v %d %d\n", i, j);




        if( M.pontosP[i+1][j-1].tipo!=EMPTY )
            etaDirBaixo = ViscosityFunction[i+1][j-1];
        else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Exrapola pra esquerda
            etaDirBaixo = ViscosityFunction[i][j-1];
        else
            DeuProblema("7: Problema Tensor DivT_v %d %d\n", i, j);


        if( M.pontosP[i+1][j].tipo!=EMPTY )
            etaDirCima = ViscosityFunction[i+1][j];
        else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapola pra esquerda
            etaDirCima = ViscosityFunction[i][j];
        else
            DeuProblema("8: Problema Tensor DivT_v %d %d\n", i, j);

        h1 = 0.5*( M.dx[i-1] + M.dx[i] );
        h2 = 0.5*( M.dx[i] + M.dx[i+1] );

        dEta_dx = ( (-h2/(h1*(h1+h2)))*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], etaEsqBaixo, etaEsqCima) +
                 ((h2-h1)/(h1*h2))*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], etaCentroBaixo, etaCentroCima) +
                 (h1/(h2*(h1+h2)))*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], etaDirBaixo, etaDirCima) );
    }

    return 2.0*dEta_dy*dvdy + dEta_dx*( dudy + dvdx );
}

double DivT_t(double **Txt, double **Tyt, int i, int j, MALHA M)
{
    double valor = 0.0;
    static double tytBaixo, tytCentro, tytCima;
    static double txtEsq, txtCentro, txtDir;
    static double rImenos1, rI, rImais1;

//    return 0.0;

    rI = M.r[i];
    rImenos1 = (i==0) ? M.r[i] - M.dx[i] : M.r[i-1];
    rImais1 = (i==M.Nx-1) ? M.r[i] + M.dx[i] : M.r[i+1];

    if( (i==0) || (M.celulas[i-1][j]!=FULL) ) { //Nao tem pontos pra esquerda. Faco diferencas progressivas
        txtCentro = rI*rI*Txt[i][j];
        txtDir = rImais1*rImais1*Txt[i+1][j];

        valor += ( txtDir - txtCentro )/( 0.5*(rI*rI)*(M.dx[i] + M.dx[i+1]) );
    }
    else if( (i==M.Nx-1) || (M.celulas[i+1][j]!=FULL) ) { //Nao tem pontos pra direita. Faco diferencas regressivas
        txtCentro = rI*rI*Txt[i][j];
        txtEsq = rImenos1*rImenos1*Txt[i-1][j];

        valor += ( txtCentro - txtEsq )/( 0.5*(rI*rI)*(M.dx[i] + M.dx[i-1]) );
    }
    else {
        txtCentro = rI*rI*Txt[i][j];
        txtEsq = rImenos1*rImenos1*Txt[i-1][j];
        txtDir = rImais1*rImais1*Txt[i+1][j];

        double h1 = 0.5*( M.dx[i-1] + M.dx[i] );
        double h2 = 0.5*( M.dx[i] + M.dx[i+1] );

        valor += (1.0/(rI*rI))*( (-h2/(h1*(h1+h2)))*txtEsq + ((h2-h1)/(h1*h2))*txtCentro + (h1/(h2*(h1+h2)))*txtDir );
    }



    if( (j==0) || (M.celulas[i][j-1]!=FULL) ) { //Nao tem pontos pra baixo. Faco diferencas progressivas
        tytCentro = Tyt[i][j];
        tytCima = Tyt[i][j+1];

//        PrintDebug("(%2d, %2d): %20.15lf %20.15lf %20.15lf\n", i, j, tytCentro, tytCima, ( tytCima - tytCentro )/( 0.5*(M.dy[j] + M.dy[j+1]) ));

        valor += ( tytCima - tytCentro )/( 0.5*(M.dy[j] + M.dy[j+1]) );
    }
    else if( (j==M.Ny-1) || (M.celulas[i][j+1]!=FULL) ) { //Nao tem pontos pra cima. Faco diferencas regressivas
        tytCentro = Tyt[i][j];
        tytBaixo = Tyt[i][j-1];

//        PrintDebug("(%2d, %2d): %20.15lf %20.15lf %20.15lf\n", i, j, tytCentro, tytBaixo, ( tytCentro - tytBaixo )/( 0.5*(M.dy[j] + M.dy[j-1]) ));

        valor += ( tytCentro - tytBaixo )/( 0.5*(M.dy[j] + M.dy[j-1]) );
    }
    else {
        tytCentro = Tyt[i][j];
        tytBaixo = Tyt[i][j-1];
        tytCima = Tyt[i][j+1];

        double h1 = 0.5*( M.dy[j-1] + M.dy[j] );
        double h2 = 0.5*( M.dy[j] + M.dy[j+1] );

        valor += (-h2/(h1*(h1+h2)))*tytBaixo + ((h2-h1)/(h1*h2))*tytCentro + (h1/(h2*(h1+h2)))*tytCima;
    }

//    if( j==0 || j==M.Ny-1 )
//        PrintDebug("(%2d, %2d): %20.15lf %20.15lf %20.15lf\n", i, j, tytCentro, tytBaixo, valor);

    return valor;
}

double Conv_Tensor(double **T, double **U, double **V, int i, int j, const char *Tipo, MALHA M)
{
    static double fluxo1, fluxo2;
    static double vel, phiR, phiU, phiD, parte1, parte2;
    static double coordR, coordU, coordD, coordF;
    static double lixo = -1e10;
    static char boundaryEsq, boundaryDir, boundaryBaixo, boundaryCima;
    static double tCima = 0.0, tBaixo = 0.0, tEsq = 0.0, tDir = 0.0;


//    return 0.0;

    //Verifica se esta no contorno superior
    if( M.pontosV[i][j+1].tipo==BOUNDARY ) {
        if( M.pontosV[i][j+1].tipoBoundary == NOSLIP || M.pontosV[i][j+1].tipoBoundary == SLIP ){
//            tCima = - T[i][j] + 2.0*TensorParede(M, i, j+1, Tipo, 'v');
            tCima = 2.0*T[i][j] - T[i][j-1];
            if( (j!=M.Ny-1) && (i!=0) && M.pontosP[i-1][j+1].tipo==FULL ) {
                tCima = 0.5*tCima + 0.5*(2.0*T[i-1][j+1] - T[i-2][j+1]);
            }
            if( (j!=M.Ny-1) && (i!=M.Nx-1) && M.pontosP[i+1][j+1].tipo==FULL )
                tCima = 0.5*tCima + 0.5*(2.0*T[i+1][j+1] - T[i+2][j+1]);
        }
        else if( M.pontosV[i][j+1].tipoBoundary == INFLOW )
            tCima = T[i][j];
        else if( M.pontosV[i][j+1].tipoBoundary == NEUMANN )
            tCima = T[i][j];
        else if( M.pontosV[i][j+1].tipoBoundary == SIMETRIA )
            tCima = T[i][j];
        else DeuProblema("Problema ConvTensor Cima %d %d\n", i, j);
    }
    else if( M.pontosP[i][j+1].tipo==EMPTY ) {
        //DeuProblema("SOH ERA PRA ENTRAR AQUI COM SUP LIVRE\n\n");
        if( M.pontosP[i][j].tipo!=EMPTY )
            tCima = T[i][j];
        else
            DeuProblema("EMPTY: Problema ConvTensor Cima %d %d\n", i, j);
    }
    else {
        tCima = T[i][j+1];
    }

    //Verifica se esta no contorno infrior
    if( M.pontosV[i][j].tipo==BOUNDARY ) {
        if( M.pontosV[i][j].tipoBoundary == NOSLIP || M.pontosV[i][j].tipoBoundary==SLIP ) {
//            tBaixo = - T[i][j] + 2.0*TensorParede(M, i, j, Tipo, 'v');
            tBaixo = 2.0*T[i][j] - T[i][j+1];
            if( (j!=0) && (i!=0) && M.pontosP[i-1][j-1].tipo==FULL )
                tBaixo = 0.5*tBaixo + 0.5*(2.0*T[i-1][j-1] - T[i-2][j-1]);
            if( (j!=0) && (i!=M.Nx-1) && M.pontosP[i+1][j-1].tipo==FULL )
                tBaixo = 0.5*tBaixo + 0.5*(2.0*T[i+1][j-1] - T[i+2][j-1]);
        }
        else if( M.pontosV[i][j].tipoBoundary == INFLOW )
            //tBaixo = T[i][j];
            tBaixo = - T[i][j] + 2.0*ValorInflowTensor(i, j, V, Tipo, M, 'h');
        else if( M.pontosV[i][j].tipoBoundary == NEUMANN )
            tBaixo = T[i][j];
        else if( M.pontosV[i][j].tipoBoundary == SIMETRIA )
            tBaixo = T[i][j];
        else DeuProblema("Problema ConvTensor Baixo %d %d\n", i, j);
    }
    else if( M.pontosP[i][j-1].tipo==EMPTY ) {
        //DeuProblema("SOH ERA PRA ENTRAR AQUI COM SUP LIVRE\n\n");
        if( M.pontosP[i][j].tipo!=EMPTY )
            tBaixo = T[i][j];
        else
            DeuProblema("EMPTY: Problema ConvTensor Baixo %d %d\n", i, j);
    }
    else {
        tBaixo = T[i][j-1];
    }



    //Verifica se esta no contorno da esquerda
    if( M.pontosU[i][j].tipo==BOUNDARY ) {
        if(M.pontosU[i][j].tipoBoundary == NEUMANN)
            tEsq = T[i][j];
        else if( (M.pontosU[i][j].tipoBoundary == NOSLIP) || (M.pontosU[i][j].tipoBoundary == SLIP) ) {
//            tEsq = - T[i][j] + 2.0*TensorParede(M, i, j, Tipo, 'u');
            tEsq = 2.0*T[i][j] - T[i+1][j];
            if( (j!=0) && (i!=0) && M.pontosP[i-1][j-1].tipo==FULL )
                tEsq = 0.5*tEsq + 0.5*(2.0*T[i-1][j-1] - T[i-1][j-2]);
            if( (j!=M.Ny-1) && (i!=0) && M.pontosP[i-1][j+1].tipo==FULL )
                tEsq = 0.5*tEsq + 0.5*(2.0*T[i-1][j+1] - T[i-1][j+2]);
        }
        else if( M.pontosU[i][j].tipoBoundary == INFLOW )
            tEsq = - T[i][j] + 2.0*ValorInflowTensor(i, j, U, Tipo, M, 'v');
            //tEsq = T[i][j];
        else if( M.pontosU[i][j].tipoBoundary == EIXO_AXISSIMETRIA ) {
            tEsq = T[i][j];
            //DeuProblema("ConvTensor: CONFIRMAR A CONDICAO SIMETRIA PRA TER CERTEZA %d %d\n", i, j);
        }
        else DeuProblema("Problema Esquerda ConvTensor: %d %d\n", i, j);
    }
    else if( M.pontosP[i-1][j].tipo==EMPTY ) {
        //DeuProblema("SOH ERA PRA ENTRAR AQUI COM SUP LIVRE\n\n");
        if( M.pontosP[i][j].tipo!=EMPTY )
            tEsq = T[i][j];
        else
            DeuProblema("EMPTY: Problema ConvTensor Esquerda %d %d\n", i, j);
    }
    else {
        tEsq = T[i-1][j];
    }

    //Verifica se esta no contorno da direita
    if( M.pontosU[i+1][j].tipo==BOUNDARY ) {
        if( M.pontosU[i+1][j].tipoBoundary==NEUMANN )
            tDir = T[i][j];
        else if( (M.pontosU[i+1][j].tipoBoundary==NOSLIP) || (M.pontosU[i+1][j].tipoBoundary==SLIP) ) {
//            tDir = - T[i][j] + 2.0*TensorParede(M, i+1, j, Tipo, 'u');
            tDir = 2.0*T[i][j] - T[i-1][j];
            if( (j!=M.Ny-1) && (i!=M.Nx-1) && M.pontosP[i+1][j+1].tipo==FULL )
                tDir = 0.5*tDir + 0.5*(2.0*T[i+1][j+1] - T[i+1][j+2]);
            if( (j!=0) && (i!=M.Nx-1) && M.pontosP[i+1][j-1].tipo==FULL )
                tDir = 0.5*tDir + 0.5*(2.0*T[i+1][j-1] - T[i+1][j-2]);
        }
        else if( M.pontosU[i+1][j].tipoBoundary==INFLOW )
            tDir = T[i][j];
        else DeuProblema("Problema Direita ConvTensor: %d %d\n", i, j);
    }
    else if( M.pontosP[i+1][j].tipo==EMPTY ) {
        //DeuProblema("SOH ERA PRA ENTRAR AQUI COM SUP LIVRE\n\n");
        if( M.pontosP[i][j].tipo!=EMPTY )
            tDir = T[i][j];
        else
            DeuProblema("EMPTY: Problema ConvTensor Direita %d %d\n", i, j);
    }
    else {
        tDir = T[i+1][j];
    }



    boundaryDir = boundaryEsq = boundaryBaixo = boundaryCima = 0;
    if( (i<2) || M.pontosU[i][j].tipo==BOUNDARY || M.pontosU[i-1][j].tipo==BOUNDARY
       || M.pontosP[i-1][j].tipo==SURFACE || M.pontosP[i-2][j].tipo==SURFACE )
        boundaryEsq = 1;
    if( (i>M.Nx-3) || M.pontosU[i+1][j].tipo==BOUNDARY || M.pontosU[i+2][j].tipo==BOUNDARY
       || M.pontosP[i+1][j].tipo==SURFACE || M.pontosP[i+2][j].tipo==SURFACE )
        boundaryDir = 1;
    if( (j<2) || M.pontosV[i][j].tipo==BOUNDARY || M.pontosV[i][j-1].tipo==BOUNDARY
       || M.pontosP[i][j-1].tipo==SURFACE || M.pontosP[i][j-2].tipo==SURFACE )
        boundaryBaixo = 1;
    if( (j>M.Ny-3) || M.pontosV[i][j+1].tipo==BOUNDARY || M.pontosV[i][j+2].tipo==BOUNDARY
       || M.pontosP[i][j+1].tipo==SURFACE || M.pontosP[i][j+2].tipo==SURFACE )
        boundaryCima = 1;

    double rImaisMeio, rImenosMeio, rI;

    if( M.tipoCoord == AXI_CILINDRICO ) {
        rImaisMeio = M.r[i] + 0.5*M.dx[i];
        rImenosMeio = M.r[i] - 0.5*M.dx[i];
        rI = M.r[i];
    }
    else
        rImaisMeio = rImenosMeio = rI = 1.0;

    vel = rImaisMeio*U[i+1][j];
    coordF = M.x[i+1];
    if( vel>=0 ) {
        coordR = (i==0) ? M.x[i] - 0.5*M.dx[i] : M.x[i-1] + 0.5*M.dx[i-1];
        coordU = M.x[i] + 0.5*M.dx[i];
        coordD = (i==M.Nx-1) ? M.x[i+1] + 0.5*M.dx[i] : M.x[i+1] + 0.5*M.dx[i+1];

        phiR = tEsq;
        phiU = T[i][j];
        phiD = tDir;
    }
    else {
        coordR = boundaryDir ? lixo : M.x[i+2] + 0.5*M.dx[i+2];
        coordU = (i==M.Nx-1) ? M.x[i+1] + 0.5*M.dx[i] : M.x[i+1] + 0.5*M.dx[i+1];
        coordD = M.x[i] + 0.5*M.dx[i];


        phiR = boundaryDir ? lixo : T[i+2][j]; ///ATENCAO
        phiU = tDir;
        phiD = T[i][j];
    }
    parte1 = boundaryDir ? vel*FuncaoUpwind(phiU) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);
    //parte1 = vel*FuncaoUpwind(phiU);

    vel = rImenosMeio*U[i][j];
    coordF = M.x[i];
    if( vel>=0 ) {
        coordR = boundaryEsq ? lixo : M.x[i-2] + 0.5*M.dx[i-2];
        coordU = (i==0) ? M.x[i] - 0.5*M.dx[i] : M.x[i-1] + 0.5*M.dx[i-1];
        coordD = M.x[i] + 0.5*M.dx[i];

        phiR = boundaryEsq ? lixo : T[i-2][j]; ///ATENCAO
        phiU = tEsq;
        phiD = T[i][j];
    }
    else {
        coordR = (i==M.Nx-1) ? M.x[i+1] + 0.5*M.dx[i] : M.x[i+1] + 0.5*M.dx[i+1];
        coordU = M.x[i] + 0.5*M.dx[i];
        coordD = (i==0) ? M.x[i] - 0.5*M.dx[i] : M.x[i-1] + 0.5*M.dx[i-1];

        phiR = tDir;
        phiU = T[i][j];
        phiD = tEsq;
    }
    parte2 = boundaryEsq ? vel*FuncaoUpwind(phiU) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);
    //parte2 = vel*FuncaoUpwind(phiU);

    fluxo1 = (1.0/rI)*(parte1 - parte2)/M.dx[i];


    vel = V[i][j+1];
    coordF = M.y[j+1];
    if( vel>=0 ) {
        coordR = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
        coordU = M.y[j] + 0.5*M.dy[j];
        coordD = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];

        phiR = tBaixo;
        phiU = T[i][j];
        phiD = tCima;
    }
    else {
        coordR = boundaryCima ? lixo : M.y[j+2] + 0.5*M.dy[j+2];
        coordU = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
        coordD = M.y[j] + 0.5*M.dy[j];

        phiR = boundaryCima ? lixo : T[i][j+2]; ///ATENCAO
        phiU = tCima;
        phiD = T[i][j];
    }
    parte1 = boundaryCima ? vel*FuncaoUpwind(phiU) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);
    //parte1 = vel*FuncaoUpwind(phiU);

    vel = V[i][j];
    coordF = M.y[j];
    if( vel>=0 ) {
        coordR = boundaryBaixo ? lixo : M.y[j-2] + 0.5*M.dy[j-2];
        coordU = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
        coordD = M.y[j] + 0.5*M.dy[j];

        phiR = boundaryBaixo ? lixo : T[i][j-2]; ///ATENCAO
        phiU = tBaixo;
        phiD = T[i][j];
    }
    else {
        coordR = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
        coordU = M.y[j] + 0.5*M.dy[j];
        coordD = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];

        phiR = tCima;
        phiU = T[i][j];
        phiD = tBaixo;
    }
    parte2 = boundaryBaixo ? vel*FuncaoUpwind(phiU) : vel*FuncaoCubista(phiR, phiU, phiD, coordF, coordR, coordU, coordD);
    //parte2 = vel*FuncaoUpwind(phiU);

    fluxo2 = (parte1 - parte2)/M.dy[j];

//    if( !strcmp(Tipo, "tyt") && (j==0 || j==M.Ny-1) )
//        PrintDebug("conv (%2d, %2d): %20.15lf %20.15lf %20.15lf %20.15lf %20.15lf\n", i, j, tBaixo, tCima, tEsq, tDir, fluxo1 + fluxo2);

    /// Atencao: No caso axissimetrico tem mais um termo que eu vou colocar direto la na funcao EquacoesConstitutivas
    /// Eh o termo multiplicando (U_theta/r)
    return fluxo1 + fluxo2;
}

// double TensorParede(MALHA M, int i, int j, const char *Tipo, char Face)
// {
//     if( Face=='u' ) {
//         if( !strcmp(Tipo, "txx") )
//             return M.pontosU[i][j].tensorParedeNovo[1][1];
//         if( !strcmp(Tipo, "txy") )
//             return M.pontosU[i][j].tensorParedeNovo[1][2];
//         if( !strcmp(Tipo, "tyy") )
//             return M.pontosU[i][j].tensorParedeNovo[2][2];
//         if( !strcmp(Tipo, "txt") )
//             return M.pontosU[i][j].tensorParedeNovo[1][3];
//         if( !strcmp(Tipo, "tyt") )
//             return M.pontosU[i][j].tensorParedeNovo[2][3];
//         if( !strcmp(Tipo, "ttt") )
//             return M.pontosU[i][j].tensorParedeNovo[3][3];

//         /// Log conformation
//         double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
//         int nrot;
//         int iT, jT;
//         for( iT=1; iT<4; iT++ )
//             for( jT=1; jT<4; jT++ )
//                 A[iT][jT] = M.pontosU[i][j].tensorParedeNovo[iT][jT];
//         A[1][1] += 1.0;
//         A[2][2] += 1.0;
//         A[3][3] += 1.0;
//         CalculaAutovalores_Jacobi(A, 3, auto_valores, aux1, aux2, O, &nrot);

//         if( !strcmp(Tipo, "PsiXX") )
//             return O[1][1]*O[1][1]*log(auto_valores[1]) + O[1][2]*O[1][2]*log(auto_valores[2]) + O[1][3]*O[1][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiXY") )
//             return O[1][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[1][2]*log(auto_valores[2]) + O[2][3]*O[1][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiYY") )
//             return O[2][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[2][2]*log(auto_valores[2]) + O[2][3]*O[2][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiXT") )
//             return O[1][1]*O[3][1]*log(auto_valores[1]) + O[1][2]*O[3][2]*log(auto_valores[2]) + O[1][3]*O[3][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiYT") )
//             return O[2][1]*O[3][1]*log(auto_valores[1]) + O[2][2]*O[3][2]*log(auto_valores[2]) + O[2][3]*O[3][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiTT") )
//             return O[3][1]*O[3][1]*log(auto_valores[1]) + O[3][2]*O[3][2]*log(auto_valores[2]) + O[3][3]*O[3][3]*log(auto_valores[3]);

//     }
//     else if( Face=='v' ) {
//         if( !strcmp(Tipo, "txx") )
//             return M.pontosV[i][j].tensorParedeNovo[1][1];
//         if( !strcmp(Tipo, "txy") )
//             return M.pontosV[i][j].tensorParedeNovo[1][2];
//         if( !strcmp(Tipo, "tyy") )
//             return M.pontosV[i][j].tensorParedeNovo[2][2];
//         if( !strcmp(Tipo, "txt") )
//             return M.pontosV[i][j].tensorParedeNovo[1][3];
//         if( !strcmp(Tipo, "tyt") )
//             return M.pontosV[i][j].tensorParedeNovo[2][3];
//         if( !strcmp(Tipo, "ttt") )
//             return M.pontosV[i][j].tensorParedeNovo[3][3];

//         /// Log conformation
//         double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
//         int nrot;
//         int iT, jT;
//         for( iT=1; iT<4; iT++ )
//             for( jT=1; jT<4; jT++ )
//                 A[iT][jT] = M.pontosV[i][j].tensorParedeNovo[iT][jT];
//         A[1][1] += 1.0;
//         A[2][2] += 1.0;
//         A[3][3] += 1.0;
//         CalculaAutovalores_Jacobi(A, 3, auto_valores, aux1, aux2, O, &nrot);

//         if( !strcmp(Tipo, "PsiXX") )
//             return O[1][1]*O[1][1]*log(auto_valores[1]) + O[1][2]*O[1][2]*log(auto_valores[2]) + O[1][3]*O[1][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiXY") )
//             return O[1][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[1][2]*log(auto_valores[2]) + O[2][3]*O[1][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiYY") )
//             return O[2][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[2][2]*log(auto_valores[2]) + O[2][3]*O[2][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiXT") )
//             return O[1][1]*O[3][1]*log(auto_valores[1]) + O[1][2]*O[3][2]*log(auto_valores[2]) + O[1][3]*O[3][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiYT") )
//             return O[2][1]*O[3][1]*log(auto_valores[1]) + O[2][2]*O[3][2]*log(auto_valores[2]) + O[2][3]*O[3][3]*log(auto_valores[3]);
//         if( !strcmp(Tipo, "PsiTT") )
//             return O[3][1]*O[3][1]*log(auto_valores[1]) + O[3][2]*O[3][2]*log(auto_valores[2]) + O[3][3]*O[3][3]*log(auto_valores[3]);
//     }

//     DeuProblema("TensorParede: tipo nao reconhecido %s\n", Tipo);
//     return 0.0;
// }

double FuncaoCubista(double PhiR, double PhiU, double PhiD, double Xf, double Xr, double Xu, double Xd)
{
    double phiChapeu, XuChapeu, XfChapeu, result = 0.0;
    double condicao1, condicao2;

    //return FuncaoConvCentral(PhiU, PhiD, Xf, Xu, Xd);
    //return FuncaoUpwind(PhiU);

    //Retornando um Upwind nesse caso
    if( fabs(PhiD - PhiR)<1e-13 )
        return PhiU;


    phiChapeu = (PhiU - PhiR)/(PhiD - PhiR);
    XuChapeu = (Xu - Xr)/(Xd - Xr);
    XfChapeu = (Xf - Xr)/(Xd - Xr);

    condicao1 = 0.75*XuChapeu;
    condicao2 = ((1.0+2.0*(XfChapeu-XuChapeu))/(2.0*XfChapeu - XuChapeu))*XuChapeu;

    if( phiChapeu<=0.0 || phiChapeu>=1.0 )
        result = phiChapeu;
    else if( phiChapeu<condicao1 )
        result = phiChapeu*(XfChapeu/XuChapeu)*( 1.0 + ((XfChapeu-XuChapeu)/(3.0*(1.0-XuChapeu))) );
    else if( phiChapeu<=condicao2 )
        result = (phiChapeu*((XfChapeu*(1.0-XfChapeu))/(XuChapeu*(1.0-XuChapeu)))) + ((XfChapeu*(XfChapeu-XuChapeu))/(1.0-XuChapeu));
    else if( phiChapeu<1.0 )
        result = 1.0 - ((1.0-phiChapeu)*((1.0-XfChapeu)/(2.0*(1.0-XuChapeu))));
    else
        DeuProblema("\n\nNAO ERA PRA TER CHEGADO AQUI: PROBLEMA NO CUBISTA %lf!!!\n\n", phiChapeu);

    return (result*(PhiD - PhiR)) + PhiR;
}

double FuncaoUpwind(double PhiU)
{
    return PhiU;
}

double FuncaoConvCentral(double PhiU, double PhiD, double Xf, double Xu, double Xd)
{
    return (PhiD*(Xf-Xu)/(Xd-Xu)) + (PhiU*(Xf-Xd)/(Xu-Xd));
}

double Interpolacao(double h1, double h2, double Phi1, double Phi2)
{
    return (Phi1*h2/(h1+h2)) + (Phi2*h1/(h1+h2)); //Interpolacao
}

double ValorInflowTensor(int i, int j, double **U, const char *Tipo, MALHA M, char Direcao)
{
    static double uCima = 0.0, uBaixo = 0.0, vDir = 0.0, vEsq = 0.0, del = 0.0, tensor = 0.0;
    static double h1 = 0.0, h2 = 0.0;

    if( Direcao=='v' ) {
        if( M.pontosV[i][j].tipo==BOUNDARY ) {
            if( M.pontosV[i][j].tipoBoundary==NOSLIP )
                uBaixo = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], 0.0, U );
            else if( M.pontosV[i][j].tipoBoundary==SLIP )
                uBaixo = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], 0.0, U );
            else if( M.pontosV[i][j].tipoBoundary==SIMETRIA )
                uBaixo = U[i][j];
            else DeuProblema("\n 1. ValorInflowTensor: Problema\n");
        }
        else
            uBaixo = U[i][j-1];

        if( M.pontosV[i][j+1].tipo==BOUNDARY ) {
            if( M.pontosV[i][j+1].tipoBoundary==NOSLIP )
                uCima = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], 0.0, U );
            else if( M.pontosV[i][j+1].tipoBoundary==SLIP )
                uCima = - U[i][j] + 2*FuncaoValorNoSlipHorizontal( M, i, j, M.x[i], 0.0, U );
            else if( M.pontosV[i][j+1].tipoBoundary==SIMETRIA )
                uCima = U[i][j];
            else DeuProblema("\n 2 ValorInflowTensor: Problema\n");
        }
        else
            uCima = U[i][j+1];

        h1 = (j==0) ? M.dy[j] : 0.5*( M.dy[j] + M.dy[j-1] );
        h2 = (j==M.Ny-1) ? M.dy[j] : 0.5*( M.dy[j] + M.dy[j+1] );
        del = ( - uBaixo*h2/(h1*(h1+h2)) ) + ( U[i][j]*(h2-h1)/(h1*h2) ) + ( uCima*h1/(h2*(h1+h2)) );
    }
    /// Neste caso o parametro "U" deve ser a velocidade vertical (que eh a V)
    else if( Direcao=='h' ) {
        if( M.pontosU[i][j].tipo==BOUNDARY ) {
            if( M.pontosU[i][j].tipoBoundary==NOSLIP )
                vEsq = - U[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, U );
            else if( M.pontosU[i][j].tipoBoundary==SLIP )
                vEsq = - U[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, U );
            else if( M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA )
                vEsq = U[i][j];
            else DeuProblema("\n 3. ValorInflowTensor: Problema %d %d\n", i, j);
        }
        else
            vEsq = U[i-1][j];

        if( M.pontosU[i+1][j].tipo==BOUNDARY ) {
            if( M.pontosU[i+1][j].tipoBoundary==NOSLIP )
                vDir = - U[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, U );
            else if( M.pontosU[i+1][j].tipoBoundary==SLIP )
                vDir = - U[i][j] + 2*FuncaoValorNoSlipVertical( M, i, j, M.y[j], 0.0, U );
            else DeuProblema("\n 4 ValorInflowTensor: Problema\n");
        }
        else
            vDir = U[i+1][j];

        h1 = (i==0) ? M.dx[i] : 0.5*( M.dx[i] + M.dx[i-1] );
        h2 = (i==M.Nx-1) ? M.dx[i] : 0.5*( M.dx[i] + M.dx[i+1] );
        del = ( - vEsq*h2/(h1*(h1+h2)) ) + ( U[i][j]*(h2-h1)/(h1*h2) ) + ( vDir*h1/(h2*(h1+h2)) );
    }
    else DeuProblema("ValorInflowTensor: direcao indefinida.\n\n");


    /// === canal horizontal...
    if( Direcao=='v' ) { /// === Canal horizontal
        if( !strcmp(Tipo, "txx") )
            return 2.0*(1.0 - M.beta)*(M.We)*del*del;
        else if( !strcmp(Tipo, "txy") )
            return (1.0 - M.beta)*del;
        else if( !strcmp(Tipo, "tyy") )
            return 0;
        else if( !strcmp(Tipo, "txt") ) /// FALTA VER QUAL COLOCAR AQUI
            return 0;
        else if( !strcmp(Tipo, "tyt") ) /// FALTA VER QUAL COLOCAR AQUI
            return 0;
        else if( !strcmp(Tipo, "ttt") ) /// FALTA VER QUAL COLOCAR AQUI
            return 0;
        else if( !strcmp(Tipo, "lambda") ) {
            tensor = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            return 1.0 + ((M.We/(1.0 - M.beta))*tensor);
        }
        else if( !strcmp(Tipo, "mu") ) {
            tensor = (1.0 - M.beta)*del;
            return (M.We/(1.0 - M.beta))*tensor;
        }
        else if( !strcmp(Tipo, "nu") ) {
            tensor = 0.0;
            return 1.0 + ((M.We/(1.0 - M.beta))*tensor);
        }
        else if( !strcmp(Tipo, "lambdaAux") ) {
            tensor = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            tensor = 1.0 + ((M.We/(1.0 - M.beta))*tensor);
            return tensor/( (U[i][j]*U[i][j]) + M.tolNSF );
        }
        else if( !strcmp(Tipo, "nuAux") ) {
            tensor = 0.0;
            tensor = 1.0 + ((M.We/(1.0 - M.beta))*tensor);
            return U[i][j]*U[i][j]*tensor;
        }
        else if( !strcmp(Tipo, "PsiXX") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            double txy = (1.0 - M.beta)*del;
            double tyy = 0.0;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;

            CalculaAutovalores_Jacobi(A, 2, auto_valores, aux1, aux2, O, &nrot);
            return O[1][1]*O[1][1]*log(auto_valores[1]) + O[1][2]*O[1][2]*log(auto_valores[2]);
        }
        else if( !strcmp(Tipo, "PsiXY") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            double txy = (1.0 - M.beta)*del;
            double tyy = 0.0;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;

            CalculaAutovalores_Jacobi(A, 2, auto_valores, aux1, aux2, O, &nrot);
            return O[1][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[1][2]*log(auto_valores[2]);
        }
        else if( !strcmp(Tipo, "PsiYY") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            double txy = (1.0 - M.beta)*del;
            double tyy = 0.0;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;

            CalculaAutovalores_Jacobi(A, 2, auto_valores, aux1, aux2, O, &nrot);
            return O[2][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[2][2]*log(auto_valores[2]);
        }
        else DeuProblema("\n 3 ValorInflowTensor: Problema %s\n", Tipo);
    }
    else if( Direcao=='h' ) { /// === Canal vertical...
        if( !strcmp(Tipo, "tyy") )
            return 2.0*(1.0 - M.beta)*(M.We)*del*del;
        else if( !strcmp(Tipo, "txy") )
            return (1.0 - M.beta)*del;
        else if( !strcmp(Tipo, "txx") )
            return 0;
        else if( !strcmp(Tipo, "txt") ) /// FALTA VER QUAL COLOCAR AQUI
            return 0;
        else if( !strcmp(Tipo, "tyt") ) /// FALTA VER QUAL COLOCAR AQUI
            return 0;
        else if( !strcmp(Tipo, "ttt") ) /// FALTA VER QUAL COLOCAR AQUI
            return 0;
        else if( !strcmp(Tipo, "lambda") ) {
            tensor = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            return 1.0 + ((M.We/(1.0 - M.beta))*tensor);
        }
        else if( !strcmp(Tipo, "mu") ) {
            tensor = (1.0 - M.beta)*del;
            return (M.We/(1.0 - M.beta))*tensor;
        }
        else if( !strcmp(Tipo, "nu") ) {
            tensor = 0.0;
            return 1.0 + ((M.We/(1.0 - M.beta))*tensor);
        }
        else if( !strcmp(Tipo, "lambdaAux") ) {
            tensor = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            tensor = 1.0 + ((M.We/(1.0 - M.beta))*tensor);
            return tensor/( (U[i][j]*U[i][j]) + M.tolNSF );
        }
        else if( !strcmp(Tipo, "nuAux") ) {
            tensor = 0.0;
            tensor = 1.0 + ((M.We/(1.0 - M.beta))*tensor);
            return U[i][j]*U[i][j]*tensor;
        }
        else if( !strcmp(Tipo, "PsiXX") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 0.0;
            double txy = (1.0 - M.beta)*del;
            double tyy = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;
            if( M.tipoCoord==AXI_CILINDRICO ) {
                A[1][3] = A[3][1] = 0.0/xi;
                A[2][3] = A[3][2] = 0.0/xi;
                A[3][3] = 1.0 + 0.0/xi;

                CalculaAutovalores_Jacobi(A, 3, auto_valores, aux1, aux2, O, &nrot);
                return O[1][1]*O[1][1]*log(auto_valores[1]) + O[1][2]*O[1][2]*log(auto_valores[2]) + O[1][3]*O[1][3]*log(auto_valores[3]);
            }
            // Caso cartesiano 2D
            CalculaAutovalores_Jacobi(A, 2, auto_valores, aux1, aux2, O, &nrot);
            return O[1][1]*O[1][1]*log(auto_valores[1]) + O[1][2]*O[1][2]*log(auto_valores[2]);
        }
        else if( !strcmp(Tipo, "PsiXY") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 0.0;
            double txy = (1.0 - M.beta)*del;
            double tyy = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;
            if( M.tipoCoord==AXI_CILINDRICO ) {
                A[1][3] = A[3][1] = 0.0/xi;
                A[2][3] = A[3][2] = 0.0/xi;
                A[3][3] = 1.0 + 0.0/xi;

                CalculaAutovalores_Jacobi(A, 3, auto_valores, aux1, aux2, O, &nrot);
                return O[1][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[1][2]*log(auto_valores[2]) + O[2][3]*O[1][3]*log(auto_valores[3]);
            }
            // Caso cartesiano 2D
            CalculaAutovalores_Jacobi(A, 2, auto_valores, aux1, aux2, O, &nrot);
            return O[1][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[1][2]*log(auto_valores[2]);
        }
        else if( !strcmp(Tipo, "PsiYY") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 0.0;
            double txy = (1.0 - M.beta)*del;
            double tyy = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;
            if( M.tipoCoord==AXI_CILINDRICO ) {
                A[1][3] = A[3][1] = 0.0/xi;
                A[2][3] = A[3][2] = 0.0/xi;
                A[3][3] = 1.0 + 0.0/xi;

                CalculaAutovalores_Jacobi(A, 3, auto_valores, aux1, aux2, O, &nrot);
                return O[2][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[2][2]*log(auto_valores[2]) + O[2][3]*O[2][3]*log(auto_valores[3]);
            }
            // Caso cartesiano 2D
            CalculaAutovalores_Jacobi(A, 2, auto_valores, aux1, aux2, O, &nrot);
            return O[2][1]*O[2][1]*log(auto_valores[1]) + O[2][2]*O[2][2]*log(auto_valores[2]);
        }
        else if( !strcmp(Tipo, "PsiXT") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 0.0;
            double txy = (1.0 - M.beta)*del;
            double tyy = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;
            if( M.tipoCoord==AXI_CILINDRICO ) {
                A[1][3] = A[3][1] = 0.0/xi;
                A[2][3] = A[3][2] = 0.0/xi;
                A[3][3] = 1.0 + 0.0/xi;

                CalculaAutovalores_Jacobi(A, 3, auto_valores, aux1, aux2, O, &nrot);
                return O[1][1]*O[3][1]*log(auto_valores[1]) + O[1][2]*O[3][2]*log(auto_valores[2]) + O[1][3]*O[3][3]*log(auto_valores[3]);
            }
        }
        else if( !strcmp(Tipo, "PsiYT") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 0.0;
            double txy = (1.0 - M.beta)*del;
            double tyy = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;
            if( M.tipoCoord==AXI_CILINDRICO ) {
                A[1][3] = A[3][1] = 0.0/xi;
                A[2][3] = A[3][2] = 0.0/xi;
                A[3][3] = 1.0 + 0.0/xi;

                CalculaAutovalores_Jacobi(A, 3, auto_valores, aux1, aux2, O, &nrot);
                return O[2][1]*O[3][1]*log(auto_valores[1]) + O[2][2]*O[3][2]*log(auto_valores[2]) + O[2][3]*O[3][3]*log(auto_valores[3]);
            }
        }
        else if( !strcmp(Tipo, "PsiTT") ) {
            double xi = (1.0 - M.beta)/(/*M.Re**/M.We);
            double A[4][4], O[4][4], auto_valores[4], aux1[4], aux2[4];
            int nrot;
            double txx = 0.0;
            double txy = (1.0 - M.beta)*del;
            double tyy = 2.0*(1.0 - M.beta)*(M.We)*del*del;
            A[1][1] = 1.0 + txx/xi;
            A[1][2] = txy/xi;
            A[2][1] = txy/xi;
            A[2][2] = 1.0 + tyy/xi;
            if( M.tipoCoord==AXI_CILINDRICO ) {
                A[1][3] = A[3][1] = 0.0/xi;
                A[2][3] = A[3][2] = 0.0/xi;
                A[3][3] = 1.0 + 0.0/xi;

                CalculaAutovalores_Jacobi(A, 3, auto_valores, aux1, aux2, O, &nrot);
                return O[3][1]*O[3][1]*log(auto_valores[1]) + O[3][2]*O[3][2]*log(auto_valores[2]) + O[3][3]*O[3][3]*log(auto_valores[3]);
            }
        }
        else DeuProblema("\n 3 ValorInflowTensor: Problemaaa %s\n", Tipo);
    }

    DeuProblema("\n 4 ValorInflowTensor: Problema\n");
    return 1e+10; //Lixo. Nao eh pra chegar aqui
}






void CalculaAutovalores_Jacobi(double a[4][4], int n, double d[4], double b[4], double z[4], double v[4][4], int *nrot)
{
//     PEGUEI ESSA FUNCAO DO FREEFLOW (QUE PEGOU DO NUMERICAL RECIPES IN C)
//     Computes all eigenvalues and eigenvectors of a real symmetric matrix a,
//     which is of size n by n, stored in a physical np by np array.
//     On output, elements of a above the diagonal are
//     destroyed. d returns the eigenvalues of a in its first n elements.
//     v is a matrix with the same logical and physical dimensions as a,
//     whose columns contain, on output, the normalized
//     eigenvectors of a.
//     nrot returns the number of Jacobi rotations that were required.
    int NMAX;
    int i, ip, iq, j;
    double c, g, h, s, sm, t, tau, theta, tresh;


    NMAX = 500;
    for (ip = 1; ip <= n; ip++) {
        for (iq = 1; iq <= n; iq++) {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
        b[ip]     = a[ip][ip];
        d[ip]     = b[ip];
        z[ip]     = 0.0;
    }
    *nrot = 0;
    for (i = 1; i <= NMAX; i++) {
        sm = 0.0;
        for (ip = 1;  ip <= n - 1; ip++)
            for (iq = ip + 1; iq <= n; iq++)
                sm += fabs (a[ip][iq]);
        if (sm == 0.0) return;
        if (i < 4)
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0.0;
        for (ip = 1; ip <= n - 1; ip++)
            for (iq = ip + 1; iq <= n; iq++) {
                g = 100.0 * fabs (a[ip][iq]);
                if ( (i > 4) && (fabs (d[ip]) + g == fabs (d[ip]) ) &&
                        (fabs (d[iq]) + g == fabs (d[iq]) ) ) {
                    a[ip][iq] = 0.0;
                } else if (fabs (a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if (fabs (h) + g ==  fabs (h) ) {
                        t = a[ip][iq] / h;
                    } else {
                        theta = 0.5 * h / a[ip][iq];
                        t     = 1.0 / (fabs (theta) + sqrt (1.0 + theta * theta) );
                        if (theta < 0.0) t = -t;
                    }
                    c = 1.0 / sqrt (1.0 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h   = t * a[ip][iq];
                    z[ip] = z[ip] - h;
                    z[iq] = z[iq] + h;
                    d[ip] = d[ip] - h;
                    d[iq] = d[iq] + h;
                    a[ip][iq] = 0.0;
                    for (j = 1; j <= ip - 1; j++) {
                        g = a[j][ip];
                        h = a[j][iq];
                        a[j][ip] = g - s * (h + g * tau);
                        a[j][iq] = h + s * (g - h * tau);
                    }
                    for (j = ip + 1; j <= iq - 1; j++) {
                        g = a[ip][j];
                        h = a[j][iq];
                        a[ip][j] = g - s * (h + g * tau);
                        a[j][iq] = h + s * (g - h * tau);
                    }
                    for (j = iq + 1; j <= n; j++) {
                        g = a[ip][j];
                        h = a[iq][j];
                        a[ip][j] = g - s * (h + g * tau);
                        a[iq][j] = h + s * (g - h * tau);
                    }
                    for (j = 1; j <= n; j++) {
                        g = v[j][ip];
                        h = v[j][iq];
                        v[j][ip] = g - s * (h + g * tau);
                        v[j][iq] = h + s * (g - h * tau);
                    }
                    *nrot = *nrot + 1;
                }
            }
        for (ip = 1; ip <= n; ip++) {
            b[ip] = b[ip] + z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    return;
}
