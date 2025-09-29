#include "ContornoSuperficie.h"

extern int velocidadeNewtoniana;

LISTA_CELULAS_SURFACE ClassificaCelulasEPontos(MALHA *M, CURVA *C, double **U, double **V, double **P)
{
    int i, j, iMedio, jMedio;
    PONTO *ponto;
    CELULA_SURFACE *celulaSurf, *proxCelula;
    LISTA_CELULAS_SURFACE listaCelulas;

    listaCelulas.prim = NULL;
    listaCelulas.ult = NULL;

    ///Primeiramente vou calcular os vetores normais em todas particulas dessa curva
    ///Soh vai usar na formulacao com normal exata
    //for( ponto=C->nodeInicial; ponto!=NULL; ponto=ponto->prox )
    //    VetorNormal(C, ponto, &(ponto->vetNormalX), &(ponto->vetNormalY));


    ///Agora vamos classificar as celulas como SURFACE
    //Colocando SURFACE em todas as celulas que tem pontos
    LOOP_CURVA(C) {
        ponto = loop_stuff.ponto;

        EncontraCelulaDaParticula(M, ponto->x, ponto->y, &i, &j);

        if( i<0 || i>=M->Nx || j<0 || j>=M->Ny )
            DeuProblema("\n\n CELULA FORA DO DOMINIO %lf %lf %d %d\n\n", ponto->x, ponto->y, i, j);

        // if( fabs(ponto->y - M->yMin)<1e-9 )
        //     continue;

        //Adicionando na lista
        AdicionaCelulaSurface(M, &listaCelulas, i, j, 0);
        M->celulas[i][j] = SURFACE;
    }

    



    //Preenchendo todo o interior com celulas FULL
    EncontraCelulaSeedFloodFill(*M, listaCelulas, &iMedio, &jMedio);
    // PrintDebug("\n\n\n RETIRAR ISTO AQUI DEPOIS!!!!!!!! INICIALIZACAO DO FLOOD_FILL!!! \n\n\n");
    // iMedio = M->Nx/2;
    // jMedio = 0;
    AlgoritmoFloodFillCelulaIterativo(*M, iMedio, jMedio, FULL);



    //Tratando casos degenerados
    celulaSurf = listaCelulas.prim;
    while( celulaSurf!=NULL ) {
        i = celulaSurf->i;
        j = celulaSurf->j;
        proxCelula = celulaSurf->prox;

        //printf("%d %d %lf %lf\n", i, j, ponto->x, ponto->y); fflush(stdout);

        if( M->celulas[i][j]!=SURFACE ) {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 1: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! %d %d\n", i, j);
        }

        //Cantinho de cima
        if( (i!=M->Nx-1 && M->celulas[i+1][j]==EMPTY) && (i!=0 && M->celulas[i-1][j]==EMPTY) && (j==M->Ny-1 || M->celulas[i][j+1]==EMPTY) ) {
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, EMPTY);
            AdicionaCelulaSurface(M, &listaCelulas, i, j-1, 0);
            M->celulas[i][j] = EMPTY;
            M->celulas[i][j-1] = SURFACE;
        }
        else if( (i!=M->Nx-1 && M->celulas[i+1][j]==FULL) && (i!=0 && M->celulas[i-1][j]==FULL) && (j==M->Ny-1 || M->celulas[i][j+1]==FULL) ) {
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, FULL);
            AdicionaCelulaSurface(M, &listaCelulas, i, j-1, 0);
            M->celulas[i][j] = FULL;
            M->celulas[i][j-1] = SURFACE;
        }

        //Cantinho de baixo
        if( (i!=M->Nx-1 && M->celulas[i+1][j]==EMPTY) && (i!=0 && M->celulas[i-1][j]==EMPTY) && (j==0 || M->celulas[i][j-1]==EMPTY) ) {
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, EMPTY);
            AdicionaCelulaSurface(M, &listaCelulas, i, j+1, 0);
            M->celulas[i][j] = EMPTY;
            M->celulas[i][j+1] = SURFACE;
        }
        else if( (i!=M->Nx-1 && M->celulas[i+1][j]==FULL) && (i!=0 && M->celulas[i-1][j]==FULL) && (j==0 || M->celulas[i][j-1]==FULL) ) {
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, FULL);
            AdicionaCelulaSurface(M, &listaCelulas, i, j+1, 0);
            M->celulas[i][j] = FULL;
            M->celulas[i][j+1] = SURFACE;
        }

        //Cantinho da direita
        if( (i==M->Nx-1 || M->celulas[i+1][j]==EMPTY) && (j!=0 && M->celulas[i][j-1]==EMPTY) && (j!=M->Ny-1 && M->celulas[i][j+1]==EMPTY) ) {


            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, EMPTY);
            AdicionaCelulaSurface(M, &listaCelulas, i-1, j, 0);
            M->celulas[i][j] = EMPTY;
            M->celulas[i-1][j] = SURFACE;
        }
        else if( (i==M->Nx-1 || M->celulas[i+1][j]==FULL) && (j==0 || M->celulas[i][j-1]==FULL) && (j==M->Ny-1 || M->celulas[i][j+1]==FULL) ) {
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, FULL);
            AdicionaCelulaSurface(M, &listaCelulas, i-1, j, 0);
            M->celulas[i][j] = FULL;
            M->celulas[i-1][j] = SURFACE;
        }

        //Cantinho da esquerda
        if( (i==0 || M->celulas[i-1][j]==EMPTY) && (j!=0 && M->celulas[i][j-1]==EMPTY) && (j!=M->Ny-1 && M->celulas[i][j+1]==EMPTY) ) {
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, EMPTY);
            AdicionaCelulaSurface(M, &listaCelulas, i+1, j, 0);
            M->celulas[i][j] = EMPTY;
            M->celulas[i+1][j] = SURFACE;
        }
        else if( (i==0 || M->celulas[i-1][j]==FULL) && (j==0 || M->celulas[i][j-1]==FULL) && (j==M->Ny-1 || M->celulas[i][j+1]==FULL) ) {
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, FULL);
            AdicionaCelulaSurface(M, &listaCelulas, i+1, j, 0);
            M->celulas[i][j] = FULL;
            M->celulas[i+1][j] = SURFACE;
        }



        //Celula SURFACE com EMPTY na direita e esquerda
        if( i!=0 && i!=M->Nx-1 && j!=0 && j!=M->Ny-1 && M->celulas[i-1][j]==EMPTY && M->celulas[i+1][j]==EMPTY && M->celulas[i][j-1]!=EMPTY && M->celulas[i][j+1]!=EMPTY ) {
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, EMPTY);
            AdicionaCelulaSurface(M, &listaCelulas, i, j-1, 0);
            M->celulas[i][j] = EMPTY;
            M->celulas[i][j-1] = SURFACE;
        }

        celulaSurf = proxCelula;
    }

    // DesenhaMalhaVTK(*M, 0);
    // DeuProblema("parou\n");


    int criterioEsq, criterioDir, criterioCima, criterioBaixo;
    ///Transformando em FULL as celulas que nao tem nenhuma empty em volta
    celulaSurf = listaCelulas.prim;
    while( celulaSurf!=NULL ) {
        i = celulaSurf->i;
        j = celulaSurf->j;
        proxCelula = celulaSurf->prox;

        if( M->celulas[i][j]!=SURFACE ) {
            DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 2: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! %d %d\n", i, j);
        }

        criterioDir = (i==M->Nx-1) || (M->celulas[i+1][j]!=EMPTY) || (M->pontosU[i+1][j].tipoBoundary==NOSLIP);
        criterioEsq = (i==0 || M->celulas[i-1][j]!=EMPTY) || (M->pontosU[i][j].tipoBoundary==NOSLIP);
        criterioCima = (j==M->Ny-1 || M->celulas[i][j+1]!=EMPTY) || (M->pontosV[i][j+1].tipoBoundary==NOSLIP);
        criterioBaixo = (j==0) || (M->celulas[i][j-1]!=EMPTY) || (M->pontosV[i][j].tipoBoundary==NOSLIP);

        if( criterioDir && criterioEsq && criterioCima && criterioBaixo ) {
            M->celulas[i][j] = FULL;
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, FULL);

            // Testando tensao superficial na celula full
            M->celulasParede[i][j] = SURFACE;
        }
        //else PrintDebug("NAO REMOVEU: %d %d\n", i, j);

        celulaSurf = proxCelula;
    }


    return listaCelulas;
}

void ReclassificaCelulasEPontos(MALHA *M, LISTA_CELULAS_SURFACE *ListaCelulas, int PassoTemporal)
{
    int i = -1, j = -1, linha = -1;
    int iVelho=-1, jVelho=-1;
    PONTO *ponto;
    CURVA *curva;
    CELULA_SURFACE *celulaSurf, *proxCelula;
    int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    //printf("M->x[5]=%f, M->x[6]=%f, M->y[119]=%f, M->y[120]=%f\n", M->x[5], M->x[6], M->y[119], M->y[120]);
    //printf("M->pontosU[5][119].tipo=%d, M->pontosU[6][119].tipo=%d, M->pontosV[5][119].tipo=%d, M->pontosV[5][120].tipo=%d\n", M->pontosU[5][119].tipo, M->pontosU[6][119].tipo, M->pontosV[5][119].tipo, M->pontosV[5][120].tipo);
    //printf("M->celulas[5][119]=%d, M->celulas[4][119]=%d, M->celulas[6][119]=%d, M->celulas[5][118]=%d, M->celulas[5][120]=%d\n", M->celulas[5][119], M->celulas[4][119], M->celulas[6][119], M->celulas[5][118], M->celulas[5][120]);

    	/*
    	int iNovo, jNovo;
	for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva) {
        	 LOOP_CURVA(curva) {
                	 ponto = loop_stuff.ponto;
 
                 //proxCurva = curva->proxCurva;
 
                 //ponto = curva->pontoInicial;
         //do {
             //proxPonto = ponto->prox;
 
             	EncontraCelulaDaParticula(M, ponto->x, ponto->y, &iNovo, &jNovo);
 
             	if (ponto->x >= 0.25 && (ponto->y >= -0.2 && ponto->y <= 0.0)){
                     	printf("\n\nVERIFICA INTERFACE\nParticula marcadora dentro de celula BOUNDARY\n");
                     	printf("iNovo=%d, jNovo=%d, xNovo=%f, yNovo=%f\n", iNovo, jNovo, ponto->x, ponto->y);
                     	exit(1);
             	}
 
             //ponto = proxPonto;
         	} //while( ponto!=curva->pontoInicial );
 
     	}
	*/	


    ///Apenas verificando se a ListaCelulas foi inicializada
    ///Vai ser preciso inicializar aqui quando a funcao CarregarEstado for usada pra debugar coisas
    if( ListaCelulas->prim==NULL && ListaCelulas->ult==NULL ) {
        return;
        //DeuProblema("\n Problema reclassificaCelulasEpontos \n");
    }


    ///Primeiramente vou calcular os vetores normais em todas particulas dessa curva
    ///Soh vai usar na formulacao com normal exata
    //for( ponto=C->nodeInicial; ponto!=NULL; ponto=ponto->prox )
    //    VetorNormal(C, ponto, &(ponto->vetNormalX), &(ponto->vetNormalY));

    /// Na formulacao com tensao superficial em celulas FULL, vou apagar as celulas SURFACE_FULL do passo anterior...
    for( linha=0; linha<M->qtdIncognitasP; linha++ ) {
        i = M->indicesP[linha].i;
        j = M->indicesP[linha].j;

        if( M->celulas[i][j]==FULL && M->celulasParede[i][j]==SURFACE )
            M->celulasParede[i][j] = EMPTY;
    }

    ///Agora vamos classificar as celulas como SURFACE
    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        //Colocando SURFACE em todas as celulas que tem pontos
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;

            //Ignora pontos que passaram pelo outflow
            if( ponto->outflow )
                continue;

            EncontraCelulaDaParticula(M, ponto->x, ponto->y, &i, &j);

	    /*
	    if ( (i>=5 && j<=119) && (i<=40 && j>=116)){
	    if ( (ponto->x >= M->x1) && (ponto->x <= M->x2) && (ponto->y <= M->y1) && (ponto->y >= M->y3) ){
		    EncontraCelulaDaParticula(M, ponto->ant->x, ponto->ant->y, &iVelho, &jVelho);
		    printf("Ha particulas marcadoras em celulas BOUNDARY\n");
		    printf("ponto->ant->x=%f, ponto->ant->y=%f\n", ponto->ant->x, ponto->ant->y);
		    printf("i=%d, j=%d, iVelho=%d, jVelho=%d\n", i, j, iVelho, jVelho);
		    printf("ponto->x=%f, ponto->y=%f\n", ponto->x, ponto->y);//getchar();
		    printf("celulas=%d\n", M->celulas[i][j]);getchar();
	    }
	    */
            if( i<0 || i>=M->Nx || j<0 || j>=M->Ny )
                DeuProblema("\n\n CELULA FORA DO DOMINIO %lf %lf %d %d\n\n", ponto->x, ponto->y, i, j);

            //Adicionando na lista
            AdicionaCelulaSurface(M, ListaCelulas, i, j, PassoTemporal);

            M->celulas[i][j] = SURFACE;
        }
    }

    ///Removendo as celulas SURFACE que nao possuem mais particulas (devido a ultima adveccao)
    celulaSurf = ListaCelulas->prim;
    while( celulaSurf!=NULL ) {
        proxCelula = celulaSurf->prox;

        if( celulaSurf->passoCelula!=PassoTemporal ) {
            i = celulaSurf->i;
            j = celulaSurf->j;
            if( (i==M->Nx-1 || M->celulas[i+1][j]==EMPTY) || (i==0 || M->celulas[i-1][j]==EMPTY) || (j==M->Ny-1 || M->celulas[i][j+1]==EMPTY) || (j==0 || M->celulas[i][j-1]==EMPTY) ) {
                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, EMPTY);
                M->celulas[i][j] = EMPTY;


                /// == FORCANDO CELULA FULL DENTRO DO CANAL NO CANTO
                /// TIRAR ISSO DEPOIS!!!!
//                if( M->pontosU[i+1][j].tipoBoundary==NOSLIP || M->pontosU[i][j].tipoBoundary==NOSLIP)
//                    M->celulas[i][j] = FULL;
            }
        }

        celulaSurf = proxCelula;
    }

    /// Corrigindo casos especificos em que particulas sairam de uma celula FULL_SURFACE
    /// Nesse caso a celula full vira EMPTY
    /// Isso acontece quando tem interacao entre frontes, tipo na colisao de duas gotas
    int continua = 1;
    while( continua ) {
        continua = 0;

        for(i=0; i<M->Nx; i++) {
            for( j=0; j<M->Ny; j++ ) {

                if( M->celulas[i][j]!=FULL )
                    continue;

                criterioEsq = (i!=0) && (M->celulas[i-1][j]==EMPTY) && (M->pontosU[i][j].tipo!=BOUNDARY);
                criterioDir = (i!=M->Nx-1) && (M->celulas[i+1][j]==EMPTY) && (M->pontosU[i+1][j].tipo!=BOUNDARY);
                criterioBaixo = (j!=0) && (M->celulas[i][j-1]==EMPTY) && (M->pontosV[i][j].tipo!=BOUNDARY);
                criterioCima = (j!=M->Ny-1) && (M->celulas[i][j+1]==EMPTY) && (M->pontosV[i][j+1].tipo!=BOUNDARY);

                if( criterioEsq || criterioDir || criterioBaixo || criterioCima ) {
                    M->celulas[i][j] = EMPTY;
                    continua = 1;
                }
            }
        }
    }

    /// === Comecei a fazer isso em casos de splashing
    /// === Vou remover gotinhas que sao TAO pequenas que nao possuem mais celulas FULL (possuem apenas surface)
    /// === Precisaria refinar mais a malha pra conseguir simular essas gotinhas muito pequenas
//    char valorInicial = 0;
//    char **celulasVisitadas = (char **)AlocaMatriz(M->Nx, M->Ny, sizeof(char), &valorInicial);
//    celulaSurf = ListaCelulas->prim;
//    while( celulaSurf!=NULL ) {
//        i = celulaSurf->i;
//        int j = celulaSurf->j;
//        proxCelula = celulaSurf->prox;
//
//        if( celulasVisitadas[i][j] ) {
//            celulaSurf = proxCelula;
//            continue;
//        }
//
//        int qtdFull = EncontraGotasPequenas_Surface(M, i, j, celulasVisitadas);
//        if( qtdFull==0 ) {
//            /// Vamos remover essa gotinha
//            PONTO *p, *proxPonto;
//            CURVA *proxCurva;
//            int iPonto, jPonto;
//
//            curva = M->interface.curvaInicial;
//            while( curva ) {
//                proxCurva = curva->proxCurva;
//                p = curva->nodeInicial;
//
//                while( p ) {
//                    proxPonto = p->prox;
//
//                    EncontraCelulaDaParticula(M, p->x, p->y, &iPonto, &jPonto);
//                    if( iPonto==i && jPonto==j ) {
//                        /// Deletando essa curva da interface
//                        if( curva->curvaAnt )
//                            curva->curvaAnt->proxCurva = proxCurva;
//                        if( proxCurva )
//                            proxCurva->curvaAnt = curva->curvaAnt;
//                        if( M->interface.curvaInicial==curva )
//                            M->interface.curvaInicial = proxCurva;
//                        LiberaMemoriaCurva(curva);
//                        break;
//                    }
//
//                    p = proxPonto;
//                }
//
//                curva = proxCurva;
//            }
//
//            RemoveGotaPequena_Surface(M, i, j, celulasVisitadas, ListaCelulas);
//            proxCelula = ListaCelulas->prim; //Comeca de novo, pois posso ter removido a proxCelula...
//        }
//        celulaSurf = proxCelula;
//    }
//    DesalocaMatriz((void**)celulasVisitadas, M->Nx, M->Ny);


    ///Tratando casos degenerados
    TratandoCasosDegenerados(M, ListaCelulas, PassoTemporal);


    ///Corrigindo um caso que pode ficar celula FULL do lado de uma celula EMPTY (ao inves de surface)
    ///Acontece com muito pouca frequencia e soh em algumas simulacoes
    ///Tentar melhorar isso depois...
    for(i=0; i<M->Nx; i++) {
        for( j=0; j<M->Ny; j++ ) {

            if( M->celulas[i][j]!=FULL )
                continue;



            criterioEsq = (i!=0) && (M->celulas[i-1][j]==EMPTY) && (M->pontosU[i][j].tipo!=BOUNDARY);
            criterioDir = (i!=M->Nx-1) && (M->celulas[i+1][j]==EMPTY) && (M->pontosU[i+1][j].tipo!=BOUNDARY);
            criterioBaixo = (j!=0) && (M->celulas[i][j-1]==EMPTY) && (M->pontosV[i][j].tipo!=BOUNDARY);
            criterioCima = (j!=M->Ny-1) && (M->celulas[i][j+1]==EMPTY) && (M->pontosV[i][j+1].tipo!=BOUNDARY);

            if( criterioEsq || criterioDir || criterioBaixo || criterioCima ) {
                AdicionaCelulaSurface(M, ListaCelulas, i, j, PassoTemporal);
                M->celulas[i][j] = SURFACE;

                DesenhaMalhaVTK(*M, 0);
                DeuProblema("FULL AO LADO DE EMPTY\n");
            }
        }
    }

    TratandoCasosDegenerados(M, ListaCelulas, PassoTemporal);


    ///Preenchendo todo o interior com pontos U, V, P FULL
    for( i=0; i<=M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = EMPTY;
            M->matIndicesU[i][j] = -1;
        }
    }

    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<=M->Ny; j++ ) {
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = EMPTY;
            M->matIndicesV[i][j] = -1;
        }
    }

    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            if( M->pontosP[i][j].tipo!=BOUNDARY )
                M->pontosP[i][j].tipo = EMPTY;
            M->matIndicesP[i][j] = -1;
        }
    }

    ///Primeiro fazendo a borda
    for( celulaSurf=ListaCelulas->prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox ) {
        i = celulaSurf->i;
        j = celulaSurf->j;

        if( M->celulas[i][j]!=SURFACE ) {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 3: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! \n");
        }

        M->pontosP[i][j].tipo = SURFACE;

        //Superficie vertical: empty na direita
        if( LADO_VERTICAL_DIREITA(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        //Superficie vertical: empty na esquerda
        else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        //Superficie horizontal: empty em cima
        else if( LADO_HORIZONTAL_CIMA(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
        }
        //Superficie horizontal: empty embaixo
        else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        //Superficie quina: empty em cima e na direita
        else if( QUINA_CIMA_DIREITA(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
        }
        //Superficie quina: empty em cima e na esquerda
        else if( QUINA_CIMA_ESQUERDA(i, j) ) {
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
        }
        //Superficie quina: empty em embaixo e na direita
        else if( QUINA_BAIXO_DIREITA(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        //Superficie quina: empty em embaixo e na esquerda
        else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        else if( DEGENERADO_CIMA_BAIXO(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
        }
        else if( DEGENERADO_ESQ_DIR(i, j) ) {
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
	//Minha modificacao - 02/09/2025
        else if( DEGENERADO_ESQ_CIMA_BAIXO(i, j) ) {
            if( M->pontosU[i+1][j].tipo==BOUNDARY ){
		//Como tratar esse caso degenerado?
		//M->celulas[i][j] = EMPTY;
		printf("\n\nCASO DEGENERADO NAO TRATADO CORRETAMENTE - vai dar algum problema!!\n");
		printf("Esq-EMPTY, cima-EMPTY, baixo-EMPTY, dir-BOUNDARY e centro-SURFACE  \n\n");
	    }
        }
        else {
            DesenhaMalhaVTK(*M, 0);
//            ImprimeInterfaceVTK(*M, 10000000);
	    printf("M->pontosU[i][j].tipo=%d, M->pontosU[i+1][j].tipo=%d, M->pontosV[i][j].tipo=%d, M->pontosV[i][j+1].tipo=%d\n", M->pontosU[i][j].tipo, M->pontosU[i+1][j].tipo, M->pontosV[i][j].tipo, M->pontosV[i][j+1].tipo);
	    printf("M->celulas[i][j]=%d, M->celulas[i-1][j]=%d, M->celulas[i+1][j]=%d, M->celulas[i][j-1]=%d, M->celulas[i][j+1]=%d\n", M->celulas[i][j], M->celulas[i-1][j], M->celulas[i+1][j], M->celulas[i][j-1], M->celulas[i][j+1]);
	    printf("M->x[i]=%f, M->x[i+1]=%f, M->y[j]=%f, M->y[j+1]=%f", M->x[i], M->x[i+1], M->y[j], M->y[j+1]);
            DeuProblema("\n\n RECLASSIFICACAO DA SUP. LIVRE: CASO INESPERADO %d %d \n\n", i, j);
        }
    }

	  
    //Classifica Bloco Paraview como BOUNDARY (2).
    //Criei essa classificação pois no problema da gota em um orificio 
    //havia a criacao de particulas em celulas BOUNDARY. Data da criacao 28/07/2025.
    double x1, x2, x3, x4;
    double y1, y2, y3, y4;
    double epsilon=1.0e-06;
    //x1 = M->PontosBlocoParaview[0][0];
    //y1 = M->PontosBlocoParaview[0][1];
    //x2 = M->PontosBlocoParaview[1][0];
    //y2 = M->PontosBlocoParaview[1][1];
    //x3 = M->PontosBlocoParaview[2][0];
    //y3 = M->PontosBlocoParaview[2][1];
    //x4 = M->PontosBlocoParaview[3][0];
    //y4 = M->PontosBlocoParaview[3][1];
    x1 = M->x1;
    x2 = M->x2;
    x3 = M->x3;
    x4 = M->x4;
    y1 = M->y1;
    y2 = M->y2;
    y3 = M->y3;
    y4 = M->y4;
    //printf("Let's print\n");
    //printf("x1=%lf, x2=%lf, x3=%lf, x4=%lf\n", x1, x2, x3, x4);
    //printf("y1=%lf, y2=%lf, y3=%lf, y4=%lf\n", y1, y2, y3, y4);getchar();
    for(i=0; i<M->Nx; i++) {
    	for( j=0; j<M->Ny; j++ ) {
        	if ( ( (M->x[i]>= (x1-epsilon)) && (M->x[i+1]<=(x2+epsilon)) ) && ( (M->y[j+1] <= (y2+epsilon)) && (M->y[j] >= (y3-epsilon)) ) ){
			//printf("x[i]=%f, y[j]=%f, x[i+1]=%f, y[j+1]=%f", M->x[i], M->y[j], M->x[i+1], M->y[j+1]);getchar();
			M->celulas[i][j] = BOUNDARY;
		}
	}
    }
    //-----		

    int iMedio = -1, jMedio = -1, iPressao = -1, jPressao = -1;
    do {
        EncontraVelPressaoSeedFloodFill(*M, &iMedio, &jMedio, &i, &j, &iPressao, &jPressao);
        AlgoritmoFloodFillUIterativo(*M, iMedio, jMedio, FULL);
        AlgoritmoFloodFillVIterativo(*M, i, j, FULL);
        AlgoritmoFloodFillPIterativo(*M, iPressao, jPressao, FULL);
    } while( iMedio!=-1 || jMedio!=-1 || i!=-1 || j!=-1 || iPressao!=-1 || jPressao!=-1 );

    return;
}

void InicializaNovasPropriedadesSURFACE(MALHA *M, LISTA_CELULAS_SURFACE *ListaCelulas, double **EVPT_Structure)
{
    CELULA_SURFACE *celula;

    return;

    for( celula=ListaCelulas->prim; celula; celula=celula->prox ) {
        if( !celula->novaCelulaSurface )
            continue;
        
        int i = celula->i;
        int j = celula->j;

        // Extrapolating from the FULL cells around
        double valor = 0.0;
        int numCelulasFULL = 0;
        valor += (M->celulas[i+1][j]==FULL) ? EVPT_Structure[i+1][j] : 0.0;
        valor += (M->celulas[i-1][j]==FULL) ? EVPT_Structure[i-1][j] : 0.0;
        valor += (M->celulas[i][j+1]==FULL) ? EVPT_Structure[i][j+1] : 0.0;
        valor += (M->celulas[i][j-1]==FULL) ? EVPT_Structure[i][j-1] : 0.0;
        valor += (M->celulas[i+1][j+1]==FULL) ? EVPT_Structure[i+1][j+1] : 0.0;
        valor += (M->celulas[i+1][j-1]==FULL) ? EVPT_Structure[i+1][j-1] : 0.0;
        valor += (M->celulas[i-1][j+1]==FULL) ? EVPT_Structure[i-1][j+1] : 0.0;
        valor += (M->celulas[i-1][j-1]==FULL) ? EVPT_Structure[i-1][j-1] : 0.0;
        numCelulasFULL += (M->celulas[i+1][j]==FULL) ? 1 : 0;
        numCelulasFULL += (M->celulas[i-1][j]==FULL) ? 1 : 0;
        numCelulasFULL += (M->celulas[i][j+1]==FULL) ? 1 : 0;
        numCelulasFULL += (M->celulas[i][j-1]==FULL) ? 1 : 0;
        numCelulasFULL += (M->celulas[i+1][j+1]==FULL) ? 1 : 0;
        numCelulasFULL += (M->celulas[i+1][j-1]==FULL) ? 1 : 0;
        numCelulasFULL += (M->celulas[i-1][j+1]==FULL) ? 1 : 0;
        numCelulasFULL += (M->celulas[i-1][j-1]==FULL) ? 1 : 0;

        EVPT_Structure[i][j] = valor/numCelulasFULL;
    }


    return;
}

void TratandoCasosDegenerados(MALHA *M, LISTA_CELULAS_SURFACE *ListaCelulas, int PassoTemporal)
{
    int i, j;
    int criterioEsq, criterioDir, criterioCima, criterioBaixo;
    CELULA_SURFACE *celulaSurf, *proxCelula;


    int continua = 1;
    while( continua ) {
        continua = 0;

        celulaSurf = ListaCelulas->prim;
        while( celulaSurf!=NULL ) {
            i = celulaSurf->i;
            j = celulaSurf->j;
            proxCelula = celulaSurf->prox;


            if( M->celulas[i][j]!=SURFACE ) {
                //DesenhaMalhaVTK(*M, 0);
                DeuProblema("\n 1: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! %d %d\n", i, j);
            }


            // Caso de quina em que aparece uma celula surface completamente cercada de celulas EMPTY
            // Eh bem raro de acontecer
            criterioEsq = (i!=0) && (M->celulas[i-1][j]==EMPTY);
            criterioDir = (i!=M->Nx-1) && (M->celulas[i+1][j]==EMPTY);
            criterioBaixo = (j!=0) && (M->celulas[i][j-1]==EMPTY);
            criterioCima = (j!=M->Ny-1) && (M->celulas[i][j+1]==EMPTY);
            if( criterioEsq && criterioDir && criterioBaixo && criterioCima ) {


                // Quina cima-esq (adiciona uma surface a direita
                if( M->celulas[i+1][j-1]==SURFACE ) {
                    if( M->pontosV[i][j].tipo==BOUNDARY ) {
                        AdicionaCelulaSurface(M, ListaCelulas, i+1, j, PassoTemporal);
                        M->celulas[i+1][j] = SURFACE;
                        continua = 1;
                    }
                    else if( M->pontosU[i+1][j].tipo==BOUNDARY ) {
                        AdicionaCelulaSurface(M, ListaCelulas, i, j-1, PassoTemporal);
                        M->celulas[i][j-1] = SURFACE;
                        continua = 1;
                    }
                    else {
//                        DesenhaMalhaVTK(*M, 0);
//                        ImprimeInterfaceVTK(*M, 1000000);
                        DeuProblema("RECLASSIFICACAO CELULAS. CASO DEGENERADO QUINAS 1!! %d %d\n\n", i, j);
                    }
                }
                // Quina cima-dir (adiciona uma surface a esquerda
                else if( M->celulas[i-1][j-1]==SURFACE ) {
                    if( M->pontosV[i][j].tipo==BOUNDARY ) {
                        AdicionaCelulaSurface(M, ListaCelulas, i-1, j, PassoTemporal);
                        M->celulas[i-1][j] = SURFACE;
                        continua = 1;
                    }
                    else if( M->pontosU[i][j].tipo==BOUNDARY ) {
                        AdicionaCelulaSurface(M, ListaCelulas, i, j-1, PassoTemporal);
                        M->celulas[i][j-1] = SURFACE;
                        continua = 1;
                    }
                    else DeuProblema("RECLASSIFICACAO CELULAS. CASO DEGENERADO QUINAS 2!!\n\n");
                }
                // Quina baixo-esq (adiciona uma surface a direita
                else if( M->celulas[i+1][j+1]==SURFACE ) {

                    if( M->pontosV[i][j+1].tipo==BOUNDARY ) {
                        AdicionaCelulaSurface(M, ListaCelulas, i+1, j, PassoTemporal);
                        M->celulas[i+1][j] = SURFACE;
                        continua = 1;
                    }
                    else if( M->pontosU[i+1][j].tipo==BOUNDARY ) {
                        AdicionaCelulaSurface(M, ListaCelulas, i, j+1, PassoTemporal);
                        M->celulas[i][j+1] = SURFACE;
                        continua = 1;
                    }
                    else DeuProblema("RECLASSIFICACAO CELULAS. CASO DEGENERADO QUINAS 3!!\n\n");
                }
                // Quina baixo-dir (adiciona uma surface a esquerda
                else if( M->celulas[i-1][j+1]==SURFACE ) {

                    if( M->pontosV[i][j+1].tipo==BOUNDARY ) {
                        AdicionaCelulaSurface(M, ListaCelulas, i-1, j, PassoTemporal);
                        M->celulas[i-1][j] = SURFACE;
                        continua = 1;
                    }
                    else if( M->pontosU[i][j].tipo==BOUNDARY ) {
                        AdicionaCelulaSurface(M, ListaCelulas, i, j+1, PassoTemporal);
                        M->celulas[i][j+1] = SURFACE;
                        continua = 1;
                    }
                    else {
                        DesenhaMalhaVTK(*M, 0);
//                        ImprimeInterfaceVTK(*M, 100000);
                        DeuProblema("RECLASSIFICACAO CELULAS. CASO DEGENERADO QUINAS 4!!\n\n");
                    }
                }


            }


            //Cantinho de cima
            if( (i!=M->Nx-1 && M->celulas[i+1][j]==EMPTY) && (i!=0 && M->celulas[i-1][j]==EMPTY) && (j==M->Ny-1 || M->celulas[i][j+1]==EMPTY) && (M->pontosV[i][j].tipo!=BOUNDARY) ) {

                /// Caso a esquerda desta celula tenha uma parede, nao eh preciso corrigir este caso (nao eh degenerado)
                if( M->pontosU[i][j].tipo!=BOUNDARY && M->pontosU[i+1][j].tipo!=BOUNDARY ) {
                    if( i==0 )
                        PrintDebug("AQUI 1\n");
                    RemoveCelulaSurface(ListaCelulas, celulaSurf, M, EMPTY);
                    AdicionaCelulaSurface(M, ListaCelulas, i, j-1, PassoTemporal);
                    M->celulas[i][j] = EMPTY;
                    M->celulas[i][j-1] = SURFACE;
                    continua = 1;
                }
            }
            /// colocar esse else if de volta (tirei pro FULL_SURFACE
//            else if( (i!=M->Nx-1 && M->celulas[i+1][j]==FULL) && (i!=0 && M->celulas[i-1][j]==FULL) && (j==M->Ny-1 || M->celulas[i][j+1]==FULL) && (M->pontosV[i][j].tipo!=BOUNDARY) ) {
//                M->celulas[i][j] = FULL;
//                M->celulas[i][j-1] = SURFACE;
//                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, FULL);
//                AdicionaCelulaSurface(ListaCelulas, i, j-1, PassoTemporal);
//                continua = 1;
//
//                M->celulasParede[i][j] = SURFACE;
//            }

            //Cantinho de baixo
            if( (i!=M->Nx-1 && M->celulas[i+1][j]==EMPTY) && (i!=0 && M->celulas[i-1][j]==EMPTY) && (j==0 || M->celulas[i][j-1]==EMPTY) && (M->pontosV[i][j+1].tipo!=BOUNDARY) ) {
                if( i==0 )
                    PrintDebug("AQUI 2\n");
                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, EMPTY);
                AdicionaCelulaSurface(M, ListaCelulas, i, j+1, PassoTemporal);
                M->celulas[i][j] = EMPTY;
                M->celulas[i][j+1] = SURFACE;
                continua = 1;
            }
            /// colocar esse else if de volta (tirei pro FULL_SURFACE
//            else if( (i!=M->Nx-1 && M->celulas[i+1][j]==FULL) && (i!=0 && M->celulas[i-1][j]==FULL) && (j==0 || M->celulas[i][j-1]==FULL) && (M->pontosV[i][j+1].tipo!=BOUNDARY) ) {
//                M->celulas[i][j] = FULL;
//                M->celulas[i][j+1] = SURFACE;
//                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, FULL);
//                AdicionaCelulaSurface(ListaCelulas, i, j+1, PassoTemporal);
//                continua = 1;
//
//                M->celulasParede[i][j] = SURFACE;
//            }

            //Cantinho da direita
            if( (i==M->Nx-1 || M->celulas[i+1][j]==EMPTY) && (j!=0 && M->celulas[i][j-1]==EMPTY) && (j!=M->Ny-1 && M->celulas[i][j+1]==EMPTY) && (M->pontosU[i][j].tipo!=BOUNDARY) ) {
                /// Caso abaixo desta celula tenha uma parede, nao eh preciso corrigir este caso (nao eh degenerado)
                if( M->pontosV[i][j].tipo!=BOUNDARY ) {
                    if( i==0 )
                        PrintDebug("AQUI 3\n");
                    RemoveCelulaSurface(ListaCelulas, celulaSurf, M, EMPTY);
                    AdicionaCelulaSurface(M, ListaCelulas, i-1, j, PassoTemporal);
                    M->celulas[i][j] = EMPTY;
                    M->celulas[i-1][j] = SURFACE;
                    continua = 1;
                }
            }
            else if( (i==M->Nx-1 || M->celulas[i+1][j]==FULL) && (j!=0 && M->celulas[i][j-1]==FULL) && (j!=M->Ny-1 && M->celulas[i][j+1]==FULL) && (M->pontosU[i][j].tipo!=BOUNDARY) ) {
                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, FULL);
                AdicionaCelulaSurface(M, ListaCelulas, i-1, j, PassoTemporal);
                M->celulas[i][j] = FULL;
                M->celulas[i-1][j] = SURFACE;
                continua = 1;
            }

            //Cantinho da esquerda
            if( (i==0 || M->celulas[i-1][j]==EMPTY) && (j!=0 && M->celulas[i][j-1]==EMPTY) && (j!=M->Ny-1 && M->celulas[i][j+1]==EMPTY) && (M->pontosU[i+1][j].tipo!=BOUNDARY) ) {

                /// Caso abaixo desta celula tenha uma parede, nao eh preciso corrigir este caso (nao eh degenerado)
                if( M->pontosU[i][j].tipo!=BOUNDARY && M->pontosV[i][j].tipo!=BOUNDARY ) {
                    RemoveCelulaSurface(ListaCelulas, celulaSurf, M, EMPTY);
                    AdicionaCelulaSurface(M, ListaCelulas, i+1, j, PassoTemporal);
                    M->celulas[i][j] = EMPTY;
                    M->celulas[i+1][j] = SURFACE;
                    continua = 1;
                }
            }
            else if( (i==0 || M->celulas[i-1][j]==FULL) && (j!=0 && M->celulas[i][j-1]==FULL) && (j!=M->Ny-1 && M->celulas[i][j+1]==FULL) && (M->pontosU[i+1][j].tipo!=BOUNDARY) ) {
                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, FULL);
                AdicionaCelulaSurface(M, ListaCelulas, i+1, j, PassoTemporal);
                M->celulas[i][j] = FULL;
                M->celulas[i+1][j] = SURFACE;
                continua = 1;
            }


            //Celula SURFACE com EMPTY na direita e esquerda
//            if( i!=0 && i!=M->Nx-1 && j!=0 && j!=M->Ny-1 && M->celulas[i-1][j]==EMPTY && M->celulas[i+1][j]==EMPTY && M->celulas[i][j-1]!=EMPTY && M->celulas[i][j+1]!=EMPTY ) {
//                M->celulas[i][j] = EMPTY;
//                M->celulas[i][j-1] = SURFACE;
//                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, EMPTY);
//                AdicionaCelulaSurface(ListaCelulas, i, j-1, PassoTemporal);
//            }
            //Celula SURFACE com EMPTY em cima e embaixo (PENSAR NUMA SOLUCAO DECENTE PRA ESSE CASO. NAO ESTA BOM ASSIM!!!)
//            else if( i!=0 && i!=M->Nx-1 && j!=0 && j!=M->Ny-1 && M->celulas[i-1][j]!=EMPTY && M->celulas[i+1][j]!=EMPTY && M->celulas[i][j-1]==EMPTY && M->celulas[i][j+1]==EMPTY ) {
//                M->celulas[i][j] = EMPTY;
//                M->celulas[i][j-1] = SURFACE;
//                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, EMPTY);
//                //AdicionaCelulaSurface(ListaCelulas, i, j-1, PassoTemporal);
//                continua = 1;
//            }

            celulaSurf = proxCelula;
        }


        ///Transformando em FULL as celulas que nao tem nenhuma empty em volta
        celulaSurf = ListaCelulas->prim;
        while( celulaSurf!=NULL ) {
            i = celulaSurf->i;
            j = celulaSurf->j;
            proxCelula = celulaSurf->prox;

            if( M->celulas[i][j]!=SURFACE ) {
                //DesenhaMalhaVTK(*M, 0);
                DeuProblema("\n 2: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! %d %d\n", i, j);
            }

            criterioDir = (i==M->Nx-1) || (M->celulas[i+1][j]!=EMPTY) || (M->pontosU[i+1][j].tipoBoundary==NOSLIP);
            criterioEsq = (i==0 || M->celulas[i-1][j]!=EMPTY) || (M->pontosU[i][j].tipoBoundary==NOSLIP);
            criterioCima = (j==M->Ny-1 || M->celulas[i][j+1]!=EMPTY) || (M->pontosV[i][j+1].tipoBoundary==NOSLIP);
            criterioBaixo = (j==0) || (M->celulas[i][j-1]!=EMPTY) || (M->pontosV[i][j].tipoBoundary==NOSLIP);

            if( criterioDir && criterioEsq && criterioCima && criterioBaixo ) {
                RemoveCelulaSurface(ListaCelulas, celulaSurf, M, FULL);
                M->celulas[i][j] = FULL;
                continua = 1;

                // Testando tensao superficial na celula full
                M->celulasParede[i][j] = SURFACE;
            }

            celulaSurf = proxCelula;
        }

    }
}

int EncontraGotasPequenas_Surface(MALHA *M, int i, int j, char **CelulasVisitadas)
{
    if( (i<0) || (j<0) || (i>=M->Nx) || (j>=M->Ny) )
        return 0;

    if( CelulasVisitadas[i][j] )
        return 0;

    CelulasVisitadas[i][j] = 1;

    if( M->celulas[i][j]==EMPTY )
        return 0;

    if( M->celulas[i][j]==SURFACE ) {
        return EncontraGotasPequenas_Surface(M, i+1, j, CelulasVisitadas) +
                EncontraGotasPequenas_Surface(M, i-1, j, CelulasVisitadas) +
                EncontraGotasPequenas_Surface(M, i, j+1, CelulasVisitadas) +
                EncontraGotasPequenas_Surface(M, i, j-1, CelulasVisitadas);
    }

    if( M->celulas[i][j]==FULL ) {
        return 1.0 +
                EncontraGotasPequenas_Surface(M, i+1, j, CelulasVisitadas) +
                EncontraGotasPequenas_Surface(M, i-1, j, CelulasVisitadas) +
                EncontraGotasPequenas_Surface(M, i, j+1, CelulasVisitadas) +
                EncontraGotasPequenas_Surface(M, i, j-1, CelulasVisitadas);
    }

    DeuProblema("RemoveGotasPequenas: nao era pra chegar aqui\n");
    return 0;
}

void RemoveGotaPequena_Surface(MALHA *M, int i, int j, char **CelulasVisitadas, LISTA_CELULAS_SURFACE *ListaSurf)
{
    if( (i<0) || (j<0) || (i>=M->Nx) || (j>=M->Ny) )
        return;

    if( CelulasVisitadas[i][j]!=1 )
        return;

    CelulasVisitadas[i][j] = 2;

    if( M->celulas[i][j]==EMPTY )
        return;

    if( M->celulas[i][j]==FULL )
        DeuProblema("RemoveGotaPequena: nao era pra ter uma full aqui...\n\n");


    if( M->celulas[i][j]==SURFACE ) {
        CELULA_SURFACE *celula_surf = EncontraCelulaSurface(ListaSurf, i, j);
        RemoveCelulaSurface(ListaSurf, celula_surf, M, -1);
        M->celulas[i][j] = EMPTY;

        RemoveGotaPequena_Surface(M, i+1, j, CelulasVisitadas, ListaSurf);
        RemoveGotaPequena_Surface(M, i-1, j, CelulasVisitadas, ListaSurf);
        RemoveGotaPequena_Surface(M, i, j+1, CelulasVisitadas, ListaSurf);
        RemoveGotaPequena_Surface(M, i, j-1, CelulasVisitadas, ListaSurf);
        return;
    }

    DeuProblema("RemoveGotasPequenas: nao era pra chegar aqui\n");
    return;
}

PetscErrorCode CondicaoContornoTangencial(int n, MALHA *M, LISTA_CELULAS_SURFACE listaCelulas, double **U, double **V, double **UNovo, double **VNovo, double **PNovo, double **Txx, double **Txy, double **Tyy, KSP *Ksp, Mat *A, Vec *Vet, Vec *Sol)
{
    int i, j, indice;
    int qtdIncognitas, qtdEquacoes, tamanhoBuffer;
    double valor;
    CELULA_SURFACE *celulaSurf;
    PetscErrorCode ierr;

    /* UMA BREVE EXPLICACAO!!!!

    As matrizes matIndicesU e matIndicesV possuem os indices de cada incognita em um dado sistema linear.
    No caso dessa funcao, elas irao conter o indice de cada incognita no sistema linear que sera construido pra superficie livre aqui.

    Primeiramente eu inicializo todos os matIndicesU/V com -1.
    Depois a funcao NovaIncognitaSuperficie vai criando as incognitas e preenchendo esses matIndices.

    A matriz pontosContinuidade contem uma flag para cada centro de celula. Ela indica em quais centros de celula nos ja aplicamos uma eq. de continuidade.
    Inicializa tudo com ZERO. Depois a funcao NovaEquacaoSuperficie vai preenchendo com UM conforme aplica as equacoes.

    A matriz pontosSupLivre contem uma flag para cada vertice da malha.
    Eh a mesma coisa da matriz pontosContinuidade, mas ela indica em quais vertices ja foi aplicado uma equacao tangencial de sup livre

    */


    /// Se nao tiver nenhuma celula SURFACE, nao tem nada pra fazer
    if( listaCelulas.prim==NULL )
        return 0;

    /// Inicializa as matrizes do sistema linear do contorno
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            M->pontosContinuidade[i][j] = 0;
            M->matIndicesP[i][j] = -1;
        }
    }
    for( i=0; i<=M->Nx; i++ ) {
        for( j=0; j<=M->Ny; j++ ) {
            M->pontosSupLivre[i][j] = 0;
        }
    }
    for( i=0; i<=M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            M->matIndicesU[i][j] = -1;
        }
    }
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<=M->Ny; j++ ) {
            M->matIndicesV[i][j] = -1;
        }
    }

    /// Criando (e ordenando) as incognitas do sistema linear da superficie livre
    qtdIncognitas = 0;
    tamanhoBuffer = 0;
    M->indicesInterface = NULL;
    for( celulaSurf=listaCelulas.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox ) {
        i = celulaSurf->i;
        j = celulaSurf->j;


        if( M->celulas[i][j]!=SURFACE ) {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 4: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! %d %d\n", i, j);
        }

        //Superficie vertical: empty na direita
        if( LADO_VERTICAL_DIREITA(i, j) ) {
            NovaIncognitaSuperficie(M, i+1, j, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        //Superficie vertical: empty na esquerda
        else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
            NovaIncognitaSuperficie(M, i, j, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i-1, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i-1, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        //Superficie horizontal: empty em cima
        else if( LADO_HORIZONTAL_CIMA(i, j) ) {
            NovaIncognitaSuperficie(M, i, j+1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j+1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        //Superficie horizontal: empty embaixo
        else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
            NovaIncognitaSuperficie(M, i, j-1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j-1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        //Superficie quina: empty em cima e na direita
        else if( QUINA_CIMA_DIREITA(i, j) ) {
            NovaIncognitaSuperficie(M, i+1, j, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j+1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);

            if( (i!=0) && M->celulas[i-1][j]==SURFACE && QUINA_CIMA_ESQUERDA(i-1, j) )
                NovaIncognitaSuperficie(M, i, j+1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);

            if( (j!=0) && M->celulas[i][j-1]==SURFACE && QUINA_BAIXO_DIREITA(i, j-1) )
                NovaIncognitaSuperficie(M, i+1, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        //Superficie quina: empty em cima e na esquerda
        else if( QUINA_CIMA_ESQUERDA(i, j) ) {
            NovaIncognitaSuperficie(M, i, j, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j+1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i-1, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);

            if( (i!=M->Nx-1) && M->celulas[i+1][j]==SURFACE && QUINA_CIMA_DIREITA(i+1, j) )
                NovaIncognitaSuperficie(M, i+1, j+1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);

            if( (j!=0) && M->celulas[i][j-1]==SURFACE && QUINA_BAIXO_ESQUERDA(i, j-1) )
                NovaIncognitaSuperficie(M, i-1, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        //Superficie quina: empty em embaixo e na direita
        else if( QUINA_BAIXO_DIREITA(i, j) ) {
            NovaIncognitaSuperficie(M, i+1, j, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j-1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);

            if( (i!=0) && M->celulas[i-1][j]==SURFACE && QUINA_BAIXO_ESQUERDA(i-1, j) ) {
                NovaIncognitaSuperficie(M, i, j-1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            }

            if( (j!=M->Ny-1) && M->celulas[i][j+1]==SURFACE && QUINA_CIMA_DIREITA(i, j+1) )
                NovaIncognitaSuperficie(M, i+1, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        //Superficie quina: empty em embaixo e na esquerda
        else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
            NovaIncognitaSuperficie(M, i, j, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j-1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i-1, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);

            if( (i!=M->Nx-1) && M->celulas[i+1][j]==SURFACE && QUINA_BAIXO_DIREITA(i+1, j) ) {
                NovaIncognitaSuperficie(M, i+1, j-1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            }

            if( (j!=M->Ny-1) && M->celulas[i][j+1]==SURFACE && QUINA_CIMA_ESQUERDA(i, j+1) )
                NovaIncognitaSuperficie(M, i-1, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        else if( DEGENERADO_CIMA_BAIXO(i, j) ) {
            NovaIncognitaSuperficie(M, i, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j+1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j+1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i, j-1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j-1, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        else if( DEGENERADO_ESQ_DIR(i, j) ) {
            NovaIncognitaSuperficie(M, i, j, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j, 'u', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i+1, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i-1, j, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
            NovaIncognitaSuperficie(M, i-1, j+1, 'v', &qtdIncognitas, 0, &tamanhoBuffer);
        }
        else {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n\n Tangencial 1: SUP LIVRE: CASO INESPERADO %d %d \n\n", i, j);
        }
    }

    /// ==== Inicializando as estruturas do sistema linear
    ierr = VecCreate(PETSC_COMM_WORLD, Vet); CHKERRQ(ierr);
	ierr = VecSetSizes(*Vet, PETSC_DECIDE, qtdIncognitas); CHKERRQ(ierr);
    ierr = VecSetFromOptions(*Vet); CHKERRQ(ierr);
    ierr = VecDuplicate(*Vet, Sol); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);//Cria a estrutura
    ierr = MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, qtdIncognitas, qtdIncognitas); CHKERRQ(ierr); //Define as dimensoes
    ierr = MatSetFromOptions(*A); CHKERRQ(ierr); //Finaliza a criacao da matriz
    ierr = MatSetUp(*A); CHKERRQ(ierr);
    //InicializaSolverKSP(Ksp, *A);

    //Criando o sistema linear (matriz e vetor independente)
    qtdEquacoes = 0;

    /// ==== Primeira passagem: aplica continuidade nas S e super_livre nas quinas
    for( celulaSurf=listaCelulas.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox ) {
        i = celulaSurf->i;
        j = celulaSurf->j;

        if( M->celulas[i][j]!=SURFACE ) {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 5: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! \n");
        }

        //Superficie vertical: empty na direita
        if( LADO_VERTICAL_DIREITA(i, j) ) {
            if( M->pontosU[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesU[i+1][j]) )
                M->pontosU[i+1][j].tipo = SURFACE;
        }
        //Superficie vertical: empty na esquerda
        else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
            if( M->pontosU[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesU[i][j]) )
                M->pontosU[i][j].tipo = SURFACE;
        }
        //Superficie horizontal: empty em cima
        else if( LADO_HORIZONTAL_CIMA(i, j) ) {
            if( M->pontosV[i][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesV[i][j+1]) )
                M->pontosV[i][j+1].tipo = SURFACE;
        }
        //Superficie horizontal: empty embaixo
        else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
            if( M->pontosV[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesV[i][j]) )
                M->pontosV[i][j].tipo = SURFACE;
        }
        //Superficie quina: empty em cima e na direita
        else if( QUINA_CIMA_DIREITA(i, j) ) {
            if( M->pontosU[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesU[i+1][j]) )
                M->pontosU[i+1][j].tipo = SURFACE;

            if( M->pontosV[i][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 1, DIR_CIMA, &qtdEquacoes, M->matIndicesV[i][j+1]) )
                M->pontosV[i][j+1].tipo = SURFACE;
        }
        //Superficie quina: empty em cima e na esquerda
        else if( QUINA_CIMA_ESQUERDA(i, j) ) {
            if( M->pontosU[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesU[i][j]) )
                M->pontosU[i][j].tipo = SURFACE;

            if( M->pontosV[i][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 1, ESQ_CIMA, &qtdEquacoes, M->matIndicesV[i][j+1]) )
                M->pontosV[i][j+1].tipo = SURFACE;
        }
        //Superficie quina: empty embaixo e na direita
        else if( QUINA_BAIXO_DIREITA(i, j) ) {
            if( M->pontosU[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesU[i+1][j]) ) {
                M->pontosU[i+1][j].tipo = SURFACE;
            }

            if( M->pontosV[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 1, DIR_BAIXO, &qtdEquacoes, M->matIndicesV[i][j]) )
                M->pontosV[i][j].tipo = SURFACE;


        }
        //Superficie quina: empty embaixo e na esquerda
        else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
            if( M->pontosU[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesU[i][j]) )
                M->pontosU[i][j].tipo = SURFACE;


            if( M->pontosV[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 1, ESQ_BAIXO, &qtdEquacoes, M->matIndicesV[i][j]) )
                M->pontosV[i][j].tipo = SURFACE;
        }
        else if( DEGENERADO_CIMA_BAIXO(i, j) ) {
            if( M->pontosV[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesV[i][j]) )
                M->pontosV[i][j].tipo = SURFACE;
            if( M->pontosV[i][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j+1, 3, 'd', &qtdEquacoes, M->matIndicesV[i][j+1]) )
                M->pontosV[i][j+1].tipo = SURFACE;
        }
        else if( DEGENERADO_ESQ_DIR(i, j) ) {
//            if( M->pontosU[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 0, 0, &qtdEquacoes, M->matIndicesU[i][j]) )
//                M->pontosU[i][j].tipo = SURFACE;
//            if( M->pontosU[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, 4, 'b', &qtdEquacoes, M->matIndicesU[i+1][j]) )
//                M->pontosU[i+1][j].tipo = SURFACE;
        }
        else DeuProblema("\n\n Tangencial 2: SUP LIVRE: CASO INESPERADO %d %d \n\n", i, j);
    }


    /// === Segunda passagem aplicando a eq de superficie livre nos vertices
    for( celulaSurf=listaCelulas.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox ) {
        i = celulaSurf->i;
        j = celulaSurf->j;

        if( M->celulas[i][j]!=SURFACE ) {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 6: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! \n");
        }


        //Superficie vertical: empty na direita
        if( LADO_VERTICAL_DIREITA(i, j) ) {
            if( (i!=M->Nx-1) && M->pontosV[i+1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j+1, 2, 'd', &qtdEquacoes, M->matIndicesV[i+1][j+1]) ) {
                M->pontosV[i+1][j+1].tipo = SURFACE;
            }

            if( (i!=M->Nx-1) && M->pontosV[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, 2, 'd', &qtdEquacoes, M->matIndicesV[i+1][j]) ) {
                M->pontosV[i+1][j].tipo = SURFACE;
            }
        }
        //Superficie vertical: empty na esquerda
        else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
            if( (i!=0) && M->pontosV[i-1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j+1, 2, 'e', &qtdEquacoes, M->matIndicesV[i-1][j+1]) ) {
                M->pontosV[i-1][j+1].tipo = SURFACE;
            }

            if( (i!=0) && M->pontosV[i-1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 2, 'e', &qtdEquacoes, M->matIndicesV[i-1][j]) ) {
                M->pontosV[i-1][j].tipo = SURFACE;
            }
        }
        //Superficie horizontal: empty em cima
        else if( LADO_HORIZONTAL_CIMA(i, j) ) {
            if( (j!=M->Ny-1) && M->pontosU[i][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j+1, 2, 'c', &qtdEquacoes, M->matIndicesU[i][j+1]) ) {
                M->pontosU[i][j+1].tipo = SURFACE;
            }

            if( (j!=M->Ny-1) && M->pontosU[i+1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j+1, 2, 'c', &qtdEquacoes, M->matIndicesU[i+1][j+1]) ) {
                M->pontosU[i+1][j+1].tipo = SURFACE;
            }
        }
        //Superficie horizontal: empty embaixo
        else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
            if( (j!=0) && M->pontosU[i][j-1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 2, 'b', &qtdEquacoes, M->matIndicesU[i][j-1]) ) {
                M->pontosU[i][j-1].tipo = SURFACE;
            }

            if( (j!=0) && M->pontosU[i+1][j-1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, 2, 'b', &qtdEquacoes, M->matIndicesU[i+1][j-1]) ) {
                M->pontosU[i+1][j-1].tipo = SURFACE;
            }
        }
        //Extrapolacao quina: empty em cima e na direita
        else if( QUINA_CIMA_DIREITA(i, j) ) {
            if( (j!=M->Ny-1) && M->pontosU[i+1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j+1, 4, 'b', &qtdEquacoes, M->matIndicesU[i+1][j+1]) ) {
                M->pontosU[i+1][j+1].tipo = SURFACE;
            }

            if( (i!=M->Nx-1) && M->pontosV[i+1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j+1, 3, 'e', &qtdEquacoes, M->matIndicesV[i+1][j+1]) ) {
                M->pontosV[i+1][j+1].tipo = SURFACE;
            }

            if( (i!=0) && (j!=M->Ny-1) && M->celulas[i-1][j]==SURFACE && QUINA_CIMA_ESQUERDA(i-1, j) ) {
                if( M->pontosU[i][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j+1, 2, 'c', &qtdEquacoes, M->matIndicesU[i][j+1]) )
                    M->pontosU[i][j+1].tipo = SURFACE;
            }

            if( (i!=M->Nx-1) && (j!=0) && (M->celulas[i][j-1]==SURFACE && QUINA_BAIXO_DIREITA(i, j-1)) ) {
                if( M->pontosV[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, 2, 'd', &qtdEquacoes, M->matIndicesV[i+1][j]) )
                    M->pontosV[i+1][j].tipo = SURFACE;
            }
        }
        //Superficie quina: empty em cima e na esquerda
        else if( QUINA_CIMA_ESQUERDA(i, j) ) {
            if( (j!=M->Ny-1) && M->pontosU[i][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j+1, 4, 'b', &qtdEquacoes, M->matIndicesU[i][j+1]) ) {
                M->pontosU[i][j+1].tipo = SURFACE;
            }

            if( (i!=0) && M->pontosV[i-1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i-1, j+1, 3, 'd', &qtdEquacoes, M->matIndicesV[i-1][j+1]) ) {
                M->pontosV[i-1][j+1].tipo = SURFACE;
            }

            if( (i!=M->Nx-1) && M->celulas[i+1][j]==SURFACE && QUINA_CIMA_DIREITA(i+1, j) ) {
                if( M->pontosU[i+1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j+1, 2, 'c', &qtdEquacoes, M->matIndicesU[i+1][j+1]) )
                    M->pontosU[i+1][j+1].tipo = SURFACE;
            }

            if( (i!=0) && (j!=0) && (M->celulas[i][j-1]==SURFACE && QUINA_BAIXO_ESQUERDA(i, j-1)) ) {
                if( M->pontosV[i-1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 2, 'e', &qtdEquacoes, M->matIndicesV[i-1][j]) )
                    M->pontosV[i-1][j].tipo = SURFACE;
            }
        }
        //Superficie quina: empty embaixo e na direita
        else if( QUINA_BAIXO_DIREITA(i, j) ) {
            if( (j!=0) && M->pontosU[i+1][j-1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j-1, 4, 'c', &qtdEquacoes, M->matIndicesU[i+1][j-1]) ) {
                M->pontosU[i+1][j-1].tipo = SURFACE;
            }

            if( (i!=M->Nx-1) && M->pontosV[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, 3, 'e', &qtdEquacoes, M->matIndicesV[i+1][j]) ) {
                M->pontosV[i+1][j].tipo = SURFACE;
            }

            if( (i!=0) && (j!=0) && M->celulas[i-1][j]==SURFACE && QUINA_BAIXO_ESQUERDA(i-1, j) ) {
                if( M->pontosU[i][j-1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, 2, 'b', &qtdEquacoes, M->matIndicesU[i][j-1]) ) {
                    M->pontosU[i][j-1].tipo = SURFACE;
                }
            }

            if( (i!=M->Nx-1) && ((j!=M->Ny-1) && M->celulas[i][j+1]==SURFACE && QUINA_CIMA_DIREITA(i, j+1)) ) {
                if( M->pontosV[i+1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j+1, 2, 'd', &qtdEquacoes, M->matIndicesV[i+1][j+1]) )
                    M->pontosV[i+1][j+1].tipo = SURFACE;
            }
        }
        //Superficie quina: empty embaixo e na esquerda
        else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
            if( (j!=0) && M->pontosU[i][j-1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j-1, 4, 'c', &qtdEquacoes, M->matIndicesU[i][j-1]) ) {
                M->pontosU[i][j-1].tipo = SURFACE;
            }

            if( (i!=0) && M->pontosV[i-1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i-1, j, 3, 'd', &qtdEquacoes, M->matIndicesV[i-1][j]) ) {
                M->pontosV[i-1][j].tipo = SURFACE;
            }


            if( (i!=M->Nx-1) && M->celulas[i+1][j]==SURFACE && QUINA_BAIXO_DIREITA(i+1, j) ) {
                if( M->pontosU[i+1][j-1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, 2, 'b', &qtdEquacoes, M->matIndicesU[i+1][j-1]) ) {
                    M->pontosU[i+1][j-1].tipo = SURFACE;
                }
            }

            if( (i!=0) && ((j!=M->Ny-1) && M->celulas[i][j+1]==SURFACE && QUINA_CIMA_ESQUERDA(i, j+1)) ) {
                if( M->pontosV[i-1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j+1, 2, 'e', &qtdEquacoes, M->matIndicesV[i-1][j+1]) )
                    M->pontosV[i-1][j+1].tipo = SURFACE;
            }
        }
        else if( DEGENERADO_CIMA_BAIXO(i, j) ) {
            if( M->pontosU[i][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j+1, 2, 'c', &qtdEquacoes, M->matIndicesU[i][j+1]) )
                M->pontosU[i][j+1].tipo = SURFACE;
            if( M->pontosU[i+1][j-1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, 2, 'b', &qtdEquacoes, M->matIndicesU[i+1][j-1]) )
                M->pontosU[i+1][j-1].tipo = SURFACE;

            //Vai precisar adicionar mais algumas incognitas aqui ainda!!!
        }
        else if( DEGENERADO_ESQ_DIR(i, j) ) {
            //Vai precisar adicionar alguns casos aqui ainda!!!
        }
        else DeuProblema("\n\n Tangencial 3: SUP LIVRE: CASO INESPERADO %d %d \n\n", i, j);
    }




    /// Fazendo uma ultima passagem apenas para uns casos degenerados horriveis por ultimo...
    for( celulaSurf=listaCelulas.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox ) {

        i = celulaSurf->i;
        j = celulaSurf->j;

        if( M->celulas[i][j]!=SURFACE ) {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 5: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! \n");
        }

        if( DEGENERADO_ESQ_DIR(i, j) ) {
            if( M->pontosU[i][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i, j, -1, 0, &qtdEquacoes, M->matIndicesU[i][j]) )
                M->pontosU[i][j].tipo = SURFACE;
            if( M->pontosU[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, -1, 0, &qtdEquacoes, M->matIndicesU[i+1][j]) )
                M->pontosU[i+1][j].tipo = SURFACE;



            if( M->pontosV[i+1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j, -2, 0, &qtdEquacoes, M->matIndicesV[i+1][j]) )
                M->pontosV[i+1][j].tipo = SURFACE;
            if( M->pontosV[i+1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i+1, j+1, -2, 0, &qtdEquacoes, M->matIndicesV[i+1][j+1]) )
                M->pontosV[i+1][j+1].tipo = SURFACE;

            if( M->pontosV[i-1][j].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i-1, j, -2, 0, &qtdEquacoes, M->matIndicesV[i-1][j]) )
                M->pontosV[i-1][j].tipo = SURFACE;
            if( M->pontosV[i-1][j+1].tipo==INCOGNITA && NovaEquacaoSuperficie(*M, *A, *Vet, U, V, Txx, Txy, Tyy, i-1, j+1, -2, 0, &qtdEquacoes, M->matIndicesV[i-1][j+1]) )
                M->pontosV[i-1][j+1].tipo = SURFACE;
        }

    }

    ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(*Vet); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*Vet); CHKERRQ(ierr);


    if( qtdEquacoes!=qtdIncognitas ) {
        DesenhaMalhaVTK(*M, 0);
//        ImprimeInterfaceVTK(*M, 1000000);
        ImprimeMatrizPraMatlab(*A, "debugs/matrizcanto.m");
        ImprimeVetorPraMatlab(*Vet, "debugs/vetorcanto.m");
        DeuProblema("\n\n PROBLEMA SUPERFICIE: QTD EQUACOES  DIFERENTE DE QTD INCOGNITAS \n Inc: %d  Eq: %d \n\n", qtdIncognitas, qtdEquacoes);
    }


    //Resolvendo o sistema linear
    InicializaSolverKSP(Ksp, *A);
    ierr = KSPSolve(*Ksp, *Vet, *Sol); CHKERRQ(ierr);

    //Jogando a solucao do sistema pra matriz de valores U e V
    for( indice=0; indice<qtdIncognitas; indice++ ) {
        if( M->indicesInterface[indice].descarta )
            continue;

        VecGetValues(*Sol, 1, &indice, &valor);
        if( M->indicesInterface[indice].incognita=='u' ) {
            UNovo[M->indicesInterface[indice].i][M->indicesInterface[indice].j] = valor;
            U[M->indicesInterface[indice].i][M->indicesInterface[indice].j] = valor;
        }
        else if( M->indicesInterface[indice].incognita=='v' ) {
            VNovo[M->indicesInterface[indice].i][M->indicesInterface[indice].j] = valor;
            V[M->indicesInterface[indice].i][M->indicesInterface[indice].j] = valor;
        }
        else
            printf("\nPROBLEMa final do sistema superficie!!!!!\n");
    }


    //Destruindo as estruturas usadas no sistema linear
    ierr = VecDestroy(Vet); CHKERRQ(ierr);
    ierr = VecDestroy(Sol); CHKERRQ(ierr);
    ierr = MatDestroy(A); CHKERRQ(ierr);
    ierr = KSPDestroy(Ksp); CHKERRQ(ierr);
    free(M->indicesInterface);

    return 0;
}

void NovaIncognitaSuperficie(MALHA *M, int i, int j, char Incognita, int *QtdIncognitas, char Descarta, int *TamanhoBuffer)
{
    POSICAO *temp;

    if( Incognita=='u' ) {
        if( i<0 || i>M->Nx || j<0 || j>=M->Ny )
            return;

        if( M->pontosU[i][j].tipo!=EMPTY && M->pontosU[i][j].tipo!=SURFACE ) {

            //Soh pra ter certeza que vai descartar se descarta=1
            int ind;
            for( ind=0; ind<*QtdIncognitas; ind++ ) {
                if( M->indicesInterface[ind].incognita=='u' && M->indicesInterface[ind].i==i && M->indicesInterface[ind].j==j && Descarta==1 )
                    M->indicesInterface[ind].descarta = 1;
            }

            return;
        }
        M->pontosU[i][j].tipo = INCOGNITA; //Incognita
        M->matIndicesU[i][j] = *QtdIncognitas;
    }
    else if( Incognita=='v' ) {
        if( i<0 || i>=M->Nx || j<0 || j>M->Ny )
            return;

        if( M->pontosV[i][j].tipo!=EMPTY && M->pontosV[i][j].tipo!=SURFACE ) {
            int ind;
            for( ind=0; ind<*QtdIncognitas; ind++ ) {
                if( M->indicesInterface[ind].incognita=='v' && M->indicesInterface[ind].i==i && M->indicesInterface[ind].j==j && Descarta==1 )
                    M->indicesInterface[ind].descarta = 1;
            }
            return;
        }
        M->pontosV[i][j].tipo = INCOGNITA; //Incognita
        M->matIndicesV[i][j] = *QtdIncognitas;
    }
    else
        printf("\n\n PROBLEMA: NOVA INCOGNITA \n\n");

    //Aloca mais memoria
    if( *TamanhoBuffer==0 ) {
        temp = (POSICAO *)realloc(M->indicesInterface, ((*QtdIncognitas)+10)*sizeof(POSICAO));
        M->indicesInterface = temp;
        *TamanhoBuffer = 10;
    }

    M->indicesInterface[*QtdIncognitas].i = i;
    M->indicesInterface[*QtdIncognitas].j = j;
    M->indicesInterface[*QtdIncognitas].incognita = Incognita;
    M->indicesInterface[*QtdIncognitas].descarta = Descarta;
    (*TamanhoBuffer)--;
    (*QtdIncognitas)++;


    return;
}

int NovaEquacaoSuperficie(MALHA M, Mat A, Vec Vet, double **U, double **V, double **Txx, double **Txy, double **Tyy, int i, int j, char TipoEquacao, char DirecaoExtrapolacao, int *QtdEquacoes, int LinhaMatriz)
{
    static double valor, ladoDireito;
    static int indice1, indice2, coluna;
    static double invDx, invDy, coefUbaixo, coefUcima, coefVesq, coefVdir;
    static double h1, h2;
    //static double nx, ny;
 
    double rImenosMeio, rImaisMeio, rI, cilindrico;
    cilindrico = (M.tipoCoord==AXI_CILINDRICO) ? 1.0 : 0.0;

    ladoDireito = 0.0;


    if( TipoEquacao==-1 ) { //Colocando ZERO nessa incognita
        coluna = LinhaMatriz;
        valor = 1.0;
        MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        ladoDireito = 0.0;
        VecSetValues(Vet, 1, &LinhaMatriz, &ladoDireito, INSERT_VALUES);
        (*QtdEquacoes)++;
    }
    else if( TipoEquacao==-2 ) { //Fazendo uma media vertical em V, apenas para um caso degenerado

        valor = -2.0;
        if( M.pontosV[i][j].tipo==SURFACE || M.pontosV[i][j].tipo==INCOGNITA )
            MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesV[i][j]), &valor, INSERT_VALUES);
        else
            ladoDireito += -valor*V[i][j];

        valor = 1.0;
        if( M.pontosV[i][j-1].tipo==SURFACE || M.pontosV[i][j-1].tipo==INCOGNITA )
            MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesV[i][j-1]), &valor, INSERT_VALUES);
        else
            ladoDireito += -valor*V[i][j-1];

        valor = 1.0;
        if( M.pontosV[i][j+1].tipo==SURFACE || M.pontosV[i][j+1].tipo==INCOGNITA )
            MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesV[i][j+1]), &valor, INSERT_VALUES);
        else
            ladoDireito += -valor*V[i][j+1];

        VecSetValues(Vet, 1, &LinhaMatriz, &ladoDireito, INSERT_VALUES);
        (*QtdEquacoes)++;

    }
    else if( TipoEquacao==0 ) { //Continuidade no centro da celula
        if( M.celulas[i][j]==FULL || M.pontosContinuidade[i][j] )
            return 0;
        M.pontosContinuidade[i][j] = 1;


        invDx = 1.0/M.dx[i];
        invDy = 1.0/M.dy[j];

	//Valores da coordenada r usados no caso axisimetrico
        if( M.tipoCoord==AXI_CILINDRICO ) {
            rI = M.r[i];
            rImenosMeio = M.r[i] - 0.5*M.dx[i];
            rImaisMeio = M.r[i] + 0.5*M.dx[i];
        }
        else
            rImenosMeio = rImaisMeio = rI = 1.0;


        /// === Criando esta equacao na matriz
        //Valor U da direita
        valor = (rImaisMeio/rI)*invDx;
        if( M.pontosU[i+1][j].tipo==SURFACE || M.pontosU[i+1][j].tipo==INCOGNITA ) {
            coluna = M.matIndicesU[i+1][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else
            ladoDireito += -valor*U[i+1][j];


        //Valor U da esquerda
        valor = -(rImenosMeio/rI)*invDx;
        if( M.pontosU[i][j].tipo==SURFACE || M.pontosU[i][j].tipo==INCOGNITA ) {
            coluna = M.matIndicesU[i][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else
            ladoDireito += -valor*U[i][j];

        //Valor V de cima
        valor = invDy;
        if( M.pontosV[i][j+1].tipo==SURFACE || M.pontosV[i][j+1].tipo==INCOGNITA ) {
            coluna = M.matIndicesV[i][j+1];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else
            ladoDireito += -valor*V[i][j+1];

        //Valor V de baixo
        valor = -invDy;
        if( M.pontosV[i][j].tipo==SURFACE || M.pontosV[i][j].tipo==INCOGNITA ) {
            coluna = M.matIndicesV[i][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else
            ladoDireito += -valor*V[i][j];

        VecSetValues(Vet, 1, &LinhaMatriz, &ladoDireito, INSERT_VALUES);
        (*QtdEquacoes)++;
    }
    else if( TipoEquacao==1 ) { //Condicao de superficie livre no centro da celula (com normal 45 graus)
        //if( M.pontosSupLivre[i][j] )
          //  return 0;
        //M.pontosSupLivre[i][j] = 1;

//        if( DirecaoExtrapolacao==DIR_CIMA ) {
//            nx = sqrt(2.0)/2.0;
//            ny = sqrt(2.0)/2.0;
//        }
//        else if( DirecaoExtrapolacao==ESQ_CIMA ) {
//            nx = -sqrt(2.0)/2.0;
//            ny = sqrt(2.0)/2.0;
//        }
//        else if( DirecaoExtrapolacao==DIR_BAIXO ) {
//            nx = sqrt(2.0)/2.0;
//            ny = -sqrt(2.0)/2.0;
//        }
//        else if( DirecaoExtrapolacao==ESQ_BAIXO ) {
//            nx = -sqrt(2.0)/2.0;
//            ny = -sqrt(2.0)/2.0;
//        }
//        else printf("\n\n PROBLEMA EQUACAO SUPERFICIE LIVRE \n\n");

        invDx = 1.0/M.dx[i];
        invDy = 1.0/M.dy[j];

        /// === Criando esta equacao na matriz
        //Valor U da direita
        valor = -invDx;
        if( M.pontosU[i+1][j].tipo==SURFACE || M.pontosU[i+1][j].tipo==INCOGNITA ) {
            coluna = M.matIndicesU[i+1][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else
            ladoDireito += -valor*U[i+1][j];

        //Valor U da esquerda
        valor = invDx;
        if( M.pontosU[i][j].tipo==SURFACE || M.pontosU[i][j].tipo==INCOGNITA ) {
            coluna = M.matIndicesU[i][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else
            ladoDireito += -valor*U[i][j];

        valor = invDy;
        if( M.pontosV[i][j+1].tipo==SURFACE || M.pontosV[i][j+1].tipo==INCOGNITA ) {
            coluna = M.matIndicesV[i][j+1];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else
            ladoDireito += -valor*V[i][j+1];

        valor = -invDy;
        if( M.pontosV[i][j].tipo==SURFACE || M.pontosV[i][j].tipo==INCOGNITA ) {
            coluna = M.matIndicesV[i][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else
            ladoDireito += -valor*V[i][j];

        if( !velocidadeNewtoniana ) {
            double beta = (M.tipo_modelo==MODELO_EVPT) ? M.eta_inf : M.beta;
            ladoDireito += (1.0/(2.0*beta))*( Txx[i][j] - Tyy[i][j] );
        }
        VecSetValues(Vet, 1, &LinhaMatriz, &ladoDireito, INSERT_VALUES);
        (*QtdEquacoes)++;
    }
    else if( TipoEquacao==2 ) { //Condicao de superficie livre no vertice (com normal vertical ou horizontal)
        if( M.pontosSupLivre[i][j] )
            return 0;
        M.pontosSupLivre[i][j] = 1;

//        if( DirecaoExtrapolacao=='c' ) {
//            nx = 0.0;
//            ny = 1.0;
//        }
//        else if( DirecaoExtrapolacao=='b' ) {
//            nx = 0.0;
//            ny = -1.0;
//        }
//        else if( DirecaoExtrapolacao=='d' ) {
//            nx = 1.0;
//            ny = 0.0;
//        }
//        else if( DirecaoExtrapolacao=='e' ) {
//            nx = -1.0;
//            ny = 0.0;
//        }
//        else printf("\n\n PROBLEMA EQUACAO SUPERFICIE LIVRE \n\n");


        h1 = (j==0) ? 0.5*M.dy[j] : 0.5*M.dy[j-1];
        h2 = (j==M.Ny) ? 0.5*M.dy[j-1] : 0.5*M.dy[j];
        coefUbaixo = ( (( h2*(h2-h1) )/( (h1+h2)*h1*h2 )) - ( h2/(h1*(h1+h2)) ) );
        coefUcima  = ( (( h1*(h2-h1) )/( (h1+h2)*h1*h2 )) + ( h1/(h2*(h1+h2)) ) );
        h1 = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        h2 = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];
        coefVesq = ( (( h2*(h2-h1) )/( (h1+h2)*h1*h2 )) - ( h2/(h1*(h1+h2)) ) );
        coefVdir = ( (( h1*(h2-h1) )/( (h1+h2)*h1*h2 )) + ( h1/(h2*(h1+h2)) ) );

        /// === Criando esta equacao na matriz
        valor = coefUcima;
        if( M.pontosU[i][j].tipo==SURFACE || M.pontosU[i][j].tipo==INCOGNITA ) {
            coluna = M.matIndicesU[i][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else if( M.pontosU[i][j].tipo==FULL || M.pontosU[i][j].tipo==BOUNDARY )
            ladoDireito += -valor*U[i][j];
        else DeuProblema("\n\n 1. PROBLEMA: NovaEquacaoSuperficie: %d %d\n", i, j);

        valor = coefUbaixo;
        if( (j!=0) && (M.pontosU[i][j-1].tipo==SURFACE || M.pontosU[i][j-1].tipo==INCOGNITA) ) {
            coluna = M.matIndicesU[i][j-1];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else if( (j!=0) &&  (M.pontosU[i][j-1].tipo==FULL || M.pontosU[i][j-1].tipo==BOUNDARY) )
            ladoDireito += -valor*U[i][j-1];
//        else if( (j==0)  )
//            ladoDireito += -valor*(-U[i][j]);
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("\n\n 2. PROBLEMA: NovaEquacaoSuperficie: %d %d\n", i, j);
        }

        valor = coefVdir;
        if( (i!=M.Nx) && (M.pontosV[i][j].tipo==SURFACE || M.pontosV[i][j].tipo==INCOGNITA) ) {
            coluna = M.matIndicesV[i][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else if( (i!=M.Nx) && (M.pontosV[i][j].tipo==FULL || M.pontosV[i][j].tipo==BOUNDARY) )
            ladoDireito += -valor*V[i][j];
        else if( M.pontosU[i][j].tipo==BOUNDARY || M.pontosU[i][j-1].tipo==BOUNDARY ) { //Ta no extremo da direita
            int tipoBoundary = (M.pontosU[i][j].tipo==BOUNDARY) ? M.pontosU[i][j].tipoBoundary : M.pontosU[i][j-1].tipoBoundary;
            if( tipoBoundary==NOSLIP )
                ladoDireito += -valor*(-V[i-1][j]);
            else if( tipoBoundary==INFLOW )
                ladoDireito += -valor*(-V[i-1][j]);
            else if( tipoBoundary==NEUMANN )
                ladoDireito += -valor*(V[i-1][j]);
            else DeuProblema("\n\n Problema NovaEquacaoSuperficie. Tipo do contorno. %d %d\n\n", i, j);
        }
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("\n\n 3. PROBLEMA: NovaEquacaoSuperficie: %d %d %d\n", i, j, M.pontosU[i][j].tipo);
        }

        valor = coefVesq;
        if( M.pontosV[i-1][j].tipo==SURFACE || M.pontosV[i-1][j].tipo==INCOGNITA ) {
            coluna = M.matIndicesV[i-1][j];
            MatSetValues(A, 1, &LinhaMatriz, 1, &coluna, &valor, INSERT_VALUES);
        }
        else if( M.pontosV[i-1][j].tipo==FULL || M.pontosV[i-1][j].tipo==BOUNDARY )
            ladoDireito += -valor*V[i-1][j];
        else DeuProblema("\n\n 4. PROBLEMA: NovaEquacaoSuperficie: %d %d\n", i, j);

        double txyCentroCima = 0.0, txyCentroBaixo = 0.0;
        double txyDirCima = 0.0, txyEsqCima = 0.0, txyDirBaixo = 0.0, txyEsqBaixo = 0.0;

        //Extrapolando os tensores que caem em celulas EMPTY
        if( (i==M.Nx) || (M.pontosP[i][j].tipo==EMPTY) ) {
            if( (i!=M.Nx) && M.pontosP[i-1][j].tipo!=EMPTY && M.pontosP[i][j-1].tipo!=EMPTY )
                txyDirCima = 0.5*( Txy[i-1][j] + Txy[i][j-1] );
            else if( M.pontosP[i-1][j].tipo!=EMPTY )
                txyDirCima = Txy[i-1][j];
            else if( (i!=M.Nx) && M.pontosP[i][j-1].tipo!=EMPTY )
                txyDirCima = Txy[i][j-1];
            else if( M.pontosP[i-1][j-1].tipo!=EMPTY ) //Extrapolando na diagonal em ultimo caso
                txyDirCima = Txy[i-1][j-1];
            else
                DeuProblema("NovaEquacaoSuperficie: Extrapolacao Txy 1. %d %d\n\n", i, j);
        }
        else
            txyDirCima = Txy[i][j];

        if( (i==M.Nx) || (M.pontosP[i][j-1].tipo==EMPTY) ) {
            if( (i!=M.Nx) && M.pontosP[i-1][j-1].tipo!=EMPTY && M.pontosP[i][j].tipo!=EMPTY )
                txyDirBaixo = 0.5*( Txy[i-1][j-1] + Txy[i][j] );
            else if( M.pontosP[i-1][j-1].tipo!=EMPTY )
                txyDirBaixo = Txy[i-1][j-1];
            else if( (i!=M.Nx) && M.pontosP[i][j].tipo!=EMPTY )
                txyDirBaixo = Txy[i][j];
            else if( M.pontosP[i-1][j].tipo!=EMPTY ) //Extrapolando na diagonal em ultimo caso
                txyDirBaixo = Txy[i-1][j];
            else DeuProblema("NovaEquacaoSuperficie: Extrapolacao Txy 2. %d %d\n\n", i, j);
        }
        else
            txyDirBaixo = Txy[i][j-1];

        if( (M.pontosP[i-1][j].tipo==EMPTY) ) {
            if( (i!=M.Nx) && M.pontosP[i][j].tipo!=EMPTY && M.pontosP[i-1][j-1].tipo!=EMPTY )
                txyEsqCima = 0.5*( Txy[i][j] + Txy[i-1][j-1] );
            else if( (i!=M.Nx) && M.pontosP[i][j].tipo!=EMPTY )
                txyEsqCima = Txy[i][j];
            else if( M.pontosP[i-1][j-1].tipo!=EMPTY )
                txyEsqCima = Txy[i-1][j-1];
            else if( M.pontosP[i][j-1].tipo!=EMPTY ) //Faz na diagonal em ultimo caso
                txyEsqCima = Txy[i][j-1];
            else
                DeuProblema("NovaEquacaoSuperficie: Extrapolacao Txy 3. %d %d\n\n", i, j);
        }
        else
            txyEsqCima = Txy[i-1][j];

        if( (M.pontosP[i-1][j-1].tipo==EMPTY) ) {
            if( (i!=M.Nx) && M.pontosP[i][j-1].tipo!=EMPTY && M.pontosP[i-1][j].tipo!=EMPTY )
                txyEsqBaixo = 0.5*( Txy[i][j-1] + Txy[i-1][j] );
            else if( (i!=M.Nx) && M.pontosP[i][j-1].tipo!=EMPTY )
                txyEsqBaixo = Txy[i][j-1];
            else if( M.pontosP[i-1][j].tipo!=EMPTY )
                txyEsqBaixo = Txy[i-1][j];
            else if( M.pontosP[i][j].tipo!=EMPTY ) //Extrapolando na diagonal em ultimo caso
                txyEsqBaixo = Txy[i][j];
            else
                DeuProblema("NovaEquacaoSuperficie: Extrapolacao Txy 4. %d %d\n\n", i, j);
        }
        else
            txyEsqBaixo = Txy[i-1][j-1];

        double dxEsq = (i==0) ? 0.5*M.dx[i] : 0.5*M.dx[i-1];
        double dxDir = (i==M.Nx) ? 0.5*M.dx[i-1] : 0.5*M.dx[i];
        txyCentroBaixo = Interpolacao(dxEsq, dxDir, txyEsqBaixo, txyDirBaixo);
        txyCentroCima = Interpolacao(dxEsq, dxDir, txyEsqCima, txyDirCima);

        if( !velocidadeNewtoniana ) {
            double beta = (M.tipo_modelo==MODELO_EVPT) ? M.eta_inf : M.beta;
            ladoDireito += (-1.0/beta)*Interpolacao(0.5*M.dy[j-1], 0.5*M.dy[j], txyCentroBaixo, txyCentroCima);
        }
        VecSetValues(Vet, 1, &LinhaMatriz, &ladoDireito, INSERT_VALUES);
        (*QtdEquacoes)++;
    }
    else if( TipoEquacao==3 ) { //Extrapolacao linear de V


        if( DirecaoExtrapolacao=='e' ) { //Esquerda
            //coord = M.x[i] + 0.5*M.dx[i];
            //coord1 = M.x[i-1] + 0.5*M.dx[i-1];
            indice1 = i-1;
            indice2 = i-2;
        }
        else if( DirecaoExtrapolacao=='d' ) { //Direita
            //coord = M.x[i] + 0.5*M.dx[i];
            //coord1 = M.x[i+1] + 0.5*M.dx[i+1];
            indice1 = i+1;
            indice2 = i+2;
        }
        else DeuProblema("\n\n EXTRAPOLACAO V: CASO INESPERADO \n\n");

        if( DirecaoExtrapolacao=='e' || DirecaoExtrapolacao=='d') {

            valor = -1.0;
            if( M.pontosV[i][j].tipo==SURFACE || M.pontosV[i][j].tipo==INCOGNITA )
                MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesV[i][j]), &valor, INSERT_VALUES);
            else
                ladoDireito += -valor*V[i][j];

            valor = 2.0;
            if( M.pontosV[indice1][j].tipo==SURFACE || M.pontosV[indice1][j].tipo==INCOGNITA )
                MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesV[indice1][j]), &valor, INSERT_VALUES);
            else
                ladoDireito += -valor*V[indice1][j];

            valor = -1.0;
            if( indice2<M.Nx &&  indice2>0 && (M.pontosV[indice2][j].tipo==SURFACE || M.pontosV[indice2][j].tipo==INCOGNITA) )
                MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesV[indice2][j]), &valor, INSERT_VALUES);
            else if( indice2<M.Nx && indice2>0 )
                ladoDireito += -valor*V[indice2][j];

        }

        VecSetValues(Vet, 1, &LinhaMatriz, &ladoDireito, INSERT_VALUES);
        (*QtdEquacoes)++;
    }
    else if( TipoEquacao==4 ) { //Extrapolacao linear de U

        if( DirecaoExtrapolacao=='b' ) { //Baixo
            //coord = M.y[j] + 0.5*M.dy[j];
            //coord1 = M.y[j-1] + 0.5*M.dy[j-1];
            indice1 = j-1;
            indice2 = j-2;
        }
        else if( DirecaoExtrapolacao=='c' ) { //cima
            //coord = M.y[j] + 0.5*M.dy[j];
            //coord1 = M.y[j+1] + 0.5*M.dy[j+1];
            indice1 = j+1;
            indice2 = j+2;
        }
        else DeuProblema("\n\nNovaEquacaoSuperficie: Tipo=4. Direcao nao prevista na extrapolacao.\n\n");

        if( DirecaoExtrapolacao=='b' || DirecaoExtrapolacao=='c' ) {
            valor = -1.0;
            if( M.pontosU[i][j].tipo==SURFACE || M.pontosU[i][j].tipo==INCOGNITA )
                MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesU[i][j]), &valor, INSERT_VALUES);
            else
                ladoDireito += -valor*U[i][j];

            valor = 2.0;
            if( M.pontosU[i][indice1].tipo==SURFACE || M.pontosU[i][indice1].tipo==INCOGNITA )
                MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesU[i][indice1]), &valor, INSERT_VALUES);
            else
                ladoDireito += -valor*U[i][indice1];

            valor = -1.0;
            if( (indice2<M.Ny) && (indice2>=0) && (M.pontosU[i][indice2].tipo==SURFACE || M.pontosU[i][indice2].tipo==INCOGNITA) )
                MatSetValues(A, 1, &LinhaMatriz, 1, &(M.matIndicesU[i][indice2]), &valor, INSERT_VALUES);
            else if( (indice2==M.Ny) && (M.pontosV[i][indice2].tipoBoundary==NOSLIP) )
                ladoDireito += valor*U[i][indice2-1];
            else if( (indice2<0) && (M.pontosV[i][0].tipoBoundary==NOSLIP) )
                ladoDireito += valor*U[i][indice2+1];
            else
                ladoDireito += -valor*U[i][indice2];
        }


        VecSetValues(Vet, 1, &LinhaMatriz, &ladoDireito, INSERT_VALUES);
        (*QtdEquacoes)++;
    }
    else
        printf("\nPROBLEMA\n");


    return 1;
}

double ValorPressaoCelulaSurface(MALHA M, int i, int j, char TipoEquacao, char Direcao, double **UVelho, double **VVelho, double **UNovo, double **VNovo)
{
    static double uDirCima, uEsqCima, uDirCentro, uEsqCentro, uDirBaixo, uEsqBaixo;
    static double vDirCima, vDirBaixo, vCentroCima, vCentroBaixo, vEsqCima, vEsqBaixo;
    static double nx, ny;

    DeuProblema("ValorPressaoCelulaSurface: NAO ERA PRA ENTRAR AQUI!\n");

    if( Direcao=='c' ) {
        nx = 0.0;
        ny = 1.0;
    }
    else if( Direcao=='b' ) {
        nx = 0.0;
        ny = -1.0;
    }
    else if( Direcao=='d' ) {
        nx = 1.0;
        ny = 0.0;
    }
    else if( Direcao=='e' ) {
        nx = -1.0;
        ny = 0.0;
    }
    else if( Direcao==DIR_CIMA ) {
        nx = sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;
    }
    else if( Direcao==ESQ_CIMA ) {
        nx = -sqrt(2.0)/2.0;
        ny = sqrt(2.0)/2.0;
    }
    else if( Direcao==DIR_BAIXO ) {
        nx = sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;
    }
    else if( Direcao==ESQ_BAIXO ) {
        nx = -sqrt(2.0)/2.0;
        ny = -sqrt(2.0)/2.0;
    }
    else printf("\n\n 1: problema pressao superficie \n\n");



    ///Pontos centrais usados tanto no caso horizontal, quanto no vertical e no 45 graus
    if( M.pontosU[i+1][j].tipo==SURFACE || M.pontosU[i+1][j].tipo==FULL || M.pontosU[i+1][j].tipoBoundary==NOSLIP || M.pontosU[i+1][j].tipoBoundary==INFLOW )
        uDirCentro = UVelho[i+1][j];
    else DeuProblema("\n 2: Problema pressao superficie \n");

    if( M.pontosU[i][j].tipo==SURFACE || M.pontosU[i][j].tipo==FULL || M.pontosU[i][j].tipoBoundary==NOSLIP || M.pontosU[i][j].tipoBoundary==INFLOW )
        uEsqCentro = UVelho[i][j];
    else DeuProblema("\n 3: Problema pressao superficie \n");

    if( M.pontosV[i][j+1].tipo==SURFACE || M.pontosV[i][j+1].tipo==FULL || M.pontosV[i][j+1].tipoBoundary==NOSLIP || M.pontosV[i][j+1].tipoBoundary==INFLOW )
        vCentroCima = VVelho[i][j+1];
    else DeuProblema("\n 4: Problema pressao superficie \n");

    if( M.pontosV[i][j].tipo==SURFACE || M.pontosV[i][j].tipo==FULL || M.pontosV[i][j].tipoBoundary==NOSLIP || M.pontosV[i][j].tipoBoundary==INFLOW )
        vCentroBaixo = VVelho[i][j];
    else DeuProblema("\n 5: Problema pressao superficie %d %d\n", i, j);



    ///Demais pontos sao usados apenas no caso 45 graus
    if( TipoEquacao==1 ) {


        if( M.pontosV[i][j].tipoBoundary==NOSLIP || M.pontosV[i][j].tipoBoundary==INFLOW )
            uDirBaixo = - UVelho[i+1][j]; //Condicao de dirichlet
        else if( M.pontosU[i+1][j-1].tipo==SURFACE || M.pontosU[i+1][j-1].tipo==FULL )
            uDirBaixo = UVelho[i+1][j-1];
        else DeuProblema("\n 6: Problema pressao superficie %d %d\n", i, j);

        if( M.pontosV[i][j].tipoBoundary==NOSLIP || M.pontosV[i][j].tipoBoundary==INFLOW )
            uEsqBaixo =  - UVelho[i][j];
        else if( M.pontosU[i][j-1].tipo==SURFACE || M.pontosU[i][j-1].tipo==FULL )
            uEsqBaixo = UVelho[i][j-1];
        else DeuProblema("\n 7: Problema pressao superficie %d %d\n", i, j);

        if( M.pontosV[i][j+1].tipoBoundary==NOSLIP || M.pontosV[i][j+1].tipoBoundary==INFLOW )
            uDirCima = - UVelho[i+1][j];
        else if( M.pontosU[i+1][j+1].tipo==SURFACE || M.pontosU[i+1][j+1].tipo==FULL )
            uDirCima = UVelho[i+1][j+1];
        else DeuProblema("\n 8: Problema pressao superficie %d %d %d\n", i, j, M.pontosU[i+1][j+1].tipo);

        if( M.pontosV[i][j+1].tipoBoundary==NOSLIP || M.pontosV[i][j+1].tipoBoundary==INFLOW )
            uEsqCima = - UVelho[i][j];
        else if( M.pontosU[i][j+1].tipo==SURFACE || M.pontosU[i][j+1].tipo==FULL )
            uEsqCima = UVelho[i][j+1];
        else DeuProblema("\n 9: Problema pressao superficie %d %d\n", i, j);

        if( M.pontosU[i+1][j].tipoBoundary==NOSLIP || M.pontosU[i+1][j].tipoBoundary==INFLOW )
            vDirCima = - VVelho[i][j+1];
        else if( M.pontosV[i+1][j+1].tipo==SURFACE || M.pontosV[i+1][j+1].tipo==FULL || M.pontosV[i+1][j+1].tipoBoundary==NOSLIP || M.pontosV[i+1][j+1].tipoBoundary==INFLOW )
            vDirCima = VVelho[i+1][j+1];
        else DeuProblema("\n 10: Problema pressao superficie %d %d\n", i, j);

        if( M.pontosU[i+1][j].tipoBoundary==NOSLIP || M.pontosU[i+1][j].tipoBoundary==INFLOW )
            vDirBaixo = - VVelho[i][j];
        else if( M.pontosV[i+1][j].tipo==SURFACE || M.pontosV[i+1][j].tipo==FULL || M.pontosV[i+1][j].tipoBoundary==NOSLIP || M.pontosV[i+1][j].tipoBoundary==INFLOW )
            vDirBaixo = VVelho[i+1][j];
        else DeuProblema("\n 11: Problema pressao superficie %d %d\n", i, j);

        if( M.pontosU[i][j].tipoBoundary==NOSLIP || M.pontosU[i][j].tipoBoundary==INFLOW )
            vEsqCima = - VVelho[i][j+1];
        else if( M.pontosV[i-1][j+1].tipo==SURFACE || M.pontosV[i-1][j+1].tipo==FULL || M.pontosV[i-1][j+1].tipoBoundary==NOSLIP || M.pontosV[i-1][j+1].tipoBoundary==INFLOW )
            vEsqCima = VVelho[i-1][j+1];
        else DeuProblema("\n 12: Problema pressao superficie %d %d\n", i, j);

        if( M.pontosU[i][j].tipoBoundary==NOSLIP || M.pontosU[i][j].tipoBoundary==INFLOW )
            vEsqBaixo = - VVelho[i][j];
        else if( M.pontosV[i-1][j].tipo==SURFACE || M.pontosV[i-1][j].tipo==FULL || M.pontosV[i-1][j].tipoBoundary==NOSLIP || M.pontosV[i-1][j].tipoBoundary==INFLOW )
            vEsqBaixo = VVelho[i-1][j];
        else DeuProblema("\n 13: Problema pressao superficie %d %d %d\n\n", i, j, M.pontosV[i-1][j].tipo);
    }

    if( TipoEquacao==0 )
        return (2.0/(M.Re*(nx*nx + ny*ny)))*( ((nx*nx/M.dx[i])*(uDirCentro - uEsqCentro)) + ((ny*ny/M.dy[j])*(vCentroCima - vCentroBaixo)) );
    else
        return (2.0/(M.Re*(nx*nx + ny*ny)))*( ( (nx*nx/M.dx[i])*(uDirCentro - uEsqCentro) ) + ( (ny*ny/M.dy[j])*(vCentroCima - vCentroBaixo) )
            + ( (0.25*nx*ny/M.dy[j])*(uDirCima + uEsqCima - uDirBaixo - uEsqBaixo) )
            + ( (0.25*nx*ny/M.dx[i])*(vDirCima + vDirBaixo - vEsqCima - vEsqBaixo) ) );
}

void VetorNormalForaDaInterface(MALHA M, double PontoX, double PontoY, double *NormalX, double *NormalY)
{
    CURVA *curva;
    PONTO *ponto;
    double dist;
    double totalDist;

    DeuProblema("\n\n VetorNormalForaDaInterface: Essa funcao precisa ser refeita pra malha nao-uniforme!!!! \n\n");

    *NormalX = *NormalY = 0.0;
    totalDist = 0.0;
    for( curva=M.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;

            dist = FuncaoDistribuicao(PontoX, PontoY, ponto, M.dx[0], M.dy[0]);
            totalDist += dist;
            (*NormalX) += dist*(ponto->vetNormalX);
            (*NormalY) += dist*(ponto->vetNormalY);
        }
    }

    (*NormalX) /= totalDist;
    (*NormalY) /= totalDist;
    totalDist = sqrt( (*NormalX)*(*NormalX) + (*NormalY)*(*NormalY) );
    (*NormalX) /= totalDist;
    (*NormalY) /= totalDist;
    return;
}

void VelocidadeEmUmaParticula(MALHA M, PONTO *P, double **U, double **V, double *VelU, double *VelV)
{
    int i, j;
    int degeneradoDir, degeneradoCima, degeneradoEsq, degeneradoBaixo;
    double meioCelula;
    double velEsqBaixo, velDirBaixo, velEsqCima, velDirCima;
    double x1, y1, x2, y2, x, y;


    x = P->x;
    y = P->y;
    EncontraCelulaDaParticula(&M, x, y, &i, &j);


    velEsqBaixo = velDirBaixo = velEsqCima = velDirCima = 0.0;
    degeneradoDir = degeneradoCima = degeneradoEsq = degeneradoBaixo = 0;

    //Interpolando o valor de U
    meioCelula = M.y[j] + 0.5*M.dy[j];
    if( y>=meioCelula ) {

        if( M.pontosU[i][j].tipo==FULL || M.pontosU[i][j].tipo==SURFACE || M.pontosU[i][j].tipo==BOUNDARY )
            velEsqBaixo = U[i][j];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i+1][j]==SURFACE )
            degeneradoEsq = 1;
        else {
            DesenhaMalhaVTK(M, 0);
            DeuProblema("\n\n 1: PROBLEMA VELOCIDADE PARTICULA %d %d %lf %lf\n\n", i, j, x, y);
        }

        if( M.pontosU[i+1][j].tipo==FULL || M.pontosU[i+1][j].tipo==SURFACE || M.pontosU[i+1][j].tipo==BOUNDARY )
            velDirBaixo = U[i+1][j];
        else if( M.celulas[i][j]==EMPTY && M.pontosU[i][j].tipo!=EMPTY ) //essa celula virou empty pq tratei um caso degenerado
            degeneradoDir = 1;
        else {
//            ImprimeInterfaceVTK(M, 100000);
            DesenhaMalhaVTK(M, 0);
            DeuProblema("\n\n 2: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);
        }

        if( M.pontosV[i][j+1].tipoBoundary==NOSLIP || M.pontosV[i][j+1].tipoBoundary==INFLOW )
            velEsqCima = - U[i][j];
        else if( M.pontosV[i][j+1].tipoBoundary==SIMETRIA )
            velEsqCima = U[i][j];
        else if( M.pontosU[i][j+1].tipo==FULL || M.pontosU[i][j+1].tipo==SURFACE || M.pontosU[i][j+1].tipo==BOUNDARY )
            velEsqCima = U[i][j+1];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i+1][j]==SURFACE )
            degeneradoEsq = 1;
        else if( M.celulas[i][j]==EMPTY && M.celulas[i][j-1]==SURFACE )
            degeneradoCima = 1;
        else if( M.celulas[i][j]==SURFACE )
            degeneradoCima = 1;
        else {
//            ImprimeInterfaceVTK(M, 100000);
            DesenhaMalhaVTK(M, 0);
            DeuProblema("\n\n 3: PROBLEMA VELOCIDADE PARTICULA %d %d %lf %lf\n\n", i, j, x, y);
        }

        if( M.pontosV[i][j+1].tipoBoundary==NOSLIP || M.pontosV[i][j+1].tipoBoundary==INFLOW )
            velDirCima = - U[i+1][j];
        else if( M.pontosV[i][j+1].tipoBoundary==SIMETRIA )
            velDirCima = U[i+1][j];
        else if( M.pontosU[i+1][j+1].tipo==FULL || M.pontosU[i+1][j+1].tipo==SURFACE || M.pontosU[i+1][j+1].tipo==BOUNDARY )
            velDirCima = U[i+1][j+1];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i-1][j]==SURFACE ) //essa celula virou empty pq tratei um caso degenerado
            degeneradoDir = 1;
        else if( M.celulas[i][j]==EMPTY && M.celulas[i][j-1]==SURFACE )
            degeneradoCima = 1;
        else if( M.pontosU[i+1][j].tipo!=EMPTY )
            degeneradoCima = 1;
        else DeuProblema("\n\n 4: PROBLEMA VELOCIDADE PARTICULA %d %d\n\n", i, j);

        x1 = M.x[i];
        x2 = M.x[i+1];
        y1 = M.y[j] + 0.5*M.dy[j];
        y2 = (j==M.Ny-1) ? M.y[j+1] + 0.5*M.dy[j] : M.y[j+1] + 0.5*M.dy[j+1];
    }
    else {

        if( (j!=0) && (M.pontosU[i][j-1].tipo==FULL || M.pontosU[i][j-1].tipo==SURFACE || M.pontosU[i][j-1].tipo==BOUNDARY) )
            velEsqBaixo = U[i][j-1];
        else if( M.pontosV[i][j].tipoBoundary==NOSLIP || M.pontosV[i][j].tipoBoundary==INFLOW )
            velEsqBaixo = - U[i][j];
        else if( M.pontosV[i][j].tipoBoundary==SIMETRIA )
            velEsqBaixo = U[i][j];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i+1][j]==SURFACE )
            degeneradoEsq = 1;
        else if( M.celulas[i][j]==EMPTY && M.pontosU[i][j].tipo!=EMPTY )
            degeneradoBaixo = 1;
        else {
            DeuProblema("\n\n 5: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);
        }

        if( (j!=0) && (M.pontosU[i+1][j-1].tipo==FULL || M.pontosU[i+1][j-1].tipo==SURFACE || M.pontosU[i+1][j-1].tipo==BOUNDARY) )
            velDirBaixo = U[i+1][j-1];
        else if( M.pontosV[i][j].tipoBoundary==NOSLIP || M.pontosV[i][j].tipoBoundary==INFLOW )
            velDirBaixo = - U[i+1][j];
        else if( M.pontosV[i][j].tipoBoundary==SIMETRIA )
            velDirBaixo = U[i+1][j];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i-1][j]==SURFACE ) //essa celula virou empty pq tratei um caso degenerado
            degeneradoDir = 1;
        else if( M.celulas[i][j]==EMPTY && M.pontosU[i+1][j].tipo!=EMPTY )
            degeneradoBaixo = 1;
        else DeuProblema("\n\n 6: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);

        if( M.pontosU[i][j].tipo==FULL || M.pontosU[i][j].tipo==SURFACE || M.pontosU[i][j].tipo==BOUNDARY )
            velEsqCima = U[i][j];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i+1][j]==SURFACE )
            degeneradoEsq = 1;
        else if( M.celulas[i][j]==EMPTY && M.pontosU[i][j-1].tipo==SURFACE )
            degeneradoCima = 1;
        else if( M.celulas[i][j]==EMPTY && M.pontosU[i][j-1].tipo==BOUNDARY )
            degeneradoCima = 1;
        else {
            DeuProblema("\n\n 7: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);
        }

        if( M.pontosU[i+1][j].tipo==FULL || M.pontosU[i+1][j].tipo==SURFACE || M.pontosU[i+1][j].tipo==BOUNDARY )
            velDirCima = U[i+1][j];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i-1][j]==SURFACE ) //essa celula virou empty pq tratei um caso degenerado
            degeneradoDir = 1;
        else if( M.celulas[i][j]==EMPTY && M.pontosU[i+1][j-1].tipo==SURFACE )
            degeneradoCima = 1;
        else if( M.pontosU[i+1][j-1].tipo!=EMPTY )
            degeneradoCima = 1;
        else {
            DeuProblema("\n\n 8: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);
        }

        x1 = M.x[i];
        x2 = M.x[i+1];
        y1 = (j==0) ? M.y[j] - 0.5*M.dy[j] : M.y[j-1] + 0.5*M.dy[j-1];
        y2 = M.y[j] + 0.5*M.dy[j];
    }

    //Corrigindo os degenerados
    if( degeneradoDir ) {
        velDirBaixo = velEsqBaixo;
        velDirCima = velEsqCima;
    }
    else if( degeneradoEsq ) {
        velEsqBaixo = velDirBaixo;
        velEsqCima = velDirCima;
    }
    else if( degeneradoCima ) {
        velEsqCima = velEsqBaixo;
        velDirCima = velDirBaixo;
    }
    else if( degeneradoBaixo ) {
        velDirBaixo = velDirCima;
        velEsqBaixo = velEsqCima;
    }

    //Formula da extrapolacao bilinear
    *VelU = (1.0/((x2-x1)*(y2-y1)))*( ((x2-x)*((y2-y)*velEsqBaixo + (y-y1)*velEsqCima)) + ((x-x1)*((y2-y)*velDirBaixo + (y-y1)*velDirCima)) );



    degeneradoDir = degeneradoCima = degeneradoEsq = degeneradoBaixo = 0;
    //Interpolando o valor de V
    meioCelula = M.x[i] + 0.5*M.dx[i];
    if( x>=meioCelula ) {

        if( M.pontosV[i][j].tipo==FULL || M.pontosV[i][j].tipo==SURFACE || M.pontosV[i][j].tipo==BOUNDARY )
            velEsqBaixo = V[i][j];
        else if( M.celulas[i][j]==EMPTY && M.pontosV[i][j+1].tipo!=EMPTY )
            degeneradoBaixo = 1;
//        else
//            velEsqBaixo = 0.0; //Nao sei se poderia fazer isso...
        else {
//            ImprimeInterfaceVTK(M, 100000);
            DesenhaMalhaVTK(M, 0);
            DeuProblema("\n\n 9: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);
        }

        if( M.pontosU[i+1][j].tipoBoundary==NOSLIP || M.pontosU[i+1][j].tipoBoundary==INFLOW )
            velDirBaixo = - V[i][j];
        else if( M.pontosU[i+1][j].tipoBoundary==NEUMANN || M.pontosU[i+1][j].tipoBoundary==SIMETRIA )
            velDirBaixo = V[i][j];
        else if( M.pontosU[i+1][j].tipoBoundary==PERIODICIDADE )
            velDirBaixo = V[0][j];
        else if( M.pontosV[i+1][j].tipo==FULL || M.pontosV[i+1][j].tipo==SURFACE || M.pontosV[i+1][j].tipo==BOUNDARY )
            velDirBaixo = V[i+1][j];
        else if( M.celulas[i][j]==EMPTY && M.pontosV[i+1][j+1].tipo!=EMPTY )
            degeneradoBaixo = 1;
        else DeuProblema("\n\n 10: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);

        if( M.pontosV[i][j+1].tipo==FULL || M.pontosV[i][j+1].tipo==SURFACE || M.pontosV[i][j+1].tipo==BOUNDARY )
            velEsqCima = V[i][j+1];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i][j-1]==SURFACE )
            degeneradoCima = 1;
//        else
//            velEsqCima = 0.0; //Nao sei nao...
        else DeuProblema("\n\n 11: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);

        if( M.pontosU[i+1][j].tipoBoundary==NOSLIP || M.pontosU[i+1][j].tipoBoundary==INFLOW )
            velDirCima = - V[i][j+1];
        else if( M.pontosU[i+1][j].tipoBoundary==NEUMANN || M.pontosU[i+1][j].tipoBoundary==SIMETRIA )
            velDirCima = V[i][j+1];
        else if( M.pontosU[i+1][j].tipoBoundary==PERIODICIDADE )
            velDirCima = V[0][j+1];
        else if( M.pontosV[i+1][j+1].tipo==FULL || M.pontosV[i+1][j+1].tipo==SURFACE || M.pontosV[i+1][j+1].tipo==BOUNDARY )
            velDirCima = V[i+1][j+1];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i][j-1]==SURFACE )
            degeneradoCima = 1;
        else {
            DeuProblema("\n\n 12: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);
        }

        x1 = M.x[i] + 0.5*M.dx[i];
        x2 = (i==M.Nx-1) ? M.x[i+1] + 0.5*M.dx[i] : M.x[i+1] + 0.5*M.dx[i+1];
        y1 = M.y[j];
        y2 = M.y[j+1];
    }
    else {

        if( M.pontosU[i][j].tipoBoundary==NOSLIP || M.pontosU[i][j].tipoBoundary==INFLOW )
            velEsqBaixo = - V[i][j];
        else if( M.pontosU[i][j].tipoBoundary==NEUMANN || M.pontosU[i][j].tipoBoundary==SIMETRIA || M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA )
            velEsqBaixo = V[i][j];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
            velEsqBaixo = V[M.Nx-1][j];
        else if( M.pontosV[i-1][j].tipo==FULL || M.pontosV[i-1][j].tipo==SURFACE || M.pontosV[i-1][j].tipo==BOUNDARY )
            velEsqBaixo = V[i-1][j];
        else if( M.celulas[i][j]==EMPTY && M.pontosV[i-1][j+1].tipo!=EMPTY )
            degeneradoBaixo = 1;
//        else
//            velEsqBaixo = 0.0; //Nao eh bom
        else DeuProblema("\n\n 13: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);

        if( M.pontosV[i][j].tipo==FULL || M.pontosV[i][j].tipo==SURFACE || M.pontosV[i][j].tipo==BOUNDARY )
            velDirBaixo = V[i][j];
        else if( M.celulas[i][j]==EMPTY && M.pontosV[i][j+1].tipo!=EMPTY )
            degeneradoBaixo = 1;
//        else
//            velDirBaixo = 0.0;
        else DeuProblema("\n\n 14: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);

        if( M.pontosU[i][j].tipoBoundary==NOSLIP || M.pontosU[i][j].tipoBoundary==INFLOW )
            velEsqCima = - V[i][j+1];
        else if( M.pontosU[i][j].tipoBoundary==NEUMANN || M.pontosU[i][j].tipoBoundary==SIMETRIA || M.pontosU[i][j].tipoBoundary==EIXO_AXISSIMETRIA )
            velEsqCima = V[i][j+1];
        else if( M.pontosU[i][j].tipoBoundary==PERIODICIDADE )
            velEsqCima = V[M.Nx-1][j+1];
        else if( M.pontosV[i-1][j+1].tipo==FULL || M.pontosV[i-1][j+1].tipo==SURFACE || M.pontosV[i-1][j+1].tipo==BOUNDARY )
            velEsqCima = V[i-1][j+1];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i][j-1]==SURFACE )
            degeneradoCima = 1;
//        else
//            velEsqCima = 0.0;
        else DeuProblema("\n\n 15: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);

        if( M.pontosV[i][j+1].tipo==FULL || M.pontosV[i][j+1].tipo==SURFACE || M.pontosV[i][j+1].tipo==BOUNDARY )
            velDirCima = V[i][j+1];
        else if( M.celulas[i][j]==EMPTY && M.celulas[i][j-1]==SURFACE )
            degeneradoCima = 1;
//        else
//            velDirCima = 0.0;
        else DeuProblema("\n\n 16: PROBLEMA VELOCIDADE PARTICULA %d %d \n\n", i, j);

        x1 = (i==0) ? M.x[i] - 0.5*M.dx[i] : M.x[i-1] + 0.5*M.dx[i-1];
        x2 = M.x[i] + 0.5*M.dx[i];
        y1 = M.y[j];
        y2 = M.y[j+1];
    }

    //Casos degenerados
    if( degeneradoCima ) {
        velEsqCima = velEsqBaixo;
        velDirCima = velDirBaixo;
    }
    else if( degeneradoBaixo ) {
        velEsqBaixo = velEsqCima;
        velDirBaixo = velDirCima;
    }

    *VelV = (1.0/((x2-x1)*(y2-y1)))*( ((x2-x)*((y2-y)*velEsqBaixo + (y-y1)*velEsqCima)) + ((x-x1)*((y2-y)*velDirBaixo + (y-y1)*velDirCima)) );


    return;
}

double FuncaoDistribuicao(double X, double Y, PONTO *P, double dx, double dy)
{
    if( (fabs(X-P->x)<2*dx) && (fabs(Y-P->y)<2*dy) )
        return (( 1.0/(4.0*dx) )*( 1.0/(4.0*dy) )*( 1+cos((0.5*M_PI/dx)*(X-P->x)) )*( 1+cos((0.5*M_PI/dy)*(Y-P->y)) ));
    return 0.0;
}

void EncontraCelulaSeedFloodFill(MALHA M, LISTA_CELULAS_SURFACE L, int *iSeed, int *jSeed)
{
    CELULA_SURFACE *celula1, *celula2;
    int achouSurface;
    int i, j;

    for( celula1=L.prim; celula1!=NULL; celula1=celula1->prox ) {
        j = celula1->j;

        //Vai indo pra direita ate achar uma empty
        for( i=celula1->i+1; i<M.Nx; i++ ) {
            if( M.celulas[i][j]==EMPTY )
                break;
        }

        *iSeed = i;
        *jSeed = j;

        //Vai continuando pra direita ate achar uma surface de novo
        for( ; i<M.Nx; i++ ) {
            if( M.celulas[i][j]==SURFACE )
                break;
        }


        //Nao encontrou surface
        if( i==M.Nx )
            continue;


        //Verifica se essa surface encontrada REALMENTE pertence a essa lista
        for( celula2=L.prim; celula2!=NULL; celula2=celula2->prox ) {
            if( i==celula2->i && j==celula2->j )
                break;
        }


        //Se pertencer a lista, PROVAVALMENTE encontramos o seed.
        //Pra evitar bugs q vem com caso degenerado. vou fazer uma passagem na direcao y pra garantir mesmo q essa celula ta dentro da regiao
        if( celula2!=NULL ) {
            i = *iSeed;
            achouSurface = 0;

            j = 0;
            while( j<M.Ny ) {

                if( achouSurface==0 ) {
                    if( M.celulas[i][j]==SURFACE ) {
                        //Verifica se esta na lista
                        for( celula2=L.prim; celula2!=NULL; celula2=celula2->prox ) {
                            if( i==celula2->i && j==celula2->j ) {
                                achouSurface = 1;

                                //Continua pra cima ate chegar numa empty
                                for( j++; j<M.Ny; j++ ) {
                                    //if( M.celulas[i][j]==SURFACE ) //Vou desconsiderar esse aso
                                      //  achouSurface = 0;

                                    if( M.celulas[i][j]==EMPTY )
                                        break;
                                }
                                break;
                            }
                        }
                    }
                    if( achouSurface )
                        continue;
                }
                else if( achouSurface==1 ) {

                    if( M.celulas[i][j]==SURFACE ) {
                        achouSurface = 0;
                        for( j++; j<M.Ny; j++ ) {
                            if( M.celulas[i][j]==EMPTY )
                                break;
                        }
                        continue;
                    }

                    //ENCONTROU O SEED! 
                    if( *iSeed==i && *jSeed==j ) {
                        return;
                    }
                }

                j++;
            }
        }
    }

    *iSeed = -1;
    *jSeed = -1;

    return;
}

void EncontraVelPressaoSeedFloodFill(MALHA M, int *iSeedU, int *jSeedU, int *iSeedV, int *jSeedV, int *iSeedP, int *jSeedP)
{
    int i, j, encontrouU, encontrouV, encontrouP;

    encontrouU = encontrouV = encontrouP = 0;
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            if( M.celulas[i][j]!=FULL )
                continue;

            if( !encontrouP && M.pontosP[i][j].tipo==EMPTY ) {
                *iSeedP = i;
                *jSeedP = j;
                encontrouP = 1;
            }

            if( !encontrouU ) {
                if( M.pontosU[i][j].tipo==EMPTY ) {
                    *iSeedU = i;
                    *jSeedU = j;
                    encontrouU = 1;
                }
                else if( M.pontosU[i+1][j].tipo==EMPTY ) {
                    *iSeedU = i+1;
                    *jSeedU = j;
                    encontrouU = 1;
                }
            }

            if( !encontrouV ) {
                if( M.pontosV[i][j].tipo==EMPTY ) {
                    *iSeedV = i;
                    *jSeedV = j;
                    encontrouV = 1;
                }
                else if( M.pontosV[i][j+1].tipo==EMPTY ) {
                    *iSeedV = i;
                    *jSeedV = j+1;
                    encontrouV = 1;
                }
            }

            if( encontrouU && encontrouV && encontrouP )
                return;
        }
    }

    if( !encontrouU ) {
        *iSeedU = -1;
        *jSeedU = -1;
    }
    if( !encontrouV ) {
        *iSeedV = -1;
        *jSeedV = -1;
    }
    if( !encontrouP ) {
        *iSeedP = -1;
        *jSeedP = -1;
    }

    return;
}

int AdicionaCelulaSurface(MALHA *M, LISTA_CELULAS_SURFACE *L, int i, int j, int Passo)
{
    CELULA_SURFACE *celula;

    //Verifica se esse (i, j) ja esta na lista
    for( celula=L->ult; celula!=NULL; celula=celula->ant ) {
        if( (celula->i==i) && (celula->j==j) ) {
            if( celula->passoCelula!=Passo )
                celula->novaCelulaSurface = 0;

            celula->passoCelula = Passo;
            return 0;
        }
    }

    celula = (CELULA_SURFACE *)malloc( sizeof(CELULA_SURFACE) );
    celula->i = i;
    celula->j = j;
    celula->passoCelula = Passo;
    celula->prox = NULL;
    celula->ant = NULL;

    // This cell is turning from EMPTY to SURFACE
    // I will set a flag indicating this which will be used by other functions sometimes
    celula->novaCelulaSurface = (M->celulas[i][j]==EMPTY) ? 1 : 0;

    

    if( L->prim==NULL ) {
        L->prim = celula;
        L->ult = celula;
    }
    else {
        L->ult->prox = celula;
        celula->ant = L->ult;
        L->ult = celula;
    }

    return 1;
}

void RemoveCelulaSurface(LISTA_CELULAS_SURFACE *L, CELULA_SURFACE *Celula, MALHA *M, TIPO_CELULA Tipo)
{
    int i, j;

    if( Celula==NULL )
        return;

    if( L->prim==Celula )
        L->prim = Celula->prox;
    if( L->ult==Celula )
        L->ult = Celula->ant;

    if( Celula->ant!=NULL )
        Celula->ant->prox = Celula->prox;
    if( Celula->prox!=NULL )
        Celula->prox->ant = Celula->ant;


    if( Tipo!=-1 ) {
        i = Celula->i;
        j = Celula->j;

        M->pontosP[i][j].tipo = Tipo;
        if( M->pontosU[i][j].tipo!=BOUNDARY )
            M->pontosU[i][j].tipo = Tipo;
        if( M->pontosU[i+1][j].tipo!=BOUNDARY )
            M->pontosU[i+1][j].tipo = Tipo;
        if( M->pontosV[i][j].tipo!=BOUNDARY )
            M->pontosV[i][j].tipo = Tipo;
        if( M->pontosV[i][j+1].tipo!=BOUNDARY )
            M->pontosV[i][j+1].tipo = Tipo;
    }

    free(Celula);
    return;
}

//void ModificaCelulaSurface(LISTA_CELULAS_SURFACE *L, CELULA_SURFACE *Celula, int i, int j)
//{
//    CELULA_SURFACE *celulaIgual;
//
//    celulaIgual = EncontraCelulaSurface(L, i, j);
//
//    if( celulaIgual!=NULL ) {
//        RemoveCelulaSurface(L, Celula);
//        return;
//    }
//
//    Celula->i = i;
//    Celula->j = j;
//    return;
//}

CELULA_SURFACE *EncontraCelulaSurface(LISTA_CELULAS_SURFACE *L, int i, int j)
{
    CELULA_SURFACE *Celula;

    for( Celula=L->prim; Celula!=NULL; Celula=Celula->prox ) {
        if( Celula->i==i && Celula->j==j )
            return Celula;
    }

    return NULL;
}

LISTA_CELULAS_SURFACE UniaoListasCelulaSurface(MALHA *M, LISTA_CELULAS_SURFACE *L, int QtdListas)
{
    int i = -1, j = -1, iMedio = -1, jMedio = -1, iPressao = -1, jPressao = -1;
    CELULA_SURFACE *celulaSurf, *proxCelula, *celula2;
    LISTA_CELULAS_SURFACE listaCelulas;

    //Removendo possiveis casos que fica uma celula nao-surface na lista
    for( i=0; i<QtdListas; i++ ) {

        celulaSurf = L[i].prim;
        while( celulaSurf!=NULL ) {
            proxCelula = celulaSurf->prox;

            if( M->celulas[celulaSurf->i][celulaSurf->j]!=SURFACE ) {
                RemoveCelulaSurface(&L[i], celulaSurf, M, -1);
            }

            celulaSurf = proxCelula;
        }

    }

    //Removendo repeticoes nas listas
    for( i=1; i<QtdListas; i++ ) {
        for( j=0; j<QtdListas; j++ ) {
            if( i==j )
                continue;

            celulaSurf = L[i].prim;
            while( celulaSurf!=NULL ) {
                proxCelula = celulaSurf->prox;

                for( celula2=L[j].prim; celula2!=NULL; celula2=celula2->prox ) {
                    if( celula2->i==celulaSurf->i && celula2->j==celulaSurf->j )
                        break;
                }

                if( celula2!=NULL ) {
                    RemoveCelulaSurface(&L[i], celulaSurf, M, -1);
                }

                celulaSurf = proxCelula;
            }
        }
    }

    //Realiza a uniao
    for( i=0; i<QtdListas-1; i++ ) {
        L[i].ult->prox = L[i+1].prim;
        L[i+1].prim->ant = L[i].ult;
    }
    L[0].ult = L[QtdListas-1].ult;
    listaCelulas = L[0];



    //Transformando em FULL as celulas que nao tem nenhuma empty em volta
    celulaSurf = listaCelulas.prim;
    while( celulaSurf!=NULL ) {
        i = celulaSurf->i;
        j = celulaSurf->j;
        proxCelula = celulaSurf->prox;

        if( M->celulas[i][j]!=SURFACE ) {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 2: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! %d %d\n", i, j);
        }

        if( (i==M->Nx-1 || M->celulas[i+1][j]!=EMPTY) && (i==0 || M->celulas[i-1][j]!=EMPTY) && (j==M->Ny-1 || M->celulas[i][j+1]!=EMPTY) && (j==0 || M->celulas[i][j-1]!=EMPTY) ) {
            M->celulas[i][j] = FULL;
            RemoveCelulaSurface(&listaCelulas, celulaSurf, M, FULL);

           // Testando tensao superficial na celula full
            M->celulasParede[i][j] = SURFACE;
        }

        celulaSurf = proxCelula;
    }



    //Preenchendo todo o interior com pontos U, V, P FULL
    //Primeiro fazendo a borda
    for( celulaSurf=listaCelulas.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox ) {
        i = celulaSurf->i;
        j = celulaSurf->j;

        if( M->celulas[i][j]!=SURFACE ) {
            //DesenhaMalhaVTK(*M, 0);
            DeuProblema("\n 3: A CELULA NAO EH SURFACE!!!!! PROBLEMA!!!! TUDO BUGADO!!!! \n");
        }

        M->pontosP[i][j].tipo = SURFACE;

        //Superficie vertical: empty na direita
        if( LADO_VERTICAL_DIREITA(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        //Superficie vertical: empty na esquerda
        else if( LADO_VERTICAL_ESQUERDA(i, j) ) {
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        //Superficie horizontal: empty em cima
        else if( LADO_HORIZONTAL_CIMA(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
        }
        //Superficie horizontal: empty embaixo
        else if( LADO_HORIZONTAL_BAIXO(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        //Superficie quina: empty em cima e na direita
        else if( QUINA_CIMA_DIREITA(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
        }
        //Superficie quina: empty em cima e na esquerda
        else if( QUINA_CIMA_ESQUERDA(i, j) ) {
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j].tipo!=BOUNDARY )
                M->pontosV[i][j].tipo = FULL;
        }
        //Superficie quina: empty em embaixo e na direita
        else if( QUINA_BAIXO_DIREITA(i, j) ) {
            if( M->pontosU[i][j].tipo!=BOUNDARY )
                M->pontosU[i][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        //Superficie quina: empty em embaixo e na esquerda
        else if( QUINA_BAIXO_ESQUERDA(i, j) ) {
            if( M->pontosU[i+1][j].tipo!=BOUNDARY )
                M->pontosU[i+1][j].tipo = FULL;
            if( M->pontosV[i][j+1].tipo!=BOUNDARY )
                M->pontosV[i][j+1].tipo = FULL;
        }
        else {
            DesenhaMalhaVTK(*M, 0);
            //ImprimeInterfaceVTK(*M, 1000000);
            DeuProblema("\n\n BBB SUP LIVRE: CASO INESPERADO %d %d \n\n", i, j);
        }
    }



    do {
        EncontraVelPressaoSeedFloodFill(*M, &iMedio, &jMedio, &i, &j, &iPressao, &jPressao);
        AlgoritmoFloodFillUIterativo(*M, iMedio, jMedio, FULL);
        AlgoritmoFloodFillVIterativo(*M, i, j, FULL);
        AlgoritmoFloodFillPIterativo(*M, iPressao, jPressao, FULL);
    } while( iMedio!=-1 || jMedio!=-1 || i!=-1 || j!=-1 || iPressao!=-1 || jPressao!=-1 );


    return listaCelulas;
}

void DestroiListaCelulasSurface(LISTA_CELULAS_SURFACE *L)
{
    CELULA_SURFACE *celula, *prox;

    celula = L->prim;
    while( celula!=NULL ) {
        prox = celula->prox;
        free(celula);
        celula = prox;
    }

    L->prim = NULL;
    L->ult = NULL;
    return;
}

int CasoDegeneradoCima(MALHA M, int i, int j)
{
    if( (i==M.Nx-1 || M.celulas[i+1][j]==EMPTY) && (i==0 || M.celulas[i-1][j]==EMPTY) && (j==M.Ny-1 || M.celulas[i][j+1]==EMPTY) )
        return 1;
    if( (i==M.Nx-1 || M.celulas[i+1][j]==FULL) && (i==0 || M.celulas[i-1][j]==FULL) && (j==M.Ny-1 || M.celulas[i][j+1]==FULL) )
        return 1;
    return 0;
}

void DeletaCurvasFULL(MALHA *M)
{
    int i, j;
    CURVA *c, *proxCurva;
    PONTO *p;
    int encontrouNaoFULL;

    c = M->interface.curvaInicial;
    while( c ) {
        proxCurva = c->proxCurva;

        encontrouNaoFULL = 0;
        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);
            if( M->celulas[i][j]!=FULL ) {
                encontrouNaoFULL = 1;
                break;
            }
        }

        // Remove essa curva
        if( !encontrouNaoFULL ) {
            if( c->curvaAnt )
                c->curvaAnt->proxCurva = proxCurva;
            if( proxCurva )
                proxCurva->curvaAnt = c->curvaAnt;
            if( M->interface.curvaInicial==c )
                M->interface.curvaInicial = proxCurva;
            LiberaMemoriaCurva(c);
        }

        c = proxCurva;
    }

    return;
}

void ArrumaMudancaTopologicaAxissimetrica(MALHA *M)
{
    PONTO *p1 = NULL, *p2 = NULL, *p;
    int i, j;
    CURVA *c;
//    int apenasFULL;

    for( c=M->interface.curvaInicial; c; c=c->proxCurva ) {

        LOOP_CURVA(c) {
            p = loop_stuff.ponto;

            // Encontrando o primeiro ponto colado no eixo de simetria
            if( (p1==NULL) && (p->x<1e-7) && (p->prox->x>=1e-7) ) {
                EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);
                int tipo1 = M->celulas[i][j];
                EncontraCelulaDaParticula(M, p->prox->x, p->prox->y, &i, &j);
                int tipo2 = M->celulas[i][j];

                if( tipo1==FULL && tipo2==FULL ) {
                    p1 = p;
                    continue;
                }
            }

            // Encontrando o segundo ponto
            if( p1 && !p2 ) {
                EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);
                int tipo = M->celulas[i][j];
                if( tipo!=FULL ) {
                    p1 = NULL;
                    continue;
                }

                if( p->x<1e-7 )
                    p2 = p;

            }

            // Eliminando os intermediarios entre p1 e p2
            if( p1 && p2 ) {
                // DeuProblema("P1 P2: %lf %lf %lf %lf\n", p1->x, p1->y, p2->x, p2->y);
                PONTO *proxPonto = NULL;
                p = p1->prox;
                while( p!=p2 ) {
                    proxPonto = p->prox;
                    free(p);
                    p = proxPonto;
                }

                p1->prox = p2;
                p2->ant = p1;

//                DeuProblema("%e %e %p ---- %e %e %p\n", p1->x, p1->y, p1, p2->x, p2->y, p2);
                break;
            }


        }
    }
}

int DeletaCurvasPequenasCompletamenteSURFACE(MALHA *M)
{
    int i, j, deletou = 0;
    CURVA *curva, *proxCurva;

    for( i=0; i<M->Nx; i++ )
        for( j=0; j<M->Ny; j++ )
            M->malhaRecursoes[i][j] = 0;

    // Verifica se cada uma das curvas precisa ser removida por esse criterio ou nao
    curva = M->interface.curvaInicial;
    while( curva ) {
        proxCurva = curva->proxCurva;

        PONTO *p = curva->pontoInicial;
        EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);
        int deleta = RecursaoDeletaCurvasSURFACE(M, curva, i, j);

//        PrintDebug("PONTO %lf %lf\n", p->x, p->y);
//        PrintDebug("DELETA: %d\n", deleta);
//        DeuProblema("PAROU!!\n\n");
        // Remove essa curva
        if( deleta ) {
            if( curva->curvaAnt )
                curva->curvaAnt->proxCurva = proxCurva;
            if( proxCurva )
                proxCurva->curvaAnt = curva->curvaAnt;
            if( M->interface.curvaInicial==curva )
                M->interface.curvaInicial = proxCurva;
            LiberaMemoriaCurva(curva);

            PrintDebug("DELETOU CURVA SURFACE\n");

            deletou = 1;
        }

        curva = proxCurva;
    }


    return deletou;
}

/// === Esta funcao retorna TRUE se soh existir celulas SURFACE dentro dessa curva
/// === Retorna FALSE se existir alguma celula FULL dentro
int RecursaoDeletaCurvasSURFACE(MALHA *M, CURVA *curva, int i, int j)
{
    int valorCentro, valorEsq, valorDir, valorBaixo, valorCima;

    if( (i<0) || (j<0) || (i>=M->Nx) || (j>=M->Ny) )
        return 1;

    // Se ja tiver visitado essa celula, retorna TRUE
    if( M->malhaRecursoes[i][j] )
        return 1;

    // Fora do dominio
    if( (i<0) && (i>=M->Nx) && (j<0) && (j>=M->Ny) )
        return 1;

    // Marca esta celula como visitada
    M->malhaRecursoes[i][j] = 1;

    // Centro dessa celula
    double xCentro = 0.5*( M->x[i] + M->x[i+1] );
    double yCentro = 0.5*( M->y[j] + M->y[j+1] );

    // Verifica se esta celula possui alguma particula dentro dela
    PONTO *p;
    int dentro = 0;
    LOOP_CURVA(curva) {
        p = loop_stuff.ponto;

        int iP, jP;

        EncontraCelulaDaParticula(M, p->x, p->y, &iP, &jP);
        if( (iP==i) && (jP==j) ) {
            dentro = 1;
            break;
        }
    }

    // Se nao tiver encontrado, verifica se o centro dessa celula esta dentro da curva
    if( (!dentro) && (is_point_in_path(xCentro, yCentro, curva)) )
        dentro = 1; // Considera essa celula

    // Ou seja, esta celula esta fora da interface. Vamos ignora-la
    if( !dentro )
        return 1;

    // Vamos usar essa celula
    valorCentro = (M->celulas[i][j]==SURFACE) ? 1 : 0;
    valorCima = RecursaoDeletaCurvasSURFACE(M, curva, i, j+1);
    valorBaixo = RecursaoDeletaCurvasSURFACE(M, curva, i, j-1);
    valorDir = RecursaoDeletaCurvasSURFACE(M, curva, i+1, j);
    valorEsq = RecursaoDeletaCurvasSURFACE(M, curva, i-1, j);

//    PrintDebug("CCCC VALOR CENTRO %d %d %d\n", i, j, valorCentro);

    return valorCentro && valorCima && valorBaixo && valorDir && valorEsq;
}

int is_point_in_path(double x, double y, CURVA *curva)
{
//    # Determine if the point is in the polygon.
//    #
//    # Args:
//    #   x -- The x coordinates of point.
//    #   y -- The y coordinates of point.
//    #   poly -- a list of tuples [(x, y), (x, y), ...]
//    #
//    # Returns:
//    #   True if the point is in the path or is a corner or on the boundary

    int c = 0; // False

    PONTO *pontoI;
    LOOP_CURVA(curva) {
        pontoI = loop_stuff.ponto;

        PONTO *pontoJ = pontoI->ant;
        if( (x==pontoI->x) && (y==pontoI->y) )
            return 1; //True. Point is exactly a corner point

        if( (pontoI->y > y) != (pontoJ->y > y) ) {
            double slope = (x - pontoI->x)*(pontoJ->y - pontoI->y)-(pontoJ->x-pontoI->x)*(y-pontoI->y);
            if( slope==0 )
                return 1; //True. Point is exactly on a boundary
            if( (slope < 0) != (pontoJ->y < pontoI->y) )
                c = !c;
        }
    }

    return c;
}

void SalvarEstado(MALHA M, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, LISTA_CELULAS_SURFACE ListaSurf, int n)
{
    FILE *arq;
    char nomeArq[500];
    int i, j, qtdCurvas, qtdPontos;
    PONTO *ponto;
    CURVA *curva;
    CELULA_SURFACE *celulaSurf;

    sprintf(nomeArq, "ArquivosPlot/VTK/%s/EstadoPasso.txt", M.nomeArquivo);
    arq = fopen(nomeArq, "wt");

    fprintf(arq, "Passo %d\n", n);
    fprintf(arq, "veloc U\n");
    for( i=0; i<=M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fprintf(arq, "%.60e ", U[i][j]);
        }
    }
    fprintf(arq, "\nveloc V\n");
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<=M.Ny; j++ ) {
            fprintf(arq, "%.60e ", V[i][j]);
        }
    }
    fprintf(arq, "Pressao\n");
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fprintf(arq, "%.60e ", P[i][j]);
        }
    }
    fprintf(arq, "Txx\n");
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fprintf(arq, "%.60e ", Txx[i][j]);
        }
    }
    fprintf(arq, "Txy\n");
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fprintf(arq, "%.60e ", Txy[i][j]);
        }
    }
    fprintf(arq, "Tyy\n");
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fprintf(arq, "%.60e ", Tyy[i][j]);
        }
    }



    //Contando qtas curvas tem
    qtdCurvas = 0;
    for( curva=M.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva )
        qtdCurvas++;

    fprintf(arq, "\nInterface qtdCurvas=%d\n", qtdCurvas);


    //Vou percorrer ate a ultima curva pra imprimir ao contrario
    for( curva=M.interface.curvaInicial; curva->proxCurva!=NULL; curva=curva->proxCurva );

    for( ; curva!=NULL; curva=curva->curvaAnt ) {
        //Busca a curva com este indice na lista (eu quero imprimir de tras pra frente)
        //curva = M.interface.curvaInicial;
        //for( j=0; j<i; j++ )
         //   curva = curva->proxCurva;

        //Conta quantos pontos tem nessa curva
        qtdPontos = 0;
        LOOP_CURVA(curva)
            qtdPontos++;

        fprintf(arq, "CURVA qtdPontos=%d regiaoEsq=%d regiaoDir=%d localEsq=%d localDir=%d\n", qtdPontos, curva->regiaoEsq->numero, curva->regiaoDir->numero, curva->regiaoEsq->localRegiao, curva->regiaoDir->localRegiao);
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;
            fprintf(arq, "[%.15e %.15e %d] ", ponto->x, ponto->y, ponto->fixo);
        }

        fprintf(arq, "\n");
    }

    fprintf(arq, "Celulas\n");
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fprintf(arq, "%u ", M.celulas[i][j]);
        }
    }

    fprintf(arq, "\nTipo U\n");
    for( i=0; i<=M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fprintf(arq, "%u ", M.pontosU[i][j].tipo);
        }
    }
    fprintf(arq, "\nTipo V\n");
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<=M.Ny; j++ ) {
            fprintf(arq, "%u ", M.pontosV[i][j].tipo);
        }
    }
    fprintf(arq, "\nTipo P\n");
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fprintf(arq, "%u ", M.pontosP[i][j].tipo);
        }
    }


    qtdPontos = 0;
    for( celulaSurf=ListaSurf.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox )
        qtdPontos++;
    fprintf(arq, "\nListaSurf %d\n", qtdPontos);
    for( celulaSurf=ListaSurf.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox )
        fprintf(arq, "[%d %d %d] ", celulaSurf->i, celulaSurf->j, celulaSurf->passoCelula);


    fclose(arq);
}

extern double *hVelho;

void SalvarEstadoBinario(MALHA M, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, double **Lambda, double **Mu, double **Nu, LISTA_CELULAS_SURFACE ListaSurf, int n)
{
    FILE *arq;
    char nomeArq[500];
    int i, j, qtdCurvas, qtdPontos;
    PONTO *ponto;
    CURVA *curva;
    CELULA_SURFACE *celulaSurf;

    sprintf(nomeArq, "ArquivosPlot/VTK/%s/EstadoPasso.txt", M.nomeArquivo);
    arq = fopen(nomeArq, "w");

    // Passo
    fwrite(&n, sizeof(int), 1, arq);
    // Veloc U
    for( i=0; i<=M.Nx; i++ )
        fwrite(U[i], sizeof(double), M.Ny, arq);
    // Veloc V
    for( i=0; i<M.Nx; i++ )
        fwrite(V[i], sizeof(double), M.Ny+1, arq);
    // Pressao
    for( i=0; i<M.Nx; i++ )
        fwrite(P[i], sizeof(double), M.Ny, arq);
    // Txx
    for( i=0; i<M.Nx; i++ )
        fwrite(Txx[i], sizeof(double), M.Ny, arq);
    // Txy
    for( i=0; i<M.Nx; i++ )
        fwrite(Txy[i], sizeof(double), M.Ny, arq);
    // Tyy
    for( i=0; i<M.Nx; i++ )
        fwrite(Tyy[i], sizeof(double), M.Ny, arq);
    // Lambda
    for( i=0; i<M.Nx; i++ )
        fwrite(Lambda[i], sizeof(double), M.Ny, arq);
    // Mu
    for( i=0; i<M.Nx; i++ )
        fwrite(Mu[i], sizeof(double), M.Ny, arq);
    // Nu
    for( i=0; i<M.Nx; i++ )
        fwrite(Nu[i], sizeof(double), M.Ny, arq);


    if( M.interface.curvaInicial!=NULL ) {
        //Contando qtas curvas tem
        qtdCurvas = 0;
        for( curva=M.interface.curvaInicial; curva!=NULL; curva=curva->proxCurva )
            qtdCurvas++;

        //qtdCurvas
        fwrite(&qtdCurvas, sizeof(int), 1, arq);


        //Vou percorrer ate a ultima curva pra imprimir ao contrario
        for( curva=M.interface.curvaInicial; curva->proxCurva!=NULL; curva=curva->proxCurva );

        for( ; curva!=NULL; curva=curva->curvaAnt ) {
            //Busca a curva com este indice na lista (eu quero imprimir de tras pra frente)
            //curva = M.interface.curvaInicial;
            //for( j=0; j<i; j++ )
             //   curva = curva->proxCurva;

            //Conta quantos pontos tem nessa curva
            qtdPontos = 0;
            LOOP_CURVA(curva)
                qtdPontos++;

            //Informacoes da curva
            fwrite(&qtdPontos, sizeof(int), 1, arq);
            fwrite(&(curva->regiaoEsq->numero), sizeof(int), 1, arq);
            fwrite(&(curva->regiaoDir->numero), sizeof(int), 1, arq);
            fwrite(&(curva->regiaoEsq->localRegiao), sizeof(TIPO_LOCAL_REGIAO), 1, arq);
            fwrite(&(curva->regiaoDir->numero), sizeof(TIPO_LOCAL_REGIAO), 1, arq);

            LOOP_CURVA(curva) {
                ponto = loop_stuff.ponto;
                fwrite(&(ponto->x), sizeof(double), 1, arq);
                fwrite(&(ponto->y), sizeof(double), 1, arq);
                fwrite(&(ponto->fixo), sizeof(int), 1, arq);
                fwrite(&(ponto->outflow), sizeof(int), 1, arq);
            }
        }
    }

    // Celulas
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fwrite(&(M.celulas[i][j]), sizeof(TIPO_CELULA), 1, arq);
        }
    }

    // TIPO U
    for( i=0; i<=M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fwrite(&(M.pontosU[i][j].tipo), sizeof(TIPO_CELULA), 1, arq);
        }
    }
    // Tipo V
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<=M.Ny; j++ ) {
            fwrite(&(M.pontosV[i][j].tipo), sizeof(TIPO_CELULA), 1, arq);
        }
    }
    // Tipo P
    for( i=0; i<M.Nx; i++ ) {
        for( j=0; j<M.Ny; j++ ) {
            fwrite(&(M.pontosP[i][j].tipo), sizeof(TIPO_CELULA), 1, arq);
        }
    }


    qtdPontos = 0;
    for( celulaSurf=ListaSurf.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox )
        qtdPontos++;
    // ListaSurf qtdPontos
    fwrite(&qtdPontos, sizeof(int), 1, arq);
    for( celulaSurf=ListaSurf.prim; celulaSurf!=NULL; celulaSurf=celulaSurf->prox ) {
        fwrite(&(celulaSurf->i), sizeof(int), 1, arq);
        fwrite(&(celulaSurf->j), sizeof(int), 1, arq);
        fwrite(&(celulaSurf->passoCelula), sizeof(int), 1, arq);
    }

    // Apenas no faraday
//    fwrite(hVelho, sizeof(double), (M.Nx+1), arq);


    fclose(arq);
    return;
}

void CarregarEstado(MALHA *M, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, LISTA_CELULAS_SURFACE *ListaSurf, int *n)
{
//    FILE *arq;
//    char nomeArq[200];
//    int i, j, idRegiaoEsq, idRegiaoDir, localEsq, localDir, qtdRegioes, fixo;
//    int qtdCurvas, qtdPontos, idCelula, passoCelula;
//    double x, y;
//    PONTO *ponto, *proxPonto, *pontoAnt, *primeiro;
//    CURVA *curva, *curva2, *proxCurva;
//    REGIAO *regioes[100], *regiao, *regiaoEsq, *regiaoDir;
//    CELULA_SURFACE *celula, *celulaAnt;
//
//
//
//    snprintf(nomeArq, 200, "debugs/EstadoPasso.txt");
//    arq = fopen(nomeArq, "rt");
//
//    fscanf(arq, "Passo %d\n", n);
//    fscanf(arq, "veloc U\n");
//    for( i=0; i<=M->Nx; i++ ) {
//        for( j=0; j<M->Ny; j++ ) {
//            fscanf(arq, "%lf ", &(U[i][j]));
//        }
//    }
//    fscanf(arq, "\nveloc V\n");
//    for( i=0; i<M->Nx; i++ ) {
//        for( j=0; j<=M->Ny; j++ ) {
//            fscanf(arq, "%lf ", &(V[i][j]));
//        }
//    }
//    fscanf(arq, "Pressao\n");
//    for( i=0; i<M->Nx; i++ ) {
//        for( j=0; j<M->Ny; j++ ) {
//            fscanf(arq, "%lf ", &(P[i][j]));
//        }
//    }
//    fscanf(arq, "Txx\n");
//    for( i=0; i<M->Nx; i++ ) {
//        for( j=0; j<M->Ny; j++ ) {
//            fscanf(arq, "%lf ", &(Txx[i][j]));
//        }
//    }
//    fscanf(arq, "Txy\n");
//    for( i=0; i<M->Nx; i++ ) {
//        for( j=0; j<M->Ny; j++ ) {
//            fscanf(arq, "%lf ", &(Txy[i][j]));
//        }
//    }
//    fscanf(arq, "Tyy\n");
//    for( i=0; i<M->Nx; i++ ) {
//        for( j=0; j<M->Ny; j++ ) {
//            fscanf(arq, "%lf ", &(Tyy[i][j]));
//        }
//    }
//
//    //Liberando as regioes da interface
//    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//
//        //Regiao da esquerda primeiro
//        regiao = curva->regiaoEsq;
//        //Coloca NULL em todos os ponteiros (de todas curvas) que apontam pra essa regiao
//        if( regiao!=NULL ) {
//            for( curva2=M->interface.curvaInicial; curva2!=NULL; curva2=curva2->proxCurva ) {
//                if( curva2->regiaoEsq==regiao )
//                    curva2->regiaoEsq = NULL;
//                if( curva2->regiaoDir==regiao )
//                    curva2->regiaoDir = NULL;
//            }
//
//            free(regiao);
//        }
//
//        //Regiao da direita agora
//        regiao = curva->regiaoDir;
//        //Coloca NULL em todos os ponteiros (de todas curvas) que apontam pra essa regiao
//        if( regiao!=NULL ) {
//            for( curva2=M->interface.curvaInicial; curva2!=NULL; curva2=curva2->proxCurva ) {
//                if( curva2->regiaoEsq==regiao )
//                    curva2->regiaoEsq = NULL;
//                if( curva2->regiaoDir==regiao )
//                    curva2->regiaoDir = NULL;
//            }
//
//            free(regiao);
//        }
//    }
//
//
//
//    //Deletando os pontos da interface atual e as curvas tb
//    curva = M->interface.curvaInicial;
//    while( curva!=NULL ) {
//        proxCurva = curva->proxCurva;
//
//        ponto = curva->nodeInicial;
//        while(ponto!=NULL) {
//            proxPonto = ponto->prox;
//            free(ponto);
//            ponto = proxPonto;
//        }
//        curva->nodeInicial=NULL;
//
//        free(curva);
//        curva = proxCurva;
//    }
//
//    M->interface.curvaInicial = NULL;
//    M->interface.qtdCurvas = 0;
//
//
//    //Temporario
//    //curva = M->interface.curvaInicial;
//
//    qtdRegioes = 0;
//    fscanf(arq, "\nInterface qtdCurvas=%d", &qtdCurvas);
//    for( i=0; i<qtdCurvas; i++ ) {
//
//        //Alocando memoria pra essa curva
//        curva = (CURVA *)malloc( sizeof(CURVA) );
//        curva->proxCurva = NULL;
//        curva->curvaAnt = NULL;
//
//        fscanf(arq, "\nCURVA qtdPontos=%d regiaoEsq=%d regiaoDir=%d localEsq=%d localDir=%d\n", &qtdPontos, &idRegiaoEsq, &idRegiaoDir, &localEsq, &localDir);
//
//        //procura se essa regiao ja foi alocada
//        regiaoEsq = regiaoDir = NULL;
//        for( j=0; j<qtdRegioes; j++ ) {
//            if( regioes[j]->numero==idRegiaoEsq )
//                regiaoEsq = regioes[j];
//            else if( regioes[j]->numero==idRegiaoDir )
//                regiaoDir = regioes[j];
//        }
//
//        //Se uma das duas regioes ainda nao existirem, vamos alocar ela
//        if( regiaoEsq == NULL ) {
//            qtdRegioes++;
//            regiao = (REGIAO *)malloc(sizeof(REGIAO));
//            regiao->localRegiao = localEsq;
//            regiao->numero = idRegiaoEsq;
//            regiao->status = FISICA;
//            regioes[qtdRegioes-1] = regiao;
//            curva->regiaoEsq = regiao;
//        }
//        else
//            curva->regiaoEsq = regiaoEsq;
//        if( regiaoDir == NULL ) {
//            qtdRegioes++;
//            regiao = (REGIAO *)malloc(sizeof(REGIAO));
//            regiao->localRegiao = localDir;
//            regiao->numero = idRegiaoDir;
//            regiao->status = FISICA;
//            regioes[qtdRegioes-1] = regiao;
//            curva->regiaoDir = regiao;
//        }
//        else
//            curva->regiaoDir = regiaoDir;
//
//
//
//        pontoAnt = NULL;
//        primeiro = NULL;
//        //while( 3==fscanf(arq, "[%lf %lf %d] ", &x, &y, &fixo) ) {
//        for( j=0; j<qtdPontos; j++ ) {
//            fscanf(arq, "[%lf %lf %d] ", &x, &y, &fixo);
//
//            ponto = (PONTO *)malloc(sizeof(PONTO));
//            ponto->x = x;
//            ponto->y = y;
//            ponto->fixo = fixo;
//            ponto->isNode = 0;
//            ponto->isNH = 0;
//            ponto->ant = pontoAnt;
//            ponto->prox = NULL;
//
//            if( pontoAnt!=NULL )
//                pontoAnt->prox = ponto;
//
//            if( primeiro==NULL )
//                primeiro = ponto;
//
//            pontoAnt = ponto;
//        }
//
//        primeiro->pontoAntCurva = ponto->ant;
//        ponto->pontoProxCurva = primeiro->prox;
//        curva->nodeInicial = primeiro;
//
//        primeiro->isNode = 1;
//        ponto->isNode = 1;
//
//        AdicionaCurvaNaInterface(&(M->interface), curva);
//        //curva = curva->proxCurva;
//    }
//
//    fscanf(arq, "\nCelulas\n");
//    for( i=0; i<M->Nx; i++ ) {
//        for( j=0; j<M->Ny; j++ ) {
//            fscanf(arq, "%u ", &(M->celulas[i][j]));
//        }
//    }
//
//    fscanf(arq, "\nTipo U\n");
//    for( i=0; i<=M->Nx; i++ ) {
//        for( j=0; j<M->Ny; j++ ) {
//            fscanf(arq, "%u ", &(M->pontosU[i][j].tipo));
//        }
//    }
//    fscanf(arq, "\nTipo V\n");
//    for( i=0; i<M->Nx; i++ ) {
//        for( j=0; j<=M->Ny; j++ ) {
//            fscanf(arq, "%u ", &(M->pontosV[i][j].tipo));
//        }
//    }
//
//    fscanf(arq, "\nTipo P\n");
//    for( i=0; i<M->Nx; i++ ) {
//        for( j=0; j<M->Ny; j++ ) {
//            fscanf(arq, "%u ", &(M->pontosP[i][j].tipo));
//        }
//    }
//
//    fscanf(arq, "\nListaSurf %d\n", &qtdPontos);
//    celulaAnt = NULL;
//    celula = NULL;
//    for( idCelula=0; idCelula<qtdPontos; idCelula++ ) {
//        fscanf(arq, "[%d %d %d] ", &i, &j, &passoCelula);
//
//        celula = (CELULA_SURFACE *)malloc(sizeof(CELULA_SURFACE));
//        celula->i = i;
//        celula->j = j;
//        celula->passoCelula = passoCelula;
//        celula->ant = celulaAnt;
//
//        if( celulaAnt!=NULL )
//            celulaAnt->prox = celula;
//
//        if( ListaSurf->prim==NULL )
//            ListaSurf->prim = celula;
//
//        celulaAnt = celula;
//    }
//
//    if( celula!=NULL )
//        celula->prox = NULL;
//    ListaSurf->ult = celula;
//
//
//    fclose(arq);
//    return;
}

void CarregarEstadoBinario(MALHA *M, double **U, double **V, double **P, double **Txx, double **Txy, double **Tyy, double **Lambda, double **Mu, double **Nu, LISTA_CELULAS_SURFACE *ListaSurf, int *n)
{
    FILE *arq;
    char nomeArq[900];
    int i, j, idRegiaoEsq, idRegiaoDir, localEsq, localDir, qtdRegioes, fixo, outflow;
    int qtdCurvas, qtdPontos, idCelula, passoCelula;
    int checkRead;
    double x, y;
    PONTO *ponto, *proxPonto, *pontoAnt, *primeiro;
    CURVA *curva, *curva2, *proxCurva;
    REGIAO *regioes[100], *regiao, *regiaoEsq, *regiaoDir;
    CELULA_SURFACE *celula, *celulaAnt;


    // sprintf(nomeArq, "ArquivosPlot/VTK/%s/EstadoPasso.txt", M->nomeArquivo);
    sprintf(nomeArq, "debugs/EstadoPasso.txt");
    arq = fopen(nomeArq, "rt");

    // Passo
    checkRead = fread(n, sizeof(int), 1, arq);
    Check_Read(checkRead, 1);
    // Veloc U
    for( i=0; i<=M->Nx; i++ ) {
        checkRead = fread(U[i], sizeof(double), M->Ny, arq);
        Check_Read(checkRead, M->Ny);
    }
    // Veloc V
    for( i=0; i<M->Nx; i++ ) {
        checkRead = fread(V[i], sizeof(double), M->Ny+1, arq);
        Check_Read(checkRead, M->Ny+1);
    }
    // Pressao
    for( i=0; i<M->Nx; i++ ) {
        checkRead = fread(P[i], sizeof(double), M->Ny, arq);
        Check_Read(checkRead, M->Ny);
    }
    // Txx
    for( i=0; i<M->Nx; i++ ) {
        checkRead = fread(Txx[i], sizeof(double), M->Ny, arq);
        Check_Read(checkRead, M->Ny);
    }
    // Txy
    for( i=0; i<M->Nx; i++ ) {
        checkRead = fread(Txy[i], sizeof(double), M->Ny, arq);
        Check_Read(checkRead, M->Ny);
    }
    // Tyy
    for( i=0; i<M->Nx; i++ ) {
        checkRead = fread(Tyy[i], sizeof(double), M->Ny, arq);
        Check_Read(checkRead, M->Ny);
    }
     // Lambda
    for( i=0; i<M->Nx; i++ ) {
        checkRead = fread(Lambda[i], sizeof(double), M->Ny, arq);
        Check_Read(checkRead, M->Ny);
    }
    // Mu
    for( i=0; i<M->Nx; i++ ) {
        checkRead = fread(Mu[i], sizeof(double), M->Ny, arq);
        Check_Read(checkRead, M->Ny);
    }
    // Nu
    for( i=0; i<M->Nx; i++ ) {
        checkRead = fread(Nu[i], sizeof(double), M->Ny, arq);
        Check_Read(checkRead, M->Ny);
    }

    //Liberando as regioes da interface
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



    //Deletando os pontos da interface atual e as curvas tb
    curva = M->interface.curvaInicial;
    while( curva!=NULL ) {
        proxCurva = curva->proxCurva;

        ponto = curva->pontoInicial;
        do {
            proxPonto = ponto->prox;
            free(ponto);
            ponto = proxPonto;
        } while(ponto!=curva->pontoInicial);
        curva->pontoInicial=NULL;

        free(curva);
        curva = proxCurva;
    }

    M->interface.curvaInicial = NULL;
    M->interface.qtdCurvas = 0;


    //Temporario
    //curva = M->interface.curvaInicial;

    qtdRegioes = 0;
    checkRead = fread(&qtdCurvas, sizeof(int), 1, arq);
    Check_Read(checkRead, 1);
    for( i=0; i<qtdCurvas; i++ ) {

        //Alocando memoria pra essa curva
        curva = (CURVA *)malloc( sizeof(CURVA) );
        curva->proxCurva = NULL;
        curva->curvaAnt = NULL;


        checkRead = fread(&qtdPontos, sizeof(int), 1, arq); Check_Read(checkRead, 1);
        checkRead = fread(&idRegiaoEsq, sizeof(int), 1, arq); Check_Read(checkRead, 1);
        checkRead = fread(&idRegiaoDir, sizeof(int), 1, arq); Check_Read(checkRead, 1);
        checkRead = fread(&localEsq, sizeof(TIPO_LOCAL_REGIAO), 1, arq); Check_Read(checkRead, 1);
        checkRead = fread(&localDir, sizeof(TIPO_LOCAL_REGIAO), 1, arq); Check_Read(checkRead, 1);

        //procura se essa regiao ja foi alocada
        regiaoEsq = regiaoDir = NULL;
        for( j=0; j<qtdRegioes; j++ ) {
            if( regioes[j]->numero==idRegiaoEsq )
                regiaoEsq = regioes[j];
            else if( regioes[j]->numero==idRegiaoDir )
                regiaoDir = regioes[j];
        }

        //Se uma das duas regioes ainda nao existirem, vamos alocar ela
        if( regiaoEsq == NULL ) {
            qtdRegioes++;
            regiao = (REGIAO *)malloc(sizeof(REGIAO));
            regiao->localRegiao = localEsq;
            regiao->numero = idRegiaoEsq;
            regiao->status = FISICA;
            regioes[qtdRegioes-1] = regiao;
            curva->regiaoEsq = regiao;
        }
        else
            curva->regiaoEsq = regiaoEsq;
        if( regiaoDir == NULL ) {
            qtdRegioes++;
            regiao = (REGIAO *)malloc(sizeof(REGIAO));
            regiao->localRegiao = localDir;
            regiao->numero = idRegiaoDir;
            regiao->status = FISICA;
            regioes[qtdRegioes-1] = regiao;
            curva->regiaoDir = regiao;
        }
        else
            curva->regiaoDir = regiaoDir;



        pontoAnt = NULL;
        primeiro = NULL;
        //while( 3==fscanf(arq, "[%lf %lf %d] ", &x, &y, &fixo) ) {
        for( j=0; j<qtdPontos; j++ ) {

            checkRead = fread(&x, sizeof(double), 1, arq); Check_Read(checkRead, 1);
            checkRead = fread(&y, sizeof(double), 1, arq); Check_Read(checkRead, 1);
            checkRead = fread(&fixo, sizeof(int), 1, arq); Check_Read(checkRead, 1);
            checkRead = fread(&outflow, sizeof(int), 1, arq); Check_Read(checkRead, 1);

            ponto = (PONTO *)malloc(sizeof(PONTO));
            ponto->x = x;
            ponto->y = y;
            ponto->fixo = fixo;
            ponto->outflow = outflow;
            ponto->isNode = 0;
            ponto->isNH = 0;
            ponto->ant = pontoAnt;
            ponto->prox = NULL;

            if( pontoAnt!=NULL )
                pontoAnt->prox = ponto;

            if( primeiro==NULL )
                primeiro = ponto;

            pontoAnt = ponto;
        }
        curva->pontoInicial = primeiro;

        primeiro->ant = pontoAnt;
        ponto->prox = primeiro;

        primeiro->isNode = 1;
        ponto->isNode = 1;

        AdicionaCurvaNaInterface(&(M->interface), curva);
        //curva = curva->proxCurva;
    }

    // Celulas
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            checkRead = fread(&(M->celulas[i][j]), sizeof(TIPO_CELULA), 1, arq); Check_Read(checkRead, 1);
        }
    }

    // TIPO U
    for( i=0; i<=M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            checkRead = fread(&(M->pontosU[i][j].tipo), sizeof(TIPO_CELULA), 1, arq); Check_Read(checkRead, 1);
        }
    }
    // Tipo V
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<=M->Ny; j++ ) {
            checkRead = fread(&(M->pontosV[i][j].tipo), sizeof(TIPO_CELULA), 1, arq); Check_Read(checkRead, 1);
        }
    }
    // Tipo P
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            checkRead = fread(&(M->pontosP[i][j].tipo), sizeof(TIPO_CELULA), 1, arq); Check_Read(checkRead, 1);
        }
    }


    checkRead = fread(&qtdPontos, sizeof(int), 1, arq); Check_Read(checkRead, 1);
    celulaAnt = NULL;
    celula = NULL;
    for( idCelula=0; idCelula<qtdPontos; idCelula++ ) {
        checkRead = fread(&(i), sizeof(int), 1, arq); Check_Read(checkRead, 1);
        checkRead = fread(&(j), sizeof(int), 1, arq); Check_Read(checkRead, 1);
        checkRead = fread(&(passoCelula), sizeof(int), 1, arq); Check_Read(checkRead, 1);


        celula = (CELULA_SURFACE *)malloc(sizeof(CELULA_SURFACE));
        celula->i = i;
        celula->j = j;
        celula->passoCelula = passoCelula;
        celula->ant = celulaAnt;

        if( celulaAnt!=NULL )
            celulaAnt->prox = celula;

        if( ListaSurf->prim==NULL )
            ListaSurf->prim = celula;

        celulaAnt = celula;
    }

    if( celula!=NULL )
        celula->prox = NULL;
    ListaSurf->ult = celula;




    /// Apenas pro caso faraday
//    if( hVelho == NULL )
//        hVelho = (double *)malloc( (M->Nx+1)*sizeof(double) );
//    checkRead = fread(hVelho, sizeof(double), M->Nx+1, arq); Check_Read(checkRead, M->Nx+1);


    fclose(arq);
    return;
}
















///CRITERIOS PARA CLASSIFICACAO DE CELULAS
int LADO_VERTICAL_DIREITA_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    criterioEsq = (i==0) || (M->celulas[i-1][j]!=EMPTY);
    criterioDir = (i==M->Nx-1) || (M->celulas[i+1][j]==EMPTY);
    criterioCima = (j==M->Ny-1) || (M->celulas[i][j+1]!=EMPTY) || (M->pontosV[i][j+1].tipo==BOUNDARY);
    criterioBaixo = (j==0) || (M->celulas[i][j-1]!=EMPTY) || (M->pontosV[i][j].tipo==BOUNDARY);

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int LADO_VERTICAL_ESQUERDA_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    criterioEsq = (i==0) || (M->celulas[i-1][j]==EMPTY);
    criterioDir = (i==M->Nx-1) || (M->celulas[i+1][j]!=EMPTY);
    criterioCima = (j==M->Ny-1) || (M->celulas[i][j+1]!=EMPTY) || (M->pontosV[i][j+1].tipo==BOUNDARY);
    criterioBaixo = (j==0) || (M->celulas[i][j-1]!=EMPTY) || (M->pontosV[i][j].tipo==BOUNDARY);

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int LADO_HORIZONTAL_CIMA_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    criterioEsq = (i==0) || (M->celulas[i-1][j]!=EMPTY) || (M->pontosU[i][j].tipo==BOUNDARY);
    criterioDir = (i==M->Nx-1) || (M->celulas[i+1][j]!=EMPTY) || (M->pontosU[i+1][j].tipo==BOUNDARY);
    criterioCima = (j==M->Ny-1) || (M->celulas[i][j+1]==EMPTY);
    criterioBaixo = (j==0) || (M->celulas[i][j-1]!=EMPTY) || (M->pontosV[i][j].tipo==BOUNDARY);

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int LADO_HORIZONTAL_BAIXO_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    criterioEsq = (i==0) || (M->celulas[i-1][j]!=EMPTY) || (M->pontosU[i][j].tipo==BOUNDARY);
    criterioDir = (i==M->Nx-1) || (M->celulas[i+1][j]!=EMPTY) || (M->pontosU[i+1][j].tipo==BOUNDARY);
    criterioCima = (j==M->Ny-1) || (M->celulas[i][j+1]!=EMPTY);
    criterioBaixo = (j==0) || (M->celulas[i][j-1]==EMPTY);

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int QUINA_CIMA_DIREITA_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    if( M->pontosV[i][j+1].tipo==BOUNDARY )
        return 0;

    criterioEsq = (i==0) || (M->pontosU[i][j].tipo==BOUNDARY) || (M->celulas[(i)-1][j]!=EMPTY);
    criterioDir = (i==M->Nx-1) || (M->celulas[(i)+1][j]==EMPTY);
    criterioCima = (j!=M->Ny-1) && (M->celulas[(i)][(j)+1]==EMPTY);
    criterioBaixo = (j==0) || (M->pontosV[i][j].tipo==BOUNDARY) || (M->celulas[(i)][(j)-1]!=EMPTY);

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int QUINA_CIMA_ESQUERDA_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    if( M->pontosV[i][j+1].tipo==BOUNDARY )
        return 0;

    criterioEsq = (i==0) || (M->celulas[(i)-1][j]==EMPTY);
    criterioDir = (i==M->Nx-1) || (M->pontosU[i+1][j].tipo==BOUNDARY) || (M->celulas[(i)+1][j]!=EMPTY);
    criterioCima = (j==M->Ny-1) || (M->celulas[(i)][(j)+1]==EMPTY);
    criterioBaixo = (j==0) || (M->pontosV[i][j].tipo==BOUNDARY) || (M->celulas[i][(j)-1]!=EMPTY);

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int QUINA_BAIXO_DIREITA_func(int i, int j, MALHA *M)
{
    if( M->pontosV[i][j].tipo==BOUNDARY )
        return 0;

    return( (i==0 || M->celulas[(i)-1][j]!=EMPTY) && (i==M->Nx-1 || M->celulas[(i)+1][j]==EMPTY) && (j==M->Ny-1 || M->celulas[i][(j)+1]!=EMPTY) && (j==0 || M->celulas[i][(j)-1]==EMPTY) );
}

int QUINA_BAIXO_ESQUERDA_func(int i,int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    if( M->pontosV[i][j].tipo==BOUNDARY )
        return 0;

    criterioEsq = (i==0) || (M->celulas[(i)-1][j]==EMPTY);
    criterioDir = (i==M->Nx-1) || (M->celulas[(i)+1][j]!=EMPTY);
    criterioCima = (j==M->Ny-1) || (M->celulas[i][(j)+1]!=EMPTY);
    criterioBaixo = (j==0) || (M->celulas[i][(j)-1]==EMPTY);

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int DEGENERADO_CIMA_BAIXO_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    if( i==M->Nx-1 || j==0 || j==M->Ny-1 )
        return 0;

    criterioEsq = (i==0) || M->celulas[(i)-1][j]!=EMPTY;
    criterioDir = M->celulas[(i)+1][j]!=EMPTY;
    criterioCima = M->celulas[i][(j)+1]==EMPTY;
    criterioBaixo = M->celulas[i][(j)-1]==EMPTY;

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int DEGENERADO_ESQ_DIR_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    if( i==0 || i==M->Nx-1 || j==0 || j==M->Ny-1 )
        return 0;

    criterioEsq = M->celulas[(i)-1][j]==EMPTY;
    criterioDir = M->celulas[(i)+1][j]==EMPTY;
    criterioCima = M->celulas[i][(j)+1]!=EMPTY;
    criterioBaixo = M->celulas[i][(j)-1]!=EMPTY;

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

int DEGENERADO_ESQ_CIMA_BAIXO_func(int i, int j, MALHA *M)
{
    static int criterioEsq, criterioDir, criterioCima, criterioBaixo;

    if( i==0 || i==M->Nx-1 || j==0 || j==M->Ny-1 )
        return 0;

    criterioEsq = M->celulas[(i)-1][j]==EMPTY;
    criterioDir = M->celulas[(i)+1][j]!=EMPTY;
    criterioCima = M->celulas[i][(j)+1]==EMPTY;
    criterioBaixo = M->celulas[i][(j)-1]==EMPTY;

    return(criterioEsq && criterioDir && criterioCima && criterioBaixo);
}

