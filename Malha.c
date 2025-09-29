#include "Malha.h"

MALHA LerModeloDoArquivo(const char *NomeArquivo)
{
    FILE *arq, *arqOriginal;
    MALHA malha;
    char nomeArqJunk[500] = "";
    char tipoModelo[500] = "";
    char tipoBloco[100], stringLida[100];
    char tipoBoundary[100] = "", tipoDirecao[100] = "";
    char perfilInflow[100] = "", movimento[100] = "";
    char tensaoSup[100] = "";
    char tipoCoord[200] = "", lixo[200], tipoStretchingX[100], tipoStretchingY[100];
    double valorDirichlet, valorDirichlet2, posicao, inicio, fim;
    double *verticesStretchingX = NULL, *paramStretchingX = NULL;
    double *verticesStretchingY = NULL, *paramStretchingY = NULL;
    int qtdVerticesX=0, qtdVerticesY=0, *Nx = NULL, *Ny = NULL;
    int id, iMin = -1, iMax = -1, jMin = -1, jMax = -1, numFaces, *faces;
    int i, j, testeScan;
    double x, y;
    REGIAO *regiaoFundo;
    int qtdRegioes;
    BLOCO_PARAVIEW *blocoParaview = NULL;

    //Criando uma regiao de fundo pra estrutura interface
    malha.interface.curvaInicial = NULL;
    malha.interface.qtdCurvas = 0;
    regiaoFundo = (REGIAO *)malloc( sizeof(REGIAO) );
    regiaoFundo->localRegiao = EXTERNA;
    regiaoFundo->numero = 0;
    regiaoFundo->status = FISICA;
    qtdRegioes = 1;

    /// === Abrindo o arquivo. Criando um arquivo temporario todo maiusculo. Pra ficar case-insensitive
    arqOriginal = fopen(NomeArquivo, "rt");
    if( arqOriginal==NULL ) DeuProblema("\n\nERRO: O nome de arquivo fornecido (%s) nao foi encontrado.\n\n", NomeArquivo);
    sprintf(nomeArqJunk, "%s_JUNK", NomeArquivo);
    arq = fopen(nomeArqJunk, "wt");
    ArquivoMaiusculoTemporario(arqOriginal, arq);
    fclose(arq);
    arq = fopen(nomeArqJunk, "rt");
    rewind(arqOriginal);


    // Lendo a linha output do arquivo original (pra nao ficar tudo maiusculo)
    testeScan = fscanf(arq, "%s VTK=[%d %d] BACKUP=%d PRINT=%d PASTA=%s\n", tipoCoord, &(malha.intervaloVTK_prop), &(malha.intervaloVTK_surf), &(malha.intervaloEstado), &(malha.intervaloPrint), malha.nomeArquivo);
    if( testeScan!=6 || (strcmp(tipoCoord, "OUTPUT")) ) DeuProblema("\n\nERRO: problema na leitura da linha OUTPUT. %d\n\n", testeScan);
    testeScan = fscanf(arqOriginal, "%s %s %s %s %s %s\n", lixo, lixo, lixo, lixo, lixo, lixo);

    //Removendo os caracteres "pasta=" da string
    for( i=6; lixo[i]!='\0'; i++ )
        malha.nomeArquivo[i-6] = lixo[i];
    malha.nomeArquivo[i-6] = '\0';

    //Le a primeira linha do arq pra ver se eh CARTESIANO ou AXISSIMETRICO_CILINDRICO
    testeScan = fscanf(arq, "%s\n", tipoCoord);
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do sistema de coordenadas (CARTESIANO OU AXI_CILINDRICO).\n\n");
    if( !strcmp(tipoCoord, "CARTESIANO") )
        malha.tipoCoord = CARTESIANO;
    else if( !strcmp(tipoCoord, "AXI_CILINDRICO") )
        malha.tipoCoord = AXI_CILINDRICO;
    else
        DeuProblema("\n\nERRO: A primeira linha %s do arquivo modelo deve ser CARTESIANO ou AXI_CILINDRICO.\n\n", tipoCoord);

    //Lendo Nt
    testeScan = fscanf(arq, "NT %d\n", &(malha.Nt));
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura de Nt.\n\n");

    testeScan = fscanf(arq, "XMIN %lf\n", &(malha.xMin));
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura de xMin.\n\n");

    testeScan = fscanf(arq, "XMAX %lf\n", &(malha.xMax));
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura de xMax.\n\n");

    testeScan = fscanf(arq, "YMIN %lf\n", &(malha.yMin));
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura de yMin.\n\n");

    testeScan = fscanf(arq, "YMAX %lf\n", &(malha.yMax));
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura de yMax.\n\n");

    testeScan = fscanf(arq, "TMAX %lf\n", &(malha.tMax));
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura de tMax.\n\n");

    testeScan = fscanf(arq, "RE %lf\n", &(malha.Re));
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura de Re.\n\n");

    testeScan = fscanf(arq, "FR %lf GRAV=[%lf %lf]\n", &(malha.Fr), &(malha.gravX), &(malha.gravY));
    if( testeScan!=3 ) DeuProblema("\n\nERRO: Problema na leitura de Fr e grav.\n\n");

    testeScan = fscanf(arq, "TIPO_MODELO %s\n", tipoModelo);
    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura dos tipo_modelo.\n\n");

    if( !strcmp(tipoModelo, "VE") ) {
        testeScan = fscanf(arq, "PARAMETROS_VE WI=%lf BETA=%lf ALPHA=%lf EPSILON=%lf FORMULACAO=%s\n", &(malha.We), &(malha.beta), &(malha.alpha), &(malha.epsilon), tipoCoord);
        if( testeScan!=5 ) DeuProblema("\n\nERRO: Problema na leitura dos param_visco.\n\n");

        malha.eta_inf = 1e+10; // Garbage value. Will not be used for the VE model
        malha.Bi = 0.0;

        malha.tipo_modelo = MODELO_VE;
    }
    else if( !strcmp(tipoModelo, "EVPT") ) {
        testeScan = fscanf(arq, "PARAMETROS_EVPT WI=%lf ETA_INF=%lf ETA_0=%lf YIELD_STRESS=%lf K=%lf N=%lf T_EQ=%lf M=%lf\n", &(malha.We), &(malha.eta_inf), &(malha.eta_0), &(malha.yield_stress), &(malha.K), &(malha.power_n), &(malha.time_eq), &(malha.evpt_m));
        if( testeScan!=8 ) DeuProblema("\n\nERRO: Problema na leitura dos parametros_evpt.\n\n");
        strcpy(tipoCoord, "CSF");
        malha.beta = 1e+10; // garbage value, since we wont use beta in the EPVT model
        malha.epsilon = 0.0; // Just keeping it as 0, so it automatically removes the PTT part from equations
        malha.alpha = 0.0; // Just keeping it as 0, so it automatically removes the Giesekus part from equations

        malha.tipo_modelo = MODELO_EVPT;
    }
    else if( !strcmp(tipoModelo, "EVP_SARAMITO") ) {
        testeScan = fscanf(arq, "PARAMETROS_EVP WI=%lf BI=%lf BETA=%lf\n", &(malha.We), &(malha.Bi), &(malha.beta));
        if( testeScan!=3 ) DeuProblema("\n\nERRO: Problema na leitura dos param_visco. %d\n\n", testeScan);
        strcpy(tipoCoord, "CSF");
        malha.epsilon = 0.0;
        malha.alpha = 0.0;

        malha.tipo_modelo = MODELO_EVP_SARAMITO;
    }
    else if( !strcmp(tipoModelo, "VP") ) {
        testeScan = fscanf(arq, "PARAMETROS_VP BI=%lf\n", &(malha.Bi));
        if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura dos param_visco. %d\n\n", testeScan);
        strcpy(tipoCoord, "CSF");
        malha.epsilon = 0.0;
        malha.alpha = 0.0;
        malha.beta = 1.0;
        malha.We = 0.0;

        malha.tipo_modelo = MODELO_VP;
    }
    else DeuProblema("\n\nERRO: O tipo de modelo fornecido nao eh um dos reconhecidos (VE ou EVPT).\n");

    testeScan = fscanf(arq, "TOLNSF %lf %lf %d\n", &(malha.tolNSF), &(malha.tolNSFconversao), &(malha.grauHibrido));
    if( testeScan!=3 ) DeuProblema("\n\nERRO: Problema na leitura de tolNSF.\n\n");

    testeScan = fscanf(arq, "TENSAO_SUPERFICIAL LIGADO=%s WEBER=%lf TSUR=%d\n", tensaoSup, &(malha.weber), &(malha.freqTSUR));
    if( testeScan!=3 ) DeuProblema("\n\nERRO: Problema na leitura de tensao_superficial.\n\n");
    if( strcmp(tensaoSup, "SIM") && strcmp(tensaoSup, "NAO") )
        DeuProblema("\n\nERRO: Problema na leitura de tensao_superficial. O parametro ligado deve ser SIM ou NAO.\n\n");
    if( !strcmp(tensaoSup, "SIM") )
        malha.tensaoSuperficial = 1;
    else
        malha.tensaoSuperficial = 0;

    testeScan = fscanf(arq, "VEL_INICIAL [%lf %lf]\n", &(malha.velInicialX), &(malha.velInicialY));
    if( testeScan!=2 ) DeuProblema("\n\nERRO: Problema na leitura de VEL_INICIAL %d.\n\n", testeScan);

    if( !strcmp(tipoCoord, "CSF") )
        malha.formulacaoVisco = VISCO_CART;
    else if( !strcmp(tipoCoord, "NSF") )
        malha.formulacaoVisco = VISCO_NSF;
    else if( !strcmp(tipoCoord, "HIB") )
        malha.formulacaoVisco = VISCO_HIBRIDO;
    else if( !strcmp(tipoCoord, "LOG") )
        malha.formulacaoVisco = VISCO_LOG;
    else DeuProblema("\n\nERRO: A FORMULACAO_VISCO deve ser CSF, NSF, HIB ou LOG.\n\n");

    //Lendo quantos vertices tem no stretching do eixo X e eixo Y
    LeituraMalhaStretching(arq, tipoStretchingX, &qtdVerticesX, &verticesStretchingX, &Nx, &paramStretchingX);
    LeituraMalhaStretching(arq, tipoStretchingY, &qtdVerticesY, &verticesStretchingY, &Ny, &paramStretchingY);

    CriaMalha(&malha, tipoStretchingX, qtdVerticesX, verticesStretchingX, Nx, paramStretchingX, tipoStretchingY, qtdVerticesY, verticesStretchingY, Ny, paramStretchingY);

    free(Nx);
    free(Ny);
    if( verticesStretchingX!=NULL )
        free(verticesStretchingX);
    if( verticesStretchingY!=NULL )
        free(verticesStretchingY);
    if( paramStretchingX!=NULL )
        free(paramStretchingX);
    if( paramStretchingY!=NULL )
        free(paramStretchingY);

    malha.dt = malha.tMax / (double)(malha.Nt);

    while( !feof(arq) ) {
        strcpy(movimento, "");

        valorDirichlet = valorDirichlet2 = 0.0;
        testeScan = fscanf(arq, "%s ", stringLida); //Le a primeira string desta linha do arquivo

        if( !strcmp(stringLida, "BOUNDARY") ) //Vamos inserir um bloco na malha
        {
            strcpy(tipoBloco, stringLida); //Lendo o tipo do bloco

            testeScan = fscanf(arq, "ID=%d ", &id);
            if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do ID de uma boundary.\n\n");

            testeScan = fscanf(arq, "DIRECAO=%s ", tipoDirecao); //Direcao normal a este contorno
            if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura da DIRECAO de uma boundary.\n\n");
//
            testeScan = fscanf(arq, "POSICAO=%lf INICIO=%lf FIM=%lf ", &posicao, &inicio, &fim); //Lendo os extremos do bloco
            if( testeScan!=3 ) DeuProblema("\n\nERRO: Problema na leitura da POSICAO, INICIO ou FIM de uma boundary.\n\n");

            testeScan = fscanf(arq, "TIPO=%s ", tipoBoundary); //Tipo deste contorno (noslip, inflow, outflow, ...)
            if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do TIPO de uma boundary.\n\n");

            if( !strcmp(tipoBoundary, "INFLOW") ) { //Lendo o valor pra ser imposto caso for dirichlet {
                testeScan = fscanf(arq, "PERFIL=%s VALORDIRICHLET=%lf\n", perfilInflow, &valorDirichlet);
                if( testeScan!=2 ) DeuProblema("\n\nERRO: Problema na leitura do PERFIL ou VALORDIRICHLET de uma boundary inflow.\n\n");
            }

            if( !strcmp(tipoBoundary, "NOSLIP") ) {
                testeScan = fscanf(arq, "MOVIMENTO=%s ", movimento);
                if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do MOVIMENTO de uma boundary noslip.\n\n");

                if( !strcmp(movimento, "DIRECAO_X") || !strcmp(movimento, "DIRECAO_Y") || !strcmp(movimento, "DIRECAO_THETA") ) {
                    testeScan = fscanf(arq, "VALORDIRICHLET=%lf\n", &valorDirichlet);
                    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do VALORDIRICHLET de uma boundary noslip com movimento.\n\n");
                }
            }
            if( !strcmp(tipoBoundary, "SLIP") ) {
                testeScan = fscanf(arq, "MOVIMENTO=%s ", movimento);
                if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do MOVIMENTO de uma boundary slip.\n\n");

                if( !strcmp(movimento, "DIRECAO_X") || !strcmp(movimento, "DIRECAO_Y") || !strcmp(movimento, "DIRECAO_THETA") ) {
                    testeScan = fscanf(arq, "COEFICIENTE=%lf\n", &valorDirichlet);
                    if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do VALORDIRICHLET de uma boundary noslip com movimento.\n\n");
                }
            }
            if( !strcmp(tipoBoundary, "KNOWN_PRESSURE") ) {
                testeScan = fscanf(arq, "PRESSURE=%lf PRESSURE_GRAD=%lf ", &valorDirichlet, &valorDirichlet2);
                if( testeScan!=2 ) DeuProblema("\n\nERRO: Problema na leitura dos parametros da condicao de contorno KNOWN_PRESSURE.\n\n");
            }

            double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;

            //Encontrando quais celulas vao ser atribuidas a este contorno
            if( !strcmp(tipoDirecao, "VERTICAL") ) {
                strcpy(tipoDirecao, "x"); //DirecaoNormal pra usar na outra funcao

                x1 = x2 = posicao;
                y1 = inicio;
                y2 = fim;

                x = malha.xMin;
                for( i=0; i<=malha.Nx; i++ ) {
                    if( ZERO(x-posicao) )
                        break;

                    if( i!=malha.Nx )
                        x += malha.dx[i];
                }
                iMin = iMax = i;

                if( !ZERO(x-posicao) ) {
                    DeuProblema("\n\nERRO 1: A malha criada nao coincide com alguma das suas BOUNDARY. \
                                \n      Verifique se as BOUNDARY estao certas ou entao ajuste a malha corretamente %lf %lf.\n\n", x, posicao);
                }

                double yMin = 1e+10, yMax = 1e+10;
                y = malha.yMin;
                for( j=0; j<=malha.Ny; j++ ) {
                    if( ZERO(y-inicio) ) {
                        jMin = j;
                        yMin = y;
                    }
                    if( ZERO(y-fim) ) {
                        jMax = j;
                        yMax = y;
                    }

                    if( j!=malha.Ny )
                        y += malha.dy[j];
                }

                if( !ZERO(yMin-inicio) || !ZERO(yMax-fim) ) {
                    DeuProblema("\n\nERRO 2: A malha criada nao coincide com alguma das suas BOUNDARY. \
                                \n      Verifique se as BOUNDARY estao certas ou entao ajuste a malha corretamente.\n\n");
                }

                //gambiarra um pouco :/
                jMax--;
            }
            else if( !strcmp(tipoDirecao, "HORIZONTAL") ) {
                strcpy(tipoDirecao, "y"); //DirecaoNormal pra usar na outra funcao

                y1 = y2 = posicao;
                x1 = inicio;
                x2 = fim;

                y = malha.yMin;
                for( j=0; j<=malha.Ny; j++ ) {
                    if( ZERO(y-posicao) )
                        break;

                    if( j!=malha.Ny )
                        y += malha.dy[j];
                }
                jMin = jMax = j;

                if( !ZERO(y-posicao) ) {
                    DeuProblema("\n\nERRO 3: A malha criada nao coincide com alguma das suas BOUNDARY. \
                                \n      Verifique se as BOUNDARY estao certas ou entao ajuste a malha corretamente.\n\n");
                }

                double xMin = 1e+10, xMax = 1e+10;
                x = malha.xMin;
                for( i=0; i<=malha.Nx; i++ ) {
                    if( ZERO(x-inicio) ) {
                        iMin = i;
                        xMin = x;
                    }
                    if( ZERO(x-fim) ) {
                        iMax = i;
                        xMax = x;
                    }

                    if( i!=malha.Nx )
                        x += malha.dx[i];
                }

                if( !ZERO(xMin-inicio) || !ZERO(xMax-fim) ) {
                    DeuProblema("\n\nERRO 4: A malha criada nao coincide com alguma das suas BOUNDARY. \
                                \n      Verifique se as BOUNDARY estao certas ou entao ajuste a malha corretamente %lf %lf.\n\n", xMin, xMax);
                }

                //gambiarra um pouco :/
                iMax--;
            }
            else DeuProblema("\n\nERRO: A direcao fornecida pra alguma BOUNDARY eh invalida. Opcoes validas: VERTICAL e HORIZONTAL.\n\n");

            InsereBloco(&malha, id, tipoBloco, tipoDirecao, tipoBoundary, valorDirichlet, valorDirichlet2, perfilInflow, movimento,
                        iMin, iMax, jMin, jMax, x1, y1, x2, y2);
        }
        else if( !strcmp(stringLida, "REGIAO") ) {
            testeScan = fscanf(arq, "TIPO=%s ", tipoBoundary);
            if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do TIPO de uma REGIAO.\n\n");
            testeScan = fscanf(arq, "NUMFACES=%d FACES=", &numFaces);
            if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do NUMFACES de uma REGIAO.\n\n");

            faces = (int *)malloc( numFaces*sizeof(int) );
            for( i=0; i<numFaces; i++ ) {
                testeScan = fscanf(arq, "%d ", &faces[i]);
                if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura nos IDs das FACES de uma REGIAO.\n\n");
            }
            PreencheRegiao(malha, numFaces, faces);
            free(faces);
        }
        else if( !strcmp(stringLida, "REGIAO_PAREDE") ) {
            testeScan = fscanf(arq, "NUMFACES=%d FACES=", &numFaces);
            if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do NUMFACES de uma REGIAO.\n\n");

            faces = (int *)malloc( numFaces*sizeof(int) );
            for( i=0; i<numFaces; i++ ) {
                testeScan = fscanf(arq, "%d ", &faces[i]);
                if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura nos IDs das FACES de uma REGIAO.\n\n");
            }

            PreencheRegiaoParede(malha, numFaces, faces);
            free(faces);
        }
        else if( !strcmp(stringLida, "REGIAO_LIVRE") ) {
            REGIAO *novaRegiao;
            novaRegiao = (REGIAO *)malloc( sizeof(REGIAO) );
            novaRegiao->localRegiao = INTERNA;
            novaRegiao->numero = qtdRegioes++;
            novaRegiao->status = FISICA;

            testeScan = fscanf(arq, "TIPO=%s ", tipoBoundary);
            if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do TIPO de uma REGIAO_LIVRE.\n\n");

            if( !strcmp(tipoBoundary, "ELIPSE") ) {
                double centroX, centroY, raioX, raioY;
                double freeVelX, freeVelY;

                testeScan = fscanf(arq, "CENTRO=[%lf %lf] RAIO=[%lf %lf] VEL=[%lf %lf]\n", &centroX, &centroY, &raioX, &raioY, &freeVelX, &freeVelY);
                if( testeScan!=6 ) DeuProblema("\n\nERRO: Problema na leitura do CENTRO ou RAIO de uma elipse (REGIAO_LIVRE).\n\n");

                //teste
//                raioX = raioX - 0.04*(malha.dx[0]);
//                raioY = raioY - 0.04*(malha.dx[0]);

                AdicionaCurvaElipse(&(malha.interface), centroX, centroY, raioX, raioY, novaRegiao, regiaoFundo, 0.0, 200);
                //AdicionaCurvaDoArquivo(&(malha.interface), "debugs/pontosIniciais.FS", novaRegiao, regiaoFundo);
                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);
                InsereBlocoFreeSurface(&malha, TIPO_FREE_ELIPSE, freeVelX, freeVelY, 0.0, 0.0, 0.0, 0.0, raioX, raioY, centroX, centroY);
            }
            else if( !strcmp(tipoBoundary, "ELIPSE_AXI") ) {
                double centroX, centroY, raioX, raioY;
                double freeVelX, freeVelY;

                testeScan = fscanf(arq, "CENTRO=[%lf %lf] RAIO=[%lf %lf] VEL=[%lf %lf]\n", &centroX, &centroY, &raioX, &raioY, &freeVelX, &freeVelY);
                if( testeScan!=6 ) DeuProblema("\n\nERRO: Problema na leitura do CENTRO ou RAIO de uma elipse (REGIAO_LIVRE).\n\n");

                //teste
//                raioX = raioX - 0.04*(malha.dx[0]);
//                raioY = raioY - 0.04*(malha.dx[0]);

                AdicionaCurvaElipseAxissimetrico(&(malha.interface), centroX, centroY, raioX, raioY, novaRegiao, regiaoFundo, 0.5*M_PI-1e-6, 200);

                // DeuProblema("AFHAUHFUA\n");
                //AdicionaCurvaDoArquivo(&(malha.interface), "debugs/pontosIniciais.FS", novaRegiao, regiaoFundo);
                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);
                InsereBlocoFreeSurface(&malha, TIPO_FREE_ELIPSE, freeVelX, freeVelY, 0.0, 0.0, 0.0, 0.0, raioX, raioY, centroX, centroY);
            }
            else if( !strcmp(tipoBoundary, "SEMI_ELIPSE_AXI") ) {
                double centroX, centroY, raioX, raioY;
                double freeVelX, freeVelY;

                testeScan = fscanf(arq, "CENTRO=[%lf %lf] RAIO=[%lf %lf] VEL=[%lf %lf]\n", &centroX, &centroY, &raioX, &raioY, &freeVelX, &freeVelY);
                if( testeScan!=6 ) DeuProblema("\n\nERRO: Problema na leitura do CENTRO ou RAIO de uma elipse (REGIAO_LIVRE).\n\n");

                //teste
//                raioX = raioX - 0.04*(malha.dx[0]);
//                raioY = raioY - 0.04*(malha.dx[0]);

                AdicionaCurvaSemiElipseAxissimetrico(&(malha.interface), centroX, centroY, raioX, raioY, novaRegiao, regiaoFundo, 0.5*M_PI-1e-6, 200);

                //AdicionaCurvaDoArquivo(&(malha.interface), "debugs/pontosIniciais.FS", novaRegiao, regiaoFundo);
                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);
                InsereBlocoFreeSurface(&malha, TIPO_FREE_ELIPSE, freeVelX, freeVelY, 0.0, 0.0, 0.0, 0.0, raioX, raioY, centroX, centroY);
            }
            else if( !strcmp(tipoBoundary, "SEMI_ELIPSE") ) {
                double centroX, centroY, raioX, raioY;

                testeScan = fscanf(arq, "CENTRO=[%lf %lf] RAIO=[%lf %lf]\n", &centroX, &centroY, &raioX, &raioY);
                if( testeScan!=4 ) DeuProblema("\n\nERRO: Problema na leitura do CENTRO ou RAIO de uma elipse (REGIAO_LIVRE).\n\n");

                //teste
//                raioX = raioX - 0.04*(malha.dx[0]);
//                raioY = raioY - 0.04*(malha.dx[0]);

                AdicionaCurvaSemiElipse(&(malha.interface), centroX, centroY, raioX, raioY, novaRegiao, regiaoFundo, 100);
                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);
                // DeuProblema("ADICIONAR O BLOCO FREE-SURFACE DA SEMI_ELIPSE");
            }
            else if( !strcmp(tipoBoundary, "SEMI_ELIPSE_ARCO") ) {
                double R0, contact_angle, translate_x, translate_y;

                testeScan = fscanf(arq, "R0=%lf CONTACT_ANGLE=%lf TRANSLATE=[%lf %lf]\n", &R0, &contact_angle, &translate_x, &translate_y);
                if( testeScan!=4 ) DeuProblema("\n\nERRO: Problema na leitura dos parametros SEMI_ELIPSE_ARCO (REGIAO_LIVRE).\n\n");

                AdicionaCurvaSemiElipseCoalescencia(&(malha.interface), R0, contact_angle, translate_x, translate_y, novaRegiao, regiaoFundo, 100);
                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);
                // InsereBlocoFreeSurface(&malha, TIPO_FREE_SEMI_ARCO, freeVelX, freeVelY, 0.0, 0.0, 0.0, 0.0, raioX, raioY, centroX, centroY);
            }
            else if( !strcmp(tipoBoundary, "SPREADING") ) {
                double R0, h_inf;

                testeScan = fscanf(arq, "R0=%lf H_INF=%lf\n", &R0, &h_inf);
                if( testeScan!=2 ) DeuProblema("\n\nERRO: Problema na leitura dos parametros SPREADING (REGIAO_LIVRE).\n\n");

                AdicionaCurvaMaziSpreading(&(malha.interface), h_inf, R0, 5.0, novaRegiao, regiaoFundo);
                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);
                // InsereBlocoFreeSurface(&malha, TIPO_FREE_SEMI_ARCO, freeVelX, freeVelY, 0.0, 0.0, 0.0, 0.0, raioX, raioY, centroX, centroY);
            }
            else if( !strcmp(tipoBoundary, "RETANGULO") ) {
                double xMinRet, xMaxRet, yMinRet, yMaxRet;
                double freeVelX, freeVelY;

                testeScan = fscanf(arq, "EXTREMOS=[%lf %lf %lf %lf] VEL=[%lf %lf]\n", &xMinRet, &xMaxRet, &yMinRet, &yMaxRet, &freeVelX, &freeVelY);
                if( testeScan!=6 ) DeuProblema("\n\nERRO: Problema na leitura do os EXTREMOS de um retangulo (REGIAO_LIVRE).\n\n");
                AdicionaCurvaRetangulo(&(malha.interface), xMinRet, xMaxRet, yMinRet, yMaxRet, novaRegiao, regiaoFundo, 0.2*malha.minDxDy);
                AdicionaERemoveParticulas(&malha, 0.2*malha.minDxDy, 1);
                InsereBlocoFreeSurface(&malha, TIPO_FREE_RETANGULO, freeVelX, freeVelY, xMinRet, yMinRet, xMaxRet, yMaxRet, 0.0, 0.0, 0.0, 0.0);
            }
            else if( !strcmp(tipoBoundary, "RETANGULO_FIXO") ) {
                double xMinRet, xMaxRet, yMinRet, yMaxRet;

                testeScan = fscanf(arq, "EXTREMOS=[%lf %lf %lf %lf]\n", &xMinRet, &xMaxRet, &yMinRet, &yMaxRet);
                if( testeScan!=4 ) DeuProblema("\n\nERRO: Problema na leitura do os EXTREMOS de um retangulo (REGIAO_LIVRE).\n\n");
                AdicionaCurvaRetangulo(&(malha.interface), xMinRet, xMaxRet, yMinRet, yMaxRet, novaRegiao, regiaoFundo, 0.2*malha.minDxDy);
                AdicionaERemoveParticulas(&malha, 0.2*malha.minDxDy, 1);

                //Faz os pontos ficarem fixos
                PONTO *p;
                LOOP_CURVA(malha.interface.curvaInicial) {
                    p = loop_stuff.ponto;
                    p->fixo = 1;
                }
            }
            else if( !strcmp(tipoBoundary, "CANAL_SUPLIVRE") ) {
                double xMinRet, xMaxRet, yMinRet, yMaxRet;

                testeScan = fscanf(arq, "EXTREMOS=[%lf %lf %lf %lf]\n", &xMinRet, &xMaxRet, &yMinRet, &yMaxRet);
                if( testeScan!=4 ) DeuProblema("\n\nERRO: Problema na leitura do os EXTREMOS de um retangulo (REGIAO_LIVRE).\n\n");
                AdicionaCurvaRetangulo(&(malha.interface), xMinRet, xMaxRet, yMinRet, yMaxRet, novaRegiao, regiaoFundo, 0.2*malha.minDxDy);
                AdicionaERemoveParticulas(&malha, 0.2*malha.minDxDy, 1);

                //Faz os pontos ficarem fixos
                PONTO *p;
                LOOP_CURVA(malha.interface.curvaInicial) {
                    p = loop_stuff.ponto;
                    p->fixo = 1;
                }
            }
            else if( !strcmp(tipoBoundary, "INFLOW") ) {
                double posicao = 0.0, coordMin = 0.0, coordMax = 0.0, comprimento = 0.0;
                char direcao = 'a';
                int id;
                BLOCO *bloco;

                testeScan = fscanf(arq, "ID=%d ", &id);
                testeScan = fscanf(arq, "COMPRIMENTO=%lf\n", &comprimento);
                if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura de um ID para criar um inflow (REGIAO_LIVRE).\n\n");
                //Encontra esse bloco inflow pra pegar as informacoes
                for( bloco=malha.primBloco; bloco!=NULL; bloco=bloco->prox ) {
                    if( bloco->id==id )
                        break;
                }
                if( bloco==NULL )
                    DeuProblema("\n\nPROBLEMA NA CRIACAO DO INFLOW.\n\n");

                if( bloco->tipoU==BOUNDARY && bloco->tipoBoundaryU==INFLOW ) {
                    direcao = 'v'; //Vertical
                    posicao = malha.x[bloco->iMin];
                    coordMin = malha.y[bloco->jMin];
                    coordMax = malha.y[bloco->jMax+1];

                    //usando o sinal da velocidade pra adicionar uma tolerancia na posicao
                    if( bloco->valorDirichlet>0 )
                        posicao += 1e-8;
                    else if( bloco->valorDirichlet<0 )
                        posicao += -1e-8;
                }
                else if( bloco->tipoV==BOUNDARY && bloco->tipoBoundaryV==INFLOW ) {
                    direcao = 'h'; //Vertical
                    posicao = malha.y[bloco->jMin];
                    coordMin = malha.x[bloco->iMin];
                    coordMax = malha.x[bloco->iMax + 1];
                }
                else
                    DeuProblema("\n\nPROBLEMA NA CRIACAO DE UM INFLOW_LIVRE. \n A ID ESPECIFICADA NAO EH INFLOW.\n %d %d\n\n", bloco->tipoU, bloco->tipoBoundaryU);

                AdicionaCurvaInflow(&(malha.interface), posicao, comprimento, coordMin, coordMax, /*(malha.minDxDy)/3.0 */ 1e-5, direcao, novaRegiao, regiaoFundo, 300);
                AdicionaERemoveParticulas(&malha, 0.05*(malha.minDxDy), 1);
            }
            else if( !strcmp(tipoBoundary, "FILAMENTO") ) {
                double Rp = 0.0, R0 = 0.0, Lz = 0.0;

                testeScan = fscanf(arq, "ALTURA=%lf ", &Lz);
                testeScan = fscanf(arq, "RAIO_CENTRO=%lf ", &R0);
                AdicionaCurvaFilamento(&(malha.interface), Rp, R0, Lz, novaRegiao, regiaoFundo);
                AdicionaERemoveParticulas(&malha, 0.05*(malha.minDxDy), 1);
            }
            else if( !strcmp(tipoBoundary, "CABER_AXI") ) {
                double Rp = 0.0, R0 = 0.0, Lz = 0.0;

                testeScan = fscanf(arq, "ALTURA=%lf ", &Lz);
                testeScan = fscanf(arq, "RAIO_CENTRO=%lf ", &R0);
                AdicionaCurvaCaberAxi(&(malha.interface), Rp, R0, Lz, novaRegiao, regiaoFundo);
                AdicionaERemoveParticulas(&malha, 0.05*(malha.minDxDy), 1);
            }
            else if( !strcmp(tipoBoundary, "COLISAO_DUAS_GOTAS") ) {
                double centroX, centroY, raio1, raio2, distanciaGotas;
                double gotas_B;

                testeScan = fscanf(arq, "CENTRO_BASE=[%lf %lf] RAIOS=[%lf %lf] B=%lf DISTANCIA=%lf MERGE_TIME=%lf\n", &centroX, &centroY, &raio1, &raio2, &gotas_B, &distanciaGotas, &(malha.merge_time));
                if( testeScan!=7 ) DeuProblema("\n\nERRO: Problema na leitura dos parametros COLISAO_DUAS_GOTAS (REGIAO_LIVRE).\n\n");

                // Calculando as velocidades de cada gota
                double b = gotas_B*(raio1 + raio2);
                double angulo = Bisseccao_ColisaoGotas(0.0, 0.5*M_PI, raio1 + raio2, 0.0, b);
                double angulo1 = - angulo;
                double angulo2 = M_PI - angulo;
                double U1[2] = {0.5*cos(angulo1), 0.5*sin(angulo1)};
                double U2[2] = {0.5*cos(angulo2), 0.5*sin(angulo2)};
                angulo = 180.0*angulo/M_PI;
//                DeuProblema("ANGULO = %lf\nVELOCIDADES \n U1 = %lf %lf\n U2 = %lf %lf\n", angulo, U1[0], U1[1], U2[0], U2[1]);

                double centroXgota1 = (centroX - raio1) - 2.0*U1[0]*distanciaGotas;
                double centroYgota1 = (centroY) - 2.0*U1[1]*distanciaGotas;
                double centroXgota2 = (centroX + raio2) - 2.0*U2[0]*distanciaGotas;
                double centroYgota2 = (centroY) - 2.0*U2[1]*distanciaGotas;

                // Uma das gotas
                AdicionaCurvaElipse(&(malha.interface), centroXgota1, centroYgota1, raio1, raio1, novaRegiao, regiaoFundo, M_PI, 200);
                InsereBlocoFreeSurface(&malha, TIPO_FREE_ELIPSE, U1[0], U1[1], 0.0, 0.0, 0.0, 0.0, raio1, raio1, centroXgota1, centroYgota1);

                // Segunda gota
                novaRegiao = (REGIAO *)malloc( sizeof(REGIAO) );
                novaRegiao->localRegiao = INTERNA;
                novaRegiao->numero = qtdRegioes++;
                novaRegiao->status = FISICA;
                AdicionaCurvaElipse(&(malha.interface), centroXgota2, centroYgota2, raio2, raio2, novaRegiao, regiaoFundo, 0.0, 200);
                InsereBlocoFreeSurface(&malha, TIPO_FREE_ELIPSE, U2[0], U2[1], 0.0, 0.0, 0.0, 0.0, raio2, raio2, centroXgota2, centroYgota2);

                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);
            }
            else if( !strcmp(tipoBoundary, "SESSILE_PENDING") ) {
                double centroX, centroY, raioBaixo, raioCima, minY, maxY;
                int baseFixa = 0;

                testeScan = fscanf(arq, "CENTRO_BASE=[%lf %lf] RAIOS=[%lf %lf] MINY=%lf MAXY=%lf\n", &centroX, &centroY, &raioBaixo, &raioCima, &minY, &maxY);
                if( testeScan!=6 ) DeuProblema("\n\nERRO: Problema na leitura dos parametros SESSILE_PENDING (REGIAO_LIVRE).\n\n");

                AdicionaCurvaSessilePendingDrops(&(malha.interface), raioCima, raioBaixo, centroX, centroY, minY, maxY, baseFixa);
                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);

            }
            else if( !strcmp(tipoBoundary, "FARADAY_WAVES") ) {
                double alfa, epsilon;

                testeScan = fscanf(arq, "ALFA=%lf EPSILON=%lf\n", &alfa, &epsilon);
                if( testeScan!=2 ) DeuProblema("\n\nERRO: Problema na leitura dos parametros FARADAY_WAVES (REGIAO_LIVRE).\n\n");

                AdicionaCurvaFaradayWaves(&(malha.interface), alfa, epsilon, novaRegiao, regiaoFundo);
                AdicionaERemoveParticulas(&malha, 0.2*malha.minDxDy, 1);
            }
            else if ( !strcmp(tipoBoundary, "EXTRAIR_CONTORNO_IMAGEM") ) {           
                char filepath[101], arqPontos[101];
                float distPontos, xMin, xMax, yMin, yMax;
                double freeVelX, freeVelY;
                PyObject *pModulo, *pFunc;
                PyObject *pArg;

		        // Define o nome do arquivo de pontos
		        arqPontos[0] = '\0';
                strcat(arqPontos, malha.nomeArquivo);
                strcat(arqPontos, ".txt");
                arqPontos[strlen(arqPontos)] = '\0';

                testeScan = fscanf(arq, "FILEPATH=%s DIST_PONTOS=%f EXTREMOS=[%f %f %f %f] VEL=[%lf %lf]\n", filepath, &distPontos, &xMin, &xMax, &yMin, &yMax, &freeVelX, &freeVelY);

                if (testeScan!=8) DeuProblema("\n\nERRO: Problema na leitura dos parametros EXTRAIR_CONTORNO_IMAGEM (REGIAO_LIVRE).\n\n");

                Py_Initialize();

                PySys_SetPath(L"/home/diego/codeflow/CodeFlowAxissimetrico/ExtracaoContorno/");
                pModulo = PyImport_Import(PyUnicode_FromString((char *) "script"));
                if (pModulo != NULL) {
                    pFunc = PyObject_GetAttrString(pModulo, (char *) "ExtracaoContorno");
                    if (pFunc && PyCallable_Check(pFunc)) {
                        pArg = PyTuple_New(7);
                        PyTuple_SetItem(pArg, 0, PyUnicode_FromString((char *) filepath));
                        PyTuple_SetItem(pArg, 1, PyUnicode_FromString((char *) malha.nomeArquivo));
                        PyTuple_SetItem(pArg, 2, PyFloat_FromDouble(distPontos));
                        PyTuple_SetItem(pArg, 3, PyFloat_FromDouble(xMin));
                        PyTuple_SetItem(pArg, 4, PyFloat_FromDouble(xMax));
                        PyTuple_SetItem(pArg, 5, PyFloat_FromDouble(yMin));
                        PyTuple_SetItem(pArg, 6, PyFloat_FromDouble(yMax));
                        PyObject_CallObject(pFunc, pArg);
                    } else {
                        PyErr_Print();
                        fprintf(stderr, "\nNão foi possível encontrar a função da extração de pontos da imagem\n");
                    }
                } else {
                    PyErr_Print();
                    fprintf(stderr, "\nFalha na extração dos pontos da imagem\n");
                }
	
                Py_Finalize();
		
		
		
                AdicionaCurvaImagemAxissimetrico(&(malha.interface), arqPontos, novaRegiao, regiaoFundo);
                AdicionaERemoveParticulas(&malha, 0.05*malha.minDxDy, 1);
                InsereBlocoFreeSurface(&malha, TIPO_FREE_IMAGEM, freeVelX, freeVelY, xMin, xMax, yMin, yMax, 0.0, 0.0, 0.0, 0.0);

		        // Deleta o arquivo de pontos
		        remove(arqPontos);  
            }
        }
        else if( !strcmp(stringLida, "PARAMETRO") ) {
            char stringParametro[500];
            testeScan = fscanf(arq, "%s\n", stringParametro);
            InsereParametroOpcional(&(malha.parametrosOpcionais), stringParametro);
        }
        else if( !strcmp(stringLida, "BLOCO_PARAVIEW") ) {

            testeScan = fscanf(arq, "NUMFACES=%d FACES=", &numFaces);
            if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do NUMFACES de um BLOCO_PARAVIEW.\n\n");

            faces = (int *)malloc( numFaces*sizeof(int) );
            for( i=0; i<numFaces; i++ ) {
                testeScan = fscanf(arq, "%d ", &faces[i]);
                if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura nos IDs das FACES de um BLOCO_PARAVIEW.\n\n");

                PrintDebug("%d\n", faces[i]);
            }

            AdicionaBlocoParaview(&malha, &blocoParaview, numFaces, faces);
	    ArmazenaQuinasBlocoParaview(&malha, blocoParaview, numFaces);//getchar();
            free(faces);

        }

    }

    //Se tiver superficie livre, cria as celulas e pontos usados como flags pra equacao do contorno
    char valorInicial = 0;
    if( malha.interface.curvaInicial!=NULL ) {
        malha.pontosContinuidade = (char **)AlocaMatriz(malha.Nx, malha.Ny, sizeof(char), &valorInicial);
        malha.pontosSupLivre = (char **)AlocaMatriz(malha.Nx+1, malha.Ny+1, sizeof(char), &valorInicial);
    }

    InicializaIndiceIncognitas(&malha);

    //Fazendo uma ultima passagem por todas as celulas pra determinar quais contornos sao cima baixo esq dir
    InicializaExtremosContorno(&malha);


    //Soh liberando aquela regiao la de cima caso nao tenha superficie livre
    if( malha.interface.curvaInicial==NULL )
        free(regiaoFundo);

    fclose(arqOriginal);
    fclose(arq);
    char stringDeletar[950];
    sprintf(stringDeletar, "rm %s", nomeArqJunk);
    testeScan = system(stringDeletar);


    ///Imprimindo os dados do modelo na tela    
    PrintDebug("\n DADOS DA SIMULACAO: \n");
    PrintDebug("Nx = %d   Ny = %d\n", malha.Nx, malha.Ny);
    PrintDebug("Extremos: xMin = %g   xMax = %g   yMin = %g   yMax = %g\n", malha.xMin, malha.xMax, malha.yMin, malha.yMax);
    PrintDebug("Espacamento Minimo: %g\n", malha.minDxDy);
    PrintDebug("Passo Temporal: dt = %e\n", malha.dt);
    PrintDebug("Parametros Adimensionais: Re = %g   Fr = %g   We = %g   Beta = %g   Alpha = %g   Epsilon = %g\n", malha.Re, malha.Fr, malha.We, malha.beta, malha.alpha, malha.epsilon);
    PrintDebug("Tolerancia NSF: %e    Tolerancia NSF-Conversao: %e\n", malha.tolNSF, malha.tolNSFconversao);
    PrintDebug("Grau Hibrido: %d\n", malha.grauHibrido);
    if( malha.formulacaoVisco==VISCO_CART )
        PrintDebug("Formulacao Viscoelastica: cartesiano\n");
    if( malha.formulacaoVisco==VISCO_NSF )
        PrintDebug("Formulacao Viscoelastica: NSF\n");
    if( malha.formulacaoVisco==VISCO_HIBRIDO )
        PrintDebug("Formulacao Viscoelastica: hibrido\n");
    if( malha.formulacaoVisco==VISCO_LOG )
        PrintDebug("Formulacao Viscoelastica: log\n");
    if( malha.tensaoSuperficial )
        PrintDebug("Tensao Superficial: ligado com Weber=%lf\n", malha.weber);
    else
        PrintDebug("Tensao Superficial: desligado\n");




    ///Criando a pasta onde ficaram salvos os arquivos VTK desta simulacao
    char nomeLocal[500], temp[100], comando[600];
    strcpy(nomeLocal, malha.nomeArquivo);
    sprintf(temp, "%g", malha.dt);
    while( SubstituiString(nomeLocal, "[dt]", temp) );
    sprintf(temp, "%g", malha.beta);
    while( SubstituiString(nomeLocal, "[beta]", temp) );
    sprintf(temp, "%g", malha.alpha);
    while( SubstituiString(nomeLocal, "[alpha]", temp) );
    sprintf(temp, "%g", malha.epsilon);
    while( SubstituiString(nomeLocal, "[epsilon]", temp) );
    sprintf(temp, "%g", malha.We);
    while( SubstituiString(nomeLocal, "[We]", temp) );
    sprintf(temp, "%1.0e", malha.tolNSF);
    while( SubstituiString(nomeLocal, "[tolNSF]", temp) );
    sprintf(temp, "%1.0e", malha.tolNSFconversao);
    while( SubstituiString(nomeLocal, "[tolNSFconv]", temp) );
    sprintf(temp, "%g", malha.Re);
    while( SubstituiString(nomeLocal, "[Re]", temp) );
    sprintf(temp, "%g", malha.weber);
    while( SubstituiString(nomeLocal, "[Weber]", temp) );
    sprintf(temp, "%g", malha.Fr);
    while( SubstituiString(nomeLocal, "[Froude]", temp) );
    sprintf(temp, "%g", malha.Bi);
    while( SubstituiString(nomeLocal, "[Bingham]", temp) );
    sprintf(temp, "%d", malha.grauHibrido);
    while( SubstituiString(nomeLocal, "[grauHIB]", temp) );
    sprintf(temp, "%g", malha.minDxDy);
    while( SubstituiString(nomeLocal, "[minDxDy]", temp) );
    sprintf(temp, "%d", malha.freqTSUR);
    while( SubstituiString(nomeLocal, "[tsur]", temp) );
    sprintf(temp, "%g", malha.gravX);
    while( SubstituiString(nomeLocal, "[gravX]", temp) );
    sprintf(temp, "%g", malha.gravY);
    while( SubstituiString(nomeLocal, "[gravY]", temp) );

    if( malha.formulacaoVisco==VISCO_CART )
        sprintf(temp, "%s", "CSF");
    else if( malha.formulacaoVisco==VISCO_NSF )
        sprintf(temp, "%s", "NSF");
    else if( malha.formulacaoVisco==VISCO_HIBRIDO )
        sprintf(temp, "%s", "HIB");
    else if( malha.formulacaoVisco==VISCO_LOG )
        sprintf(temp, "%s", "LOG");
    else DeuProblema("LerModelo: formulacao visco invalida.\n\n");
    while( SubstituiString(nomeLocal, "[formVisco]", temp) );

    PrintDebug("Pasta de Arquivos: %s\n", nomeLocal);
    PrintDebug("\n\n");

    sprintf(comando, "mkdir ArquivosPlot/VTK/%s", nomeLocal);
    testeScan = system(comando);
    puts("\n\n");
    strcpy(malha.nomeArquivo, nomeLocal);

    EscreveArquivoBlocosParaview(blocoParaview, malha);


    // Copiando este arquivo .sim para a pasta, caso eu precise debugar depois
    char comandoCopia[500];
    int indiceBarra = -1;
    for( i=0; NomeArquivo[i]!='\0'; i++ ) {
        if( NomeArquivo[i]=='/' )
            indiceBarra = i;
    }

    sprintf(comandoCopia, "cp %s ArquivosPlot/VTK/%s/%s", NomeArquivo, malha.nomeArquivo, NomeArquivo + indiceBarra + 1);
    testeScan = system(comandoCopia);

    return malha;
}

void InsereBloco(MALHA *M, int id, char *TipoCelula, char *DirecaoNormal, char *TipoBoundary, double ValorDirichlet, double ValorDirichlet2, char *PerfilInflow,
                 char *Movimento, int iMin, int iMax, int jMin, int jMax,
                 double X1, double Y1, double X2, double Y2)
{
    TIPO_CELULA tipoBloco;
    TIPO_BOUNDARY tipoBoundary = NAO_BOUNDARY;
    PERFIL_INFLOW perfilInflow = NAO_INFLOW;
    int i, j;
    DIRECAO novaDirecao;
    MOVIMENTO_PAREDE movimento = -1;
    BLOCO *novoBloco;

    tipoBloco = ReconheceTipoDaStringCelula(TipoCelula);

    if( tipoBloco == BOUNDARY )
        tipoBoundary = ReconheceTipoDaStringBoundary(TipoBoundary);
    if( tipoBoundary == INFLOW )
        perfilInflow = ReconheceTipoDaStringPerfilInflow(PerfilInflow);

    novaDirecao = ReconheceTipoDaStringDirecao(DirecaoNormal);

    //Usado no caso NO_SLIP e SLIP, caso tenha velocidade na parede
    if( tipoBoundary==NOSLIP || tipoBoundary==SLIP )
        movimento = ReconheceTipoDaStringMovimento(Movimento);

    for(i=iMin; i<=iMax; i++) {
        for(j=jMin; j<=jMax; j++) {
//            direcaoAtual = M->celulas[i][j].direcaoNormal;

            if( novaDirecao == DIRECAO_X ) {
                M->pontosU[i][j].tipo = tipoBloco;
                M->pontosU[i][j].tipoBoundary = tipoBoundary;

                if( tipoBoundary==INFLOW )
                    M->pontosU[i][j].valorDirichlet = ValorDirichlet;
                else if( tipoBoundary==KNOWN_PRESSURE ) {
                    M->pontosU[i][j].valorDirichlet = ValorDirichlet;
                    M->pontosU[i][j].valorDirichlet2 = ValorDirichlet2;
                }
                else if( tipoBoundary==NOSLIP ) {
                    M->pontosU[i][j].valorDirichlet = (movimento!=SEM_MOVIMENTO) ? ValorDirichlet : 0.0;
                    M->pontosU[i][j].movimento = movimento;

                    // double iniciaDouble = 0.0;
                    // M->pontosU[i][j].tensorParedeVelho = (double **)AlocaMatriz(4, 4, sizeof(double), &iniciaDouble);
                    // M->pontosU[i][j].tensorParedeNovo = (double **)AlocaMatriz(4, 4, sizeof(double), &iniciaDouble);
                }
                else if( tipoBoundary==SLIP ) {
                    M->pontosU[i][j].valorDirichlet = ValorDirichlet;
                    M->pontosU[i][j].movimento = movimento;
                }
            }
            else {
                M->pontosV[i][j].tipo = tipoBloco;
                M->pontosV[i][j].tipoBoundary = tipoBoundary;

                if( tipoBoundary==INFLOW )
                    M->pontosV[i][j].valorDirichlet = ValorDirichlet;
                else if( tipoBoundary==KNOWN_PRESSURE ) {
                    M->pontosV[i][j].valorDirichlet = ValorDirichlet;
                    M->pontosV[i][j].valorDirichlet2 = ValorDirichlet2;
                }
                else if( tipoBoundary==NOSLIP ) {
                    M->pontosV[i][j].valorDirichlet = (movimento!=SEM_MOVIMENTO) ? ValorDirichlet : 0.0;
                    M->pontosV[i][j].movimento = movimento;

                    // double iniciaDouble = 0.0;
                    // M->pontosV[i][j].tensorParedeVelho = (double **)AlocaMatriz(4, 4, sizeof(double), &iniciaDouble);
                    // M->pontosV[i][j].tensorParedeNovo = (double **)AlocaMatriz(4, 4, sizeof(double), &iniciaDouble);
                }
                else if( tipoBoundary==SLIP ) {
                    M->pontosV[i][j].valorDirichlet = ValorDirichlet;
                    M->pontosV[i][j].movimento = movimento;
                }
            }

        }
    }

    //Inserindo o bloco na lista de blocos
    novoBloco = (BLOCO *)malloc( sizeof(BLOCO) );
    novoBloco->id = id;
    novoBloco->direcaoNormal = novaDirecao;
    novoBloco->tipoU = novoBloco->tipoV = EMPTY;
    novoBloco->tipoBoundaryU = novoBloco->tipoBoundaryV = NAO_BOUNDARY;
    if( novaDirecao==DIRECAO_X ) {
        novoBloco->tipoU = tipoBloco;
        novoBloco->tipoBoundaryU = tipoBoundary;
    }
    else {
        novoBloco->tipoV = tipoBloco;
        novoBloco->tipoBoundaryV = tipoBoundary;
    }
    novoBloco->x1 = X1;
    novoBloco->y1 = Y1;
    novoBloco->x2 = X2;
    novoBloco->y2 = Y2;

    novoBloco->movimentoParede = movimento;
    novoBloco->valorDirichlet = ValorDirichlet;
    novoBloco->valorDirichlet2 = ValorDirichlet2;
    novoBloco->iMax = iMax;
    novoBloco->iMin = iMin;
    novoBloco->jMax = jMax;
    novoBloco->jMin = jMin;
    novoBloco->perfilInflow = perfilInflow;
    novoBloco->prox = M->primBloco;
    M->primBloco = novoBloco;

    return;
}

void CriaMalha(MALHA *Malha, const char *TipoStretchingX, int QtdVertX, double *VerticesX, int *Nx, double *ParametrosX, const char *TipoStretchingY, int QtdVertY, double *VerticesY, int *Ny, double *ParametrosY)
{
    double h0, q, Lx, Ly, *pol, espacamento;
    int i, j, indice, celulaBloco, bloco, sentido;

    Malha->minDxDy = 1e+10;

    /// ==== Criando a malha na direcao X
    if( !strcmp(TipoStretchingX, "UNIFORME") ) {
        Malha->Nx = Nx[0];

        Malha->dx = (double *)malloc( Malha->Nx*sizeof(double) );
        Malha->r = (double *)malloc( Malha->Nx*sizeof(double) );
        Malha->x = (double *)malloc( (Malha->Nx+1)*sizeof(double) );

        i = 0;
        Malha->x[0] = Malha->xMin;
        espacamento = (Malha->xMax - Malha->xMin)/(Malha->Nx);
        for( i++; i<=Malha->Nx; i++ ) {
            Malha->dx[i-1] = espacamento;
            Malha->x[i] = Malha->x[i-1] + espacamento;

            if( Malha->dx[i-1] < Malha->minDxDy )
                Malha->minDxDy = Malha->dx[i-1];

            if( Malha->tipoCoord==AXI_CILINDRICO )
                Malha->r[i-1] = 0.5*(Malha->x[i-1] + Malha->x[i]);
            else if( Malha->tipoCoord==CARTESIANO )
                Malha->r[i-1] = 1.0;
            else
                DeuProblema("\nProblema. Sistema de coordenadas inesperado.\n");
        }
        Malha->x[Malha->Nx] = Malha->xMax; //Soh pra eliminar arredondamentos das contas acima
    }
    else if( !strcmp(TipoStretchingX, "GEOMETRICO") ) {
        Malha->Nx = 0;
        for( bloco=0; bloco<QtdVertX-1; bloco++ )
            (Malha->Nx) += Nx[bloco];

        Malha->dx = (double *)malloc( Malha->Nx*sizeof(double) );
        Malha->r = (double *)malloc( Malha->Nx*sizeof(double) );
        Malha->x = (double *)malloc( (Malha->Nx+1)*sizeof(double) );

        Malha->x[0] = VerticesX[0];
        i = 1;
        for( bloco=0; bloco<QtdVertX-1; bloco++ ) {
            Malha->x[i-1+Nx[bloco]] = VerticesX[bloco+1];
            indice = i+Nx[bloco] - 2;

            Lx = (VerticesX[bloco+1] - VerticesX[bloco]);
            h0 = ParametrosX[2*bloco];
            sentido = (ParametrosX[2*bloco+1]>0) ? 1 : 0;
            pol = (double *)malloc(Nx[bloco]*sizeof(double));
            pol[0] = (h0-Lx)/h0;
            for(j=1; j<Nx[bloco]; j++)
                pol[j] = 1.0;
            q = AlgoritmoBisseccao(Nx[bloco]-1, pol, 0, Lx/h0);
            free(pol);

            //Colocando os dx, os x, e os r
            for( celulaBloco=0; celulaBloco<Nx[bloco]; celulaBloco++ ) {

                if( sentido ) {
                    indice = i;


                    Malha->dx[indice-1] = h0;
                    Malha->x[indice] = Malha->x[indice-1] + h0;

                    if( Malha->dx[indice-1] < Malha->minDxDy )
                        Malha->minDxDy = Malha->dx[indice-1];

                    //printf("x[%d] = %lf\n", indice, Malha->x[indice]);

                    if( Malha->tipoCoord==AXI_CILINDRICO )
                        Malha->r[indice-1] = 0.5*(Malha->x[indice] + Malha->x[indice-1]);
                    else if( Malha->tipoCoord==CARTESIANO )
                        Malha->r[indice-1] = 1.0;
                    else
                        DeuProblema("\nProblema. Sistema de coordenadas inesperado.\n");
                }
                else {

                    Malha->dx[indice] = h0;
                    Malha->x[indice] = Malha->x[indice+1] - h0;

                    if( Malha->dx[indice] < Malha->minDxDy )
                        Malha->minDxDy = Malha->dx[indice];

                    //printf("x[%d] = %lf\n", indice, Malha->x[indice]);

                    if( Malha->tipoCoord==AXI_CILINDRICO )
                        Malha->r[indice] = 0.5*(Malha->x[indice] + Malha->x[indice+1]);
                    else if( Malha->tipoCoord==CARTESIANO )
                        Malha->r[indice] = 1.0;
                    else
                        DeuProblema("\nProblema. Sistema de coordenadas inesperado.\n");
                }

                h0 = h0*q;
                i++;
                indice--;
            }

            PrintDebug("Razao X: %lf --- dx_final: %lf\n", q, h0/q);
        }

    }
    else
        DeuProblema("\n\n ERRO: O Tipo de malha fornecido eh invalido. Opcoes validas: UNIFORME e GEOMETRICO. \n\n");





    /// ==== Criando a malha na direcao Y
    if( !strcmp(TipoStretchingY, "UNIFORME") ) {
        Malha->Ny = Ny[0];

        Malha->dy = (double *)malloc( Malha->Ny*sizeof(double) );
        Malha->y = (double *)malloc( (Malha->Ny+1)*sizeof(double) );


        j = 0;
        Malha->y[0] = Malha->yMin;
        espacamento = (Malha->yMax - Malha->yMin)/(Malha->Ny);
        for( j++; j<=Malha->Ny; j++ ) {
            Malha->dy[j-1] = espacamento;
            Malha->y[j] = Malha->y[j-1] + espacamento;

            if( Malha->dy[j-1] < Malha->minDxDy )
                Malha->minDxDy = Malha->dy[j-1];
        }
        Malha->y[Malha->Ny] = Malha->yMax; //Soh pra eliminar arredondamentos das contas acima
    }
    else if( !strcmp(TipoStretchingY, "GEOMETRICO") ) {
        Malha->Ny = 0;
        for( bloco=0; bloco<QtdVertY-1; bloco++ )
            (Malha->Ny) += Ny[bloco];

        Malha->dy = (double *)malloc( Malha->Ny*sizeof(double) );
        Malha->y = (double *)malloc( (Malha->Ny+1)*sizeof(double) );




        Malha->y[0] = VerticesY[0];
        j = 1;
        for( bloco=0; bloco<QtdVertY-1; bloco++ ) {
            Malha->y[j-1+Ny[bloco]] = VerticesY[bloco+1];
            indice = j+Ny[bloco] - 2;

            Ly = (VerticesY[bloco+1] - VerticesY[bloco]);
            h0 = ParametrosY[2*bloco];
            sentido = (ParametrosY[2*bloco+1]>0) ? 1 : 0;
            pol = (double *)malloc(Ny[bloco]*sizeof(double));
            pol[0] = (h0-Ly)/h0;
            for(i=1; i<Ny[bloco]; i++)
                pol[i] = 1.0;
            q = AlgoritmoBisseccao(Ny[bloco]-1, pol, 0, Ly/h0);
            free(pol);

            for( celulaBloco=0; celulaBloco<Ny[bloco]; celulaBloco++ ) {

                if( sentido ) {
                    Malha->dy[j-1] = h0;
                    Malha->y[j] = Malha->y[j-1] + h0;

                    if( Malha->dy[j-1] < Malha->minDxDy )
                        Malha->minDxDy = Malha->dy[j-1];

                    h0 = q*h0;
                    j++;
                }
                else {
                    Malha->dy[indice] = h0;
                    Malha->y[indice] = Malha->y[indice+1] - h0;

                    if( Malha->dy[indice] < Malha->minDxDy )
                        Malha->minDxDy = Malha->dy[indice];

                    h0 = q*h0;
                    j++;
                    indice--;
                }

            }

            PrintDebug("Razao Y: %lf --- dx_final: %lf\n", q, h0/q);
        }
    }
    else
        DeuProblema("\n\n ERRO: O Tipo de malha fornecido eh invalido. Opcoes validas: UNIFORME e GEOMETRICO. \n\n");



    INFO_PONTO info;
    int valorInicial = 0;
    info.tipo = EMPTY;
    info.tipoBoundary = NAO_BOUNDARY;
    info.movimento = SEM_MOVIMENTO;
    info.extremo = NAO_EXTREMO;
    info.valorDirichlet = 0.0;
    // info.tensorParedeVelho = NULL;
    // info.tensorParedeNovo = NULL;

    Malha->celulas = (TIPO_CELULA **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(TIPO_CELULA), &valorInicial);
    Malha->celulasParede = (TIPO_CELULA **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(TIPO_CELULA), &valorInicial);
    Malha->pontosU = (INFO_PONTO **)AlocaMatriz(Malha->Nx+1, Malha->Ny, sizeof(INFO_PONTO), &info);
    Malha->pontosV = (INFO_PONTO **)AlocaMatriz(Malha->Nx, Malha->Ny+1, sizeof(INFO_PONTO), &info);
    if( Malha->tipoCoord == AXI_CILINDRICO )
        Malha->pontosW = (INFO_PONTO **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(INFO_PONTO), &info);
    Malha->pontosP = (INFO_PONTO **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(INFO_PONTO), &info);

    double valorDouble = 0.0;
    //Armazena pontos bloco paraview
    //Malha->PontosBlocoParaview = (double **)AlocaMatriz(4, 2, sizeof(double), &valorDouble);
    //--
    PONTO *valorNULL = NULL;
    Malha->curvaturaCelulas = (double **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(double), &valorDouble);
    Malha->pontoCelula = (PONTO ***)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(PONTO *), &valorNULL);
    Malha->vetNx = (double **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(double), &valorDouble);
    Malha->vetNy = (double **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(double), &valorDouble);


    valorInicial = -1;
    Malha->matIndicesU = (int **)AlocaMatriz(Malha->Nx+1, Malha->Ny, sizeof(int), &valorInicial);
    Malha->matIndicesV = (int **)AlocaMatriz(Malha->Nx, Malha->Ny+1, sizeof(int), &valorInicial);
    Malha->matIndicesP = (int **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(int), &valorInicial);
    if( Malha->tipoCoord==AXI_CILINDRICO )
        Malha->matIndicesW = (int **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(int), &valorInicial);

    Malha->primBloco = NULL;
    Malha->blocoLivre = NULL;

    Malha->indicesU = NULL;
    Malha->indicesV = NULL;
    Malha->indicesW = NULL;
    Malha->indicesP = NULL;

    Malha->parametrosOpcionais.prim = NULL;

    Malha->malhaRecursoes = (int **)AlocaMatriz(Malha->Nx, Malha->Ny, sizeof(int), &valorInicial);

    return;
}

void CriaMalhaUniforme(MALHA *Malha, int Nx, int Ny)
{
    double espacamento;
    int i, j;

    Malha->minDxDy = 1e+10;

    /// ==== Primeiro vamos discretizar o eixo X
    Malha->Nx = Nx;

    Malha->dx = (double *)malloc( Malha->Nx*sizeof(double) );
    Malha->r = (double *)malloc( Malha->Nx*sizeof(double) );
    Malha->x = (double *)malloc( (Malha->Nx+1)*sizeof(double) );

    i = 0;
    Malha->x[0] = Malha->xMin;
    espacamento = (Malha->xMax - Malha->xMin)/(Malha->Nx);
    for( i++; i<=Malha->Nx; i++ ) {
        Malha->dx[i-1] = espacamento;
        Malha->x[i] = Malha->x[i-1] + espacamento;

        if( Malha->dx[i-1] < Malha->minDxDy )
            Malha->minDxDy = Malha->dx[i-1];

        if( Malha->tipoCoord==AXI_CILINDRICO )
            Malha->r[i-1] = 0.5*(Malha->x[i-1] + Malha->x[i]);
        else if( Malha->tipoCoord==CARTESIANO )
            Malha->r[i-1] = 1.0;
        else {
            DeuProblema("\nProblema. Caso inesperado.\n");
        }
    }



    /// ==== Agora  vamos discretizar o eixo Y
    Malha->Ny = Ny;

    Malha->dy = (double *)malloc( Malha->Ny*sizeof(double) );
    Malha->y = (double *)malloc( (Malha->Ny+1)*sizeof(double) );


    j = 0;
    Malha->y[0] = Malha->yMin;
    espacamento = (Malha->yMax - Malha->yMin)/(Malha->Ny);
    for( j++; j<=Malha->Ny; j++ ) {
        Malha->dy[j-1] = espacamento;
        Malha->y[j] = Malha->y[j-1] + espacamento;

        if( Malha->dy[j-1] < Malha->minDxDy )
            Malha->minDxDy = Malha->dy[j-1];
    }

    return;
}

void CriaMalhaStretchingTrygvasson(MALHA *Malha, int QtdVertX, double *VerticesX, int *Nx, double *ParametrosX, int QtdVertY, double *VerticesY, int *Ny, double *ParametrosY)
{
    double Lx, Ly, Ax, Ay, Xc, Yc, ds, s;
    double pontoAnt, pontoNovo;
    int i, j, bloco;

    Malha->minDxDy = 1e+8;

    /// ==== Primeiro vamos discretizar o eixo X
    Malha->Nx = 0;
    for( bloco=0; bloco<QtdVertX-1; bloco++ )
        (Malha->Nx) += Nx[bloco];

    Malha->dx = (double *)malloc( Malha->Nx*sizeof(double) );
    Malha->r = (double *)malloc( Malha->Nx*sizeof(double) );
    Malha->x = (double *)malloc( (Malha->Nx+1)*sizeof(double) );

    i = 0;
    Malha->x[0] = VerticesX[0];
    for( bloco=0; bloco<QtdVertX-1; bloco++ ) {
        Lx = (VerticesX[bloco+1] - VerticesX[bloco]);
        Ax = ParametrosX[2*bloco];
        Xc = ParametrosX[2*bloco + 1];

        //Colocando os dx, os x, e os r
        pontoAnt = VerticesX[bloco];
        ds = 1.0/Nx[bloco];
        s = ds;
        for( i++; i<=Nx[bloco]; i++ ) {
            pontoNovo = VerticesX[bloco] + Lx*s + Ax*(Xc - Lx*s)*s*(1-s);
            Malha->dx[i-1] = pontoNovo - pontoAnt;
            Malha->x[i] = pontoNovo;

            if( Malha->dx[i-1] < Malha->minDxDy )
                Malha->minDxDy = Malha->dx[i-1];

            if( Malha->tipoCoord==AXI_CILINDRICO )
                Malha->r[i-1] = 0.5*(pontoAnt + pontoNovo);
            else if( Malha->tipoCoord==CARTESIANO )
                Malha->r[i-1] = 1.0;
            else {
                DeuProblema("\nProblema. Caso inesperado.\n");
            }

            pontoAnt = pontoNovo;
            s += ds;
        }
    }



    /// ==== Agora  vamos discretizar o eixo Y
    Malha->Ny = 0;
    for( bloco=0; bloco<QtdVertY-1; bloco++ )
        (Malha->Ny) += Ny[bloco];

    Malha->dy = (double *)malloc( Malha->Ny*sizeof(double) );
    Malha->y = (double *)malloc( (Malha->Ny+1)*sizeof(double) );



    j = 0;
    Malha->y[0] = VerticesY[0];
    for( bloco=0; bloco<QtdVertY-1; bloco++ ) {
        Ly = (VerticesY[bloco+1] - VerticesY[bloco]);
        Ay = ParametrosY[2*bloco];
        Yc = ParametrosY[2*bloco + 1];


        pontoAnt = VerticesY[bloco];
        ds = 1.0/Ny[bloco];
        s = ds;
        for( j=1; j<=Ny[bloco]; j++ ) {
            pontoNovo = VerticesY[bloco] + Ly*s + Ay*(Yc - Ly*s)*s*(1-s);
            Malha->dy[j-1] = pontoNovo - pontoAnt;
            Malha->y[j] = pontoNovo;

            if( Malha->dy[j-1] < Malha->minDxDy )
                Malha->minDxDy = Malha->dy[j-1];

            pontoAnt = pontoNovo;
            s += ds;
        }
    }

    return;
}

void CriaMalhaStretchingGeometrico(MALHA *Malha, int QtdVertX, double *VerticesX, int *Nx, double *ParametrosX, int QtdVertY, double *VerticesY, int *Ny, double *ParametrosY)
{
    double h0, q, Lx, Ly, *pol;
    int i, j, indice, celulaBloco, bloco, sentido;

    Malha->minDxDy = 1e+8;

    /// ==== Primeiro vamos discretizar o eixo X
    Malha->Nx = 0;
    for( bloco=0; bloco<QtdVertX-1; bloco++ )
        (Malha->Nx) += Nx[bloco];

    Malha->dx = (double *)malloc( Malha->Nx*sizeof(double) );
    Malha->r = (double *)malloc( Malha->Nx*sizeof(double) );
    Malha->x = (double *)malloc( (Malha->Nx+1)*sizeof(double) );

    Malha->x[0] = VerticesX[0];
    i = 1;
    for( bloco=0; bloco<QtdVertX-1; bloco++ ) {
        Malha->x[i-1+Nx[bloco]] = VerticesX[bloco+1];
        indice = i+Nx[bloco] - 2;

        Lx = (VerticesX[bloco+1] - VerticesX[bloco]);
        h0 = ParametrosX[2*bloco];
        sentido = (ParametrosX[2*bloco+1]>0) ? 1 : 0;
        pol = (double *)malloc(Nx[bloco]*sizeof(double));
        pol[0] = (h0-Lx)/h0;
        for(j=1; j<Nx[bloco]; j++)
            pol[j] = 1.0;
        q = AlgoritmoBisseccao(Nx[bloco]-1, pol, 0, Lx/h0);
        free(pol);

        PrintDebug("Razao X: %lf\n", q);

        //Colocando os dx, os x, e os r
        for( celulaBloco=0; celulaBloco<Nx[bloco]; celulaBloco++ ) {

            if( sentido ) {
                indice = i;


                Malha->dx[indice-1] = h0;
                Malha->x[indice] = Malha->x[indice-1] + h0;

                if( Malha->dx[indice-1] < Malha->minDxDy )
                    Malha->minDxDy = Malha->dx[indice-1];

                //printf("x[%d] = %lf\n", indice, Malha->x[indice]);

                if( Malha->tipoCoord==AXI_CILINDRICO )
                    Malha->r[indice-1] = 0.5*(Malha->x[indice] + Malha->x[indice-1]);
                else if( Malha->tipoCoord==CARTESIANO )
                    Malha->r[indice-1] = 1.0;
                else {
                    DeuProblema("\nProblema. Caso inesperado.\n");
                }
            }
            else {

                Malha->dx[indice] = h0;
                Malha->x[indice] = Malha->x[indice+1] - h0;

                if( Malha->dx[indice] < Malha->minDxDy )
                    Malha->minDxDy = Malha->dx[indice];

                //printf("x[%d] = %lf\n", indice, Malha->x[indice]);

                if( Malha->tipoCoord==AXI_CILINDRICO )
                    Malha->r[indice] = 0.5*(Malha->x[indice] + Malha->x[indice+1]);
                else if( Malha->tipoCoord==CARTESIANO )
                    Malha->r[indice] = 1.0;
                else {
                    DeuProblema("\nProblema. Caso inesperado.\n");
                }
            }

            h0 = h0*q;
            i++;
            indice--;
        }
    }

    /// ==== Agora  vamos discretizar o eixo Y
    Malha->Ny = 0;
    for( bloco=0; bloco<QtdVertY-1; bloco++ )
        (Malha->Ny) += Ny[bloco];

    Malha->dy = (double *)malloc( Malha->Ny*sizeof(double) );
    Malha->y = (double *)malloc( (Malha->Ny+1)*sizeof(double) );




    Malha->y[0] = VerticesY[0];
    j = 1;
    for( bloco=0; bloco<QtdVertY-1; bloco++ ) {
        Malha->y[j-1+Ny[bloco]] = VerticesY[bloco+1];
        indice = j+Ny[bloco] - 2;

        Ly = (VerticesY[bloco+1] - VerticesY[bloco]);
        h0 = ParametrosY[2*bloco];
        sentido = (ParametrosY[2*bloco+1]>0) ? 1 : 0;
        pol = (double *)malloc(Ny[bloco]*sizeof(double));
        pol[0] = (h0-Ly)/h0;
        for(i=1; i<Ny[bloco]; i++)
            pol[i] = 1.0;
        q = AlgoritmoBisseccao(Ny[bloco]-1, pol, 0, Ly/h0);
        free(pol);

        PrintDebug("Razao Y: %lf\n", q);


        for( celulaBloco=0; celulaBloco<Ny[bloco]; celulaBloco++ ) {

            if( sentido ) {
                Malha->dy[j-1] = h0;
                Malha->y[j] = Malha->y[j-1] + h0;

                if( Malha->dy[j-1] < Malha->minDxDy )
                    Malha->minDxDy = Malha->dy[j-1];

                h0 = q*h0;
                j++;
            }
            else {
                Malha->dy[indice] = h0;
                Malha->y[indice] = Malha->y[indice+1] - h0;

                if( Malha->dy[indice] < Malha->minDxDy )
                    Malha->minDxDy = Malha->dy[indice];

                h0 = q*h0;
                j++;
                indice--;
            }

        }
    }

    return;
}



void InicializaMalhaComEmpty(MALHA *M)
{
    int i, j;

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
            if( M->celulas[i][j]!=BOUNDARY )
                M->celulas[i][j] = EMPTY;
            M->matIndicesP[i][j] = -1;
        }
    }

    if( M->tipoCoord==AXI_CILINDRICO ) {
        for( i=0; i<M->Nx; i++ ) {
            for( j=0; j<M->Ny; j++ ) {
                if( M->pontosW[i][j].tipo!=BOUNDARY )
                    M->pontosW[i][j].tipo = EMPTY;
                M->matIndicesW[i][j] = -1;
            }
        }
    }
    return;
}

void InicializaIndiceIncognitas(MALHA *M)
{
    TIPO_CELULA tipo;
    //TIPO_BOUNDARY tipoBoundary;
    int i, j, linha = 0;

    M->qtdIncognitasU = 0;
    M->qtdIncognitasV = 0;
    M->qtdIncognitasW = 0;
    M->qtdIncognitasP = 0;

    //Contando quantas incognitas vai ter no sistema
    for( i=0; i<=M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            tipo = M->pontosU[i][j].tipo;
            //tipoBoundary = M->pontosU[i][j].tipoBoundary;
            if( tipo==FULL || tipo==SURFACE )
                (M->qtdIncognitasU)++;
            // boundary da esquerda
            else if( (tipo==BOUNDARY) && (i!=M->Nx) && (M->pontosU[i+1][j].tipo!=EMPTY) )
                (M->qtdIncognitasU)++;
            // boundary da direita
            else if( (tipo==BOUNDARY) && (i!=0) && (M->pontosU[i-1][j].tipo!=EMPTY) )
                (M->qtdIncognitasU)++;
        }
    }

    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<=M->Ny; j++ ) {
            tipo = M->pontosV[i][j].tipo;
            //tipoBoundary = M->pontosV[i][j].tipoBoundary;
            if( tipo==FULL || tipo==BOUNDARY || tipo==SURFACE ) {
                (M->qtdIncognitasV)++;
            }
        }
    }

    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            tipo = M->pontosP[i][j].tipo;

            //tipoBoundary = M->pontosP[i][j].tipoBoundary;
            if( tipo==FULL || tipo==SURFACE ) {
                (M->qtdIncognitasP)++;
            }
        }
    }

    if( M->tipoCoord==AXI_CILINDRICO ) {
        for( i=0; i<M->Nx; i++ ) {
            for( j=0; j<M->Ny; j++ ) {
                // W classifica igual a P
                M->pontosW[i][j].tipo = M->pontosP[i][j].tipo;

                tipo = M->pontosW[i][j].tipo;
                //tipoBoundary = M->pontosP[i][j].tipoBoundary;
                if( tipo==FULL || tipo==SURFACE ) {
                    (M->qtdIncognitasW)++;
                }
            }
        }
    }

    //Liberando memoria que estiver previamente alocada (se tiver)
    free(M->indicesU);
    free(M->indicesV);
    free(M->indicesP);
    free(M->indicesW);

    //Alocando memoria
    M->indicesU = (POSICAO *)malloc( M->qtdIncognitasU*sizeof(POSICAO) );
    M->indicesV = (POSICAO *)malloc( M->qtdIncognitasV*sizeof(POSICAO) );
    M->indicesP = (POSICAO *)malloc( M->qtdIncognitasP*sizeof(POSICAO) );
    if( M->tipoCoord==AXI_CILINDRICO )
        M->indicesW = (POSICAO *)malloc( M->qtdIncognitasW*sizeof(POSICAO) );

    //Colocando os indices
    for( i=0; i<=M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            tipo = M->pontosU[i][j].tipo;
            //tipoBoundary = M->pontosU[i][j].tipoBoundary;
            if( tipo==FULL || tipo==SURFACE ) {
                M->indicesU[linha].i = i;
                M->indicesU[linha].j = j;
                M->matIndicesU[i][j] = linha;
                linha++;
            }
            // boundary da esquerda
            else if( (tipo==BOUNDARY) && (i!=M->Nx) && (M->pontosU[i+1][j].tipo!=EMPTY) ) {
                M->indicesU[linha].i = i;
                M->indicesU[linha].j = j;
                M->matIndicesU[i][j] = linha;
                linha++;
            }
            // boundary da direita
            else if( (tipo==BOUNDARY) && (i!=0) && (M->pontosU[i-1][j].tipo!=EMPTY) ) {
                M->indicesU[linha].i = i;
                M->indicesU[linha].j = j;
                M->matIndicesU[i][j] = linha;
                linha++;
            }
        }
    }


    //Colocando os indices da velocidade V
    linha = 0;
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<=M->Ny; j++ ) {
            tipo = M->pontosV[i][j].tipo;
            //tipoBoundary = M->pontosV[i][j].tipoBoundary;
            if( tipo==FULL || tipo==BOUNDARY || tipo==SURFACE ) {
                M->indicesV[linha].i = i;
                M->indicesV[linha].j = j;
                M->matIndicesV[i][j] = linha;
                linha++;
            }
        }
    }

    //Colocando os indices da pressao P
    linha = 0;
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            tipo = M->pontosP[i][j].tipo;
            //tipoBoundary = M->pontosP[i][j].tipoBoundary;
            if( tipo==FULL || tipo==SURFACE ) {
                M->indicesP[linha].i = i;
                M->indicesP[linha].j = j;
                M->matIndicesP[i][j] = linha;
                linha++;
            }
        }
    }

    if( M->tipoCoord==AXI_CILINDRICO ) {
        //Colocando os indices da velocidade W
        linha = 0;
        for( i=0; i<M->Nx; i++ ) {
            for( j=0; j<M->Ny; j++ ) {
                tipo = M->pontosW[i][j].tipo;
                //tipoBoundary = M->pontosP[i][j].tipoBoundary;
                if( tipo==FULL || tipo==SURFACE ) {
                    M->indicesW[linha].i = i;
                    M->indicesW[linha].j = j;
                    M->matIndicesW[i][j] = linha;
                    linha++;
                }
            }
        }
    }


    return;
}

void InicializaExtremosContorno(MALHA *M)
{
    int i, j;

    //Colocando os extremos nos contornos tipoU
    for( i=0; i<=M->Nx; i++ ) {
        for( j=0; j<M->Ny; j++ ) {
            M->pontosU[i][j].extremo = NAO_EXTREMO;

            //Contorno vertical. Vai ser esquerda ou direita
            if( M->pontosU[i][j].tipo == BOUNDARY ) {
                if( i==0 ) //Eh esquerda
                    M->pontosU[i][j].extremo = ESQUERDA;
                else if( i==M->Nx ) //Eh direita
                    M->pontosU[i][j].extremo = DIREITA;
                else if( M->pontosU[i-1][j].tipo == EMPTY )
                    M->pontosU[i][j].extremo = ESQUERDA;
                else if( M->pontosU[i+1][j].tipo == EMPTY )
                    M->pontosU[i][j].extremo = DIREITA;
            }
        }
    }

    //Colocando os extremos nos contornos tipoU
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<=M->Ny; j++ ) {
            M->pontosV[i][j].extremo = NAO_EXTREMO;

            //Contorno horizontal. Vai ser esquerda ou direita
            if( M->pontosV[i][j].tipo == BOUNDARY ) {
                if( j==0 )
                    M->pontosV[i][j].extremo = BAIXO;
                else if( j==M->Ny )
                    M->pontosV[i][j].extremo = CIMA;
                else if( M->pontosV[i][j-1].tipo == EMPTY )
                    M->pontosV[i][j].extremo = BAIXO;
                else if( M->pontosV[i][j+1].tipo == EMPTY )
                    M->pontosV[i][j].extremo = CIMA;
            }
        }
    }
}

void LeituraMalhaStretching(FILE *Arq, char *Tipo, int *QtdVertices, double **Vertices, int **Espacamentos, double **Parametros)
{
    char malhaEixo[30];
    char *linha = NULL;
    size_t len = 0;
    ssize_t comprimentoLinha;
    double vertice;
    int indice, indice2, testeScan;


    comprimentoLinha = getline(&linha, &len, Arq);
    fseek(Arq, -comprimentoLinha, SEEK_CUR);
    free(linha);
    linha = NULL;


    testeScan = fscanf(Arq, "%s TIPO=%s ", malhaEixo, Tipo);
    if( testeScan!=2 ) DeuProblema("\n\nERRO: Problema na leitura de algum dos eixos da malha.\n\n");

    if( !strcmp(Tipo, "UNIFORME") ) {
        *Espacamentos = (int *)malloc( 1*sizeof(int) );
        testeScan = fscanf(Arq, "N=%d\n", &((*Espacamentos)[0]));
        if( testeScan!=1 ) DeuProblema("\n\nERRO: Problema na leitura do numero de celulas.\n\n");
        return;
    }


    //Lendo a quantidade de vertices
    testeScan = fscanf(Arq, "VERTICES=(");
    *QtdVertices = 0;
    while(1==fscanf(Arq, "%lf", &vertice))
        (*QtdVertices)++;

    //Alocando memoria
    *Vertices = (double *)malloc( (*QtdVertices)*sizeof(double) );
    *Espacamentos = (int *)malloc( ((*QtdVertices)-1)*sizeof(int) );
    *Parametros = (double *)malloc( 2*((*QtdVertices) - 1)*sizeof(double) );

    testeScan = getline(&linha, &len, Arq);
    fseek(Arq, -comprimentoLinha, SEEK_CUR);
    free(linha);
    linha = NULL;

    //lendo os vertices e parametros
    indice = 0;
    testeScan = fscanf(Arq, "%s TIPO=%s VERTICES=(", malhaEixo, Tipo);
    if( testeScan!=2 ) DeuProblema("\n\nERRO: Problema na leitura de algum dos eixos da malha.\n\n");

    while( 1==fscanf(Arq, "%lf", &((*Vertices)[indice++])) );

    indice = 0;
    indice2 = 0;
    testeScan = fscanf(Arq, ") PARAM_STRETCHING=(");
    while(3==fscanf(Arq, "[%d %lf %lf] ", &((*Espacamentos)[indice2++]), &((*Parametros)[indice]), &((*Parametros)[indice+1])))
        indice = indice + 2;


    testeScan = fscanf(Arq, ")\n");
    return;
}

TIPO_CELULA ReconheceTipoDaStringCelula(char *Cadeia)
{
    if( !strcmp(Cadeia, "FULL") )
        return FULL;
    else if( !strcmp(Cadeia, "EMPTY") )
        return EMPTY;
    else if( !strcmp(Cadeia, "BOUNDARY") )
        return BOUNDARY;
    DeuProblema("\n\n Erro: Tipo de celula nao reconhecido. \n\n");

    return ERRO_CELULA;
}

TIPO_BOUNDARY ReconheceTipoDaStringBoundary(char *Cadeia)
{
    if( !strcmp(Cadeia, "INFLOW") )
        return INFLOW;
    else if( !strcmp(Cadeia, "NEUMANN") )
        return NEUMANN;
    else if( !strcmp(Cadeia, "NOSLIP") )
        return NOSLIP;
    else if( !strcmp(Cadeia, "SLIP") )
        return SLIP;
    else if( !strcmp(Cadeia, "EIXO_AXISSIMETRIA") )
        return EIXO_AXISSIMETRIA;
    else if( !strcmp(Cadeia, "SIMETRIA") )
        return SIMETRIA;
    else if( !strcmp(Cadeia, "PERIODICIDADE") )
        return PERIODICIDADE;
    else if( !strcmp(Cadeia, "KNOWN_PRESSURE") )
        return KNOWN_PRESSURE;
    else DeuProblema("\n\nERRO: O TIPO fornecido pra uma das boundaries nao eh conhecido. Opcoes validas: INFLOW, NEUMANN, NOSLIP.\n\n");

    return ERRO_BOUNDARY;
}

DIRECAO ReconheceTipoDaStringDirecao(char *Cadeia)
{
    if( !strcmp(Cadeia, "x") )
        return DIRECAO_X;
    else if( !strcmp(Cadeia, "y") )
        return DIRECAO_Y;

    return ERRO_DIRECAO;
}

MOVIMENTO_PAREDE ReconheceTipoDaStringMovimento(char *Cadeia)
{
    if( !strcmp(Cadeia, "SEM_MOVIMENTO") )
        return SEM_MOVIMENTO;
    else if( !strcmp(Cadeia, "DIRECAO_X") )
        return MOVIMENTO_X;
    else if( !strcmp(Cadeia, "DIRECAO_Y") )
        return MOVIMENTO_Y;
    else if( !strcmp(Cadeia, "DIRECAO_THETA") )
        return MOVIMENTO_THETA;
    else DeuProblema("\n\nERRO: O MOVIMENTO %s fornecido pra uma das boundaries noslip nao eh conhecido. Opcoes validas: SEM_MOVIMENTO, DIRECAO_X, DIRECAO_Y.\n\n", Cadeia);

    return SEM_MOVIMENTO;
}

PERFIL_INFLOW ReconheceTipoDaStringPerfilInflow(char *Cadeia)
{
    if( !strcmp(Cadeia, "PARABOLICO") )
        return PARABOLICO;
    else if( !strcmp(Cadeia, "PARABOLICO_SIMETRICO_INICIO") )
        return PARABOLICO_SIMETRICO_INICIO;
    else if( !strcmp(Cadeia, "PARABOLICO_SIMETRICO_FIM") )
        return PARABOLICO_SIMETRICO_FIM;
    else if( !strcmp(Cadeia, "RETO") )
        return RETO;

    return ERRO_DIRECAO;
}

void PreencheRegiao(MALHA M, int NumFaces, int *Faces)
{
    int faceEsq = -1;
    double xMin = M.xMax + 1, xMax = M.xMin-1, ponto;
    int i, j, iCelulaInicial, jCelulaInicial;
    BLOCO *bloco;

    for( i=0; i<NumFaces; i++ ) {
        bloco = BuscaBloco( M, Faces[i] );
        if(bloco->direcaoNormal==DIRECAO_X) {
            ponto = M.xMin;
            int i2;
            for( i2=0; i2<bloco->iMin; i2++ )
                ponto += M.dx[i2];

            if( ponto<xMin ) {
                xMin = ponto;
                faceEsq = i;
            }
            if( ponto>xMax ) {
                xMax = ponto;
                //faceDir = i;
            }
        }
    }

    bloco = BuscaBloco(M, Faces[faceEsq]);
    iCelulaInicial = bloco->iMin + 1;
    jCelulaInicial = (bloco->jMin + bloco->jMax)/2;


    AlgoritmoFloodFillUIterativo(M, iCelulaInicial, jCelulaInicial, FULL);
    AlgoritmoFloodFillVIterativo(M, iCelulaInicial, jCelulaInicial, FULL);
    AlgoritmoFloodFillPIterativo(M, iCelulaInicial, jCelulaInicial, FULL);

    //Copiando o pontosP pro celulas, apenas pra visualizacao do VTK qdo for confinado
    for( i=0; i<M.Nx; i++ )
        for( j=0; j<M.Ny; j++ )
            M.celulas[i][j] = M.pontosP[i][j].tipo;

    if( M.tipoCoord == AXI_CILINDRICO )
        AlgoritmoFloodFillWIterativo(M, iCelulaInicial, jCelulaInicial, FULL);
    return;
}

void PreencheRegiaoParede(MALHA M, int NumFaces, int *Faces)
{
    int faceEsq = -1;
    double xMin = M.xMax + 1, xMax = M.xMin-1, ponto;
    int i, iCelulaInicial, jCelulaInicial;
    BLOCO *bloco;

    for( i=0; i<NumFaces; i++ ) {
        bloco = BuscaBloco( M, Faces[i] );
        if(bloco->direcaoNormal==DIRECAO_X) {
            ponto = M.xMin + M.dx[0]*bloco->iMin;
            if( ponto<xMin ) {
                xMin = ponto;
                faceEsq = i;
            }
            if( ponto>xMax ) {
                xMax = ponto;
                //faceDir = i;
            }
        }
    }

    bloco = BuscaBloco(M, Faces[faceEsq]);
    iCelulaInicial = bloco->iMin + 1;
    jCelulaInicial = (bloco->jMin + bloco->jMax)/2;

    AlgoritmoFloodFillCelulaParedeIterativo(M, iCelulaInicial, jCelulaInicial, BOUNDARY);

    return;
}

// void PreencheRegiaoPoligono(MALHA M, CURVA *C)
// {
//     int i, j, qtdIntersec, trocou, indice, dentro;
//     double x, y, intersecX[100], temp;
//     PONTO *p1, *p2;

//     if( C==NULL || C->pontoInicial==NULL )
//         return;

//     ///Fazendo o preenchimento dos pontos onde estao as velocidades U
//     //Passando por cada reta horizontal
//     for( j=0; j<M.Ny; j++ ) {
//         qtdIntersec = 0;
//         y = M.y[j] + 0.5*M.dy[j];

//         for( p1=C->nodeInicial; p1->prox!=NULL; p1=p1->prox ) {
//             p2 = p1->prox;

//             //Encontrando as interseccoes do poligono com esta reta horizontal
//             if( (p1->y<y && p2->y>=y) || (p2->y<y && p1->y>=y) )
//                 intersecX[qtdIntersec++]= ( ((p1->x - p2->x)*y) + (p2->x*p1->y) - (p1->x*p2->y) )/(p1->y - p2->y);
//         }

//         //Ordenando as interseccoes com um bubble-sort
//         trocou = 1;
//         while(trocou) {
//             trocou = 0;
//             for( i=1; i<qtdIntersec; i++ ) {
//                 if( intersecX[i-1]>intersecX[i] ) {
//                     temp = intersecX[i-1];
//                     intersecX[i-1] = intersecX[i];
//                     intersecX[i] = temp;
//                     trocou = 1;
//                 }
//             }
//         }

//         //Preenchendo os pontos da malha desta linha horizontal
//         dentro = 0;
//         indice = 0;
//         for( i=0; i<=M.Nx; i++ ) {
//             x = M.x[i];

//             if( (indice<qtdIntersec) && (intersecX[indice]<x) ) {
//                 dentro = !dentro;
//                 indice++;
//             }

//             if( dentro )
//                 M.pontosU[i][j].tipo = FULL;
//         }
//     }


//     ///Fazendo o preenchimento dos pontos onde estao as velocidades V
//     //Passando por cada reta horizontal
//     for( j=0; j<=M.Ny; j++ ) {
//         qtdIntersec = 0;
//         y = M.y[j];

//         for( p1=C->nodeInicial; p1->prox!=NULL; p1=p1->prox ) {
//             p2 = p1->prox;

//             //Encontrando as interseccoes do poligono com esta reta horizontal
//             if( (p1->y<y && p2->y>=y) || (p2->y<y && p1->y>=y) )
//                 intersecX[qtdIntersec++]= ( ((p1->x - p2->x)*y) + (p2->x*p1->y) - (p1->x*p2->y) )/(p1->y - p2->y);
//         }

//         //Ordenando as interseccoes com um bubble-sort
//         trocou = 1;
//         while(trocou) {
//             trocou = 0;
//             for( i=1; i<qtdIntersec; i++ ) {
//                 if( intersecX[i-1]>intersecX[i] ) {
//                     temp = intersecX[i-1];
//                     intersecX[i-1] = intersecX[i];
//                     intersecX[i] = temp;
//                     trocou = 1;
//                 }
//             }
//         }

//         //Preenchendo os pontos da malha desta linha horizontal
//         dentro = 0;
//         indice = 0;
//         for( i=0; i<M.Nx; i++ ) {
//             x = M.x[i] + 0.5*M.dx[i];

//             if( (indice<qtdIntersec) && (intersecX[indice]<x) ) {
//                 dentro = !dentro;
//                 indice++;
//             }

//             if( dentro )
//                 M.pontosV[i][j].tipo = FULL;
//         }
//     }


//     ///Fazendo o preenchimento dos pontos onde esta a pressao P
//     //Passando por cada reta horizontal
//     for( j=0; j<M.Ny; j++ ) {
//         qtdIntersec = 0;
//         y = M.y[j] + 0.5*M.dy[j];

//         for( p1=C->nodeInicial; p1->prox!=NULL; p1=p1->prox ) {
//             p2 = p1->prox;

//             //Encontrando as interseccoes do poligono com esta reta horizontal
//             if( (p1->y<y && p2->y>=y) || (p2->y<y && p1->y>=y) )
//                 intersecX[qtdIntersec++]= ( ((p1->x - p2->x)*y) + (p2->x*p1->y) - (p1->x*p2->y) )/(p1->y - p2->y);
//         }

//         //Ordenando as interseccoes com um bubble-sort
//         trocou = 1;
//         while(trocou) {
//             trocou = 0;
//             for( i=1; i<qtdIntersec; i++ ) {
//                 if( intersecX[i-1]>intersecX[i] ) {
//                     temp = intersecX[i-1];
//                     intersecX[i-1] = intersecX[i];
//                     intersecX[i] = temp;
//                     trocou = 1;
//                 }
//             }
//         }

//         //Preenchendo os pontos da malha desta linha horizontal
//         dentro = 0;
//         indice = 0;
//         for( i=0; i<M.Nx; i++ ) {
//             x = M.x[i] + 0.5*M.dx[i];

//             if( (indice<qtdIntersec) && (intersecX[indice]<x) ) {
//                 dentro = !dentro;
//                 indice++;
//             }

//             if( dentro )
//                 M.pontosP[i][j].tipo = FULL;
//         }
//     }

// }

BLOCO *BuscaBloco(MALHA M, int id)
{
    BLOCO *bloco;

    for( bloco=M.primBloco; bloco!=NULL; bloco=bloco->prox ) {
        if( bloco->id==id )
            return bloco;
    }

    return NULL;
}

void InsereBlocoFreeSurface(MALHA *M, TIPO_REGIAO_FREE TipoBloco, double VelX, double VelY, double xMin, double yMin, double xMax, double yMax, double RaioX, double RaioY, double CentroX, double CentroY)
{
    BLOCO *novoBloco;

    novoBloco = (BLOCO *)malloc( sizeof(BLOCO) );
    novoBloco->xMin = xMin;
    novoBloco->yMin = yMin;
    novoBloco->xMax = xMax;
    novoBloco->yMax = yMax;
    novoBloco->centroX = CentroX;
    novoBloco->centroY = CentroY;
    novoBloco->raioX = RaioX;
    novoBloco->raioY = RaioY;
    novoBloco->tipoRegiaoFree = TipoBloco;
    novoBloco->freeVelX = VelX;
    novoBloco->freeVelY = VelY;

    novoBloco->prox = M->blocoLivre;
    M->blocoLivre = novoBloco;
    return;
}

BLOCO *BuscaBlocoFreeSurface(MALHA *M, double X, double Y)
{
    BLOCO *b;

    for( b=M->blocoLivre; b!=NULL; b=b->prox ) {

        if( b->tipoRegiaoFree == TIPO_FREE_ELIPSE ) {
            double valor = 0.0;
            valor += (X - b->centroX)*(X - b->centroX)/(b->raioX*b->raioX);
            valor += (Y - b->centroY)*(Y - b->centroY)/(b->raioY*b->raioY);

            if( valor<=(1.0 + 1e-8) )
                return b;
        }
        else if( b->tipoRegiaoFree == TIPO_FREE_RETANGULO ) {
            if( (X>=b->xMin-(1e-8)) && (X<=b->xMax+(1e-8)) &&
                (Y>=b->yMin-(1e-8)) && (Y<=b->yMax+(1e-8)) ) {
                return b;
            }
        }

    }

    return NULL;
}

void AlgoritmoFloodFillU(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    //printf("FLOW FILL %d %d", i, j);

    if( i>M.Nx || i<0 )
        return;
    if( j>=M.Ny || j<0 )
        return;

    if( M.pontosU[i][j].tipo==Tipo )
        return;
    if( M.pontosU[i][j].tipo==BOUNDARY )
        return;

    M.pontosU[i][j].tipo = Tipo;

    //Chamadas recursivas
    //if( i!=M.Nx && M.pontosU[i+1][j].tipo!=Tipo && M.pontosU[i+1][j].tipo!=BOUNDARY )
    AlgoritmoFloodFillU(M, i+1, j, Tipo);
    //if( i!=0 && M.pontosU[i-1][j].tipo!=Tipo && M.pontosU[i-1][j].tipo!=BOUNDARY )
    AlgoritmoFloodFillU(M, i-1, j, Tipo);
    //if( j!=M.Ny-1 && M.pontosU[i][j+1].tipo!=Tipo && M.pontosU[i][j+1].tipo!=BOUNDARY )
    if( j!=M.Ny-1 && M.pontosV[i][j+1].tipo!=BOUNDARY )
        AlgoritmoFloodFillU(M, i, j+1, Tipo);
    //if( j!=0 && M.pontosU[i][j-1].tipo!=Tipo && M.pontosU[i][j-1].tipo!=BOUNDARY )
    if( j!=0 && M.pontosV[i][j].tipo!=BOUNDARY )
        AlgoritmoFloodFillU(M, i, j-1, Tipo);


    return;
}

void AlgoritmoFloodFillV(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    if( i>=M.Nx || i<0 )
        return;
    if( j>M.Ny || j<0 )
        return;

    if( M.pontosV[i][j].tipo==Tipo )
        return;
    if( M.pontosV[i][j].tipo==BOUNDARY )
        return;

    M.pontosV[i][j].tipo = Tipo;

    //Chamadas recursivas
    //if( i!=M.Nx && M.pontosU[i+1][j].tipo!=Tipo && M.pontosU[i+1][j].tipo!=BOUNDARY )
    if( i!=M.Nx-1 && M.pontosU[i+1][j].tipo!=BOUNDARY )
        AlgoritmoFloodFillV(M, i+1, j, Tipo);
    //if( i!=0 && M.pontosU[i-1][j].tipo!=Tipo && M.pontosU[i-1][j].tipo!=BOUNDARY )
    if( i!=0 && M.pontosU[i][j].tipo!=BOUNDARY )
        AlgoritmoFloodFillV(M, i-1, j, Tipo);
    //if( j!=M.Ny-1 && M.pontosU[i][j+1].tipo!=Tipo && M.pontosU[i][j+1].tipo!=BOUNDARY )
    AlgoritmoFloodFillV(M, i, j+1, Tipo);
    //if( j!=0 && M.pontosU[i][j-1].tipo!=Tipo && M.pontosU[i][j-1].tipo!=BOUNDARY )
    AlgoritmoFloodFillV(M, i, j-1, Tipo);

    return;
}

void AlgoritmoFloodFillP(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    //printf("FLOW FILL %d %d", i, j);

    if( i>=M.Nx || i<0 )
        return;
    if( j>=M.Ny || j<0 )
        return;

    if( M.pontosP[i][j].tipo==Tipo )
        return;
    if( M.pontosP[i][j].tipo==BOUNDARY )
        return;

    M.pontosP[i][j].tipo = Tipo;

    //Chamadas recursivas
    if( i!=M.Nx-1 && M.pontosU[i+1][j].tipo!=BOUNDARY )
        AlgoritmoFloodFillP(M, i+1, j, Tipo);
    if( i!=0 && M.pontosU[i][j].tipo!=BOUNDARY )
        AlgoritmoFloodFillP(M, i-1, j, Tipo);
    if( j!=M.Ny-1 && M.pontosV[i][j+1].tipo!=BOUNDARY )
        AlgoritmoFloodFillP(M, i, j+1, Tipo);
    if( j!=0 && M.pontosV[i][j].tipo!=BOUNDARY )
        AlgoritmoFloodFillP(M, i, j-1, Tipo);


    return;
}

void AlgoritmoFloodFillUIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    PILHA pilha;
    int continua;

    pilha.topo = NULL;

    //Insere o elemento inicial na pilha
    PilhaPush(&pilha, i, j);

    while( pilha.topo!=NULL ) {
        PilhaPop(&pilha, &i, &j);

        if( i>M.Nx || i<0 )
            continue;
        if( j>=M.Ny || j<0 )
            continue;

        if( M.pontosU[i][j].tipo==Tipo )
            continue;
        if( M.pontosU[i][j].tipo==BOUNDARY )
            continue;

        M.pontosU[i][j].tipo = Tipo;


        PilhaPush(&pilha, i+1, j);
        PilhaPush(&pilha, i-1, j);

        //Verifica se continua pra cima
        continua = (i!=M.Nx) ? (M.pontosV[i][j+1].tipo!=BOUNDARY) : 1;
        continua = (continua && i!=0) ? (M.pontosV[i-1][j+1].tipo!=BOUNDARY) : 0;
        if( continua && j!=M.Ny-1 )
            PilhaPush(&pilha, i, j+1);

        //Verifica se continua pra baixo
        continua = (i!=M.Nx) ? (M.pontosV[i][j].tipo!=BOUNDARY) : 1;
        continua = (continua && i!=0) ? (M.pontosV[i-1][j].tipo!=BOUNDARY) : 0;
        if( j!=0 && continua )
            PilhaPush(&pilha, i, j-1);
    }


    return;
}

void AlgoritmoFloodFillVIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    PILHA pilha;

    pilha.topo = NULL;

    //Insere o elemento inicial na pilha
    PilhaPush(&pilha, i, j);

    while( pilha.topo!=NULL ) {
        PilhaPop(&pilha, &i, &j);

        if( i>=M.Nx || i<0 )
            continue;
        if( j>M.Ny || j<0 )
            continue;

        if( M.pontosV[i][j].tipo==Tipo )
            continue;
        if( M.pontosV[i][j].tipo==BOUNDARY )
            continue;

        M.pontosV[i][j].tipo = Tipo;

        //Chamadas recursivas
        if( i!=M.Nx-1 && M.pontosU[i+1][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i+1, j);
        if( i!=0 && M.pontosU[i][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i-1, j);
        PilhaPush(&pilha, i, j+1);
        PilhaPush(&pilha, i, j-1);
    }

    return;
}

void AlgoritmoFloodFillPIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    PILHA pilha;

    pilha.topo = NULL;

    //Insere o elemento inicial na pilha
    PilhaPush(&pilha, i, j);

    while( pilha.topo!=NULL ) {
        PilhaPop(&pilha, &i, &j);

        if( i>=M.Nx || i<0 )
            continue;
        if( j>=M.Ny || j<0 )
            continue;

        if( M.pontosP[i][j].tipo==Tipo )
            continue;
        if( M.pontosP[i][j].tipo==SURFACE )
            continue;
        if( M.pontosP[i][j].tipo==BOUNDARY )
            continue;

        M.pontosP[i][j].tipo = Tipo;

        //Chamadas recursivas
        if( i!=M.Nx-1 && M.pontosU[i+1][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i+1, j);
        if( i!=0 && M.pontosU[i][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i-1, j);
        if( j!=M.Ny-1 && M.pontosV[i][j+1].tipo!=BOUNDARY )
            PilhaPush(&pilha, i, j+1);
        if( j!=0 && M.pontosV[i][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i, j-1);
    }

    return;
}

void AlgoritmoFloodFillCelulaParedeIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    PILHA pilha;

    pilha.topo = NULL;



    //Insere o elemento inicial na pilha
    PilhaPush(&pilha, i, j);

    while( pilha.topo!=NULL ) {
        PilhaPop(&pilha, &i, &j);

        if( i>=M.Nx || i<0 )
            continue;
        if( j>=M.Ny || j<0 )
            continue;

        if( M.celulasParede[i][j]==Tipo )
            continue;
        if( M.celulasParede[i][j]==SURFACE )
            continue;
        if( M.celulasParede[i][j]==BOUNDARY )
            continue;

        M.celulasParede[i][j] = Tipo;

        //Chamadas recursivas
        if( i!=M.Nx-1 && M.pontosU[i+1][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i+1, j);
        if( i!=0 && M.pontosU[i][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i-1, j);
        if( j!=M.Ny-1 && M.pontosV[i][j+1].tipo!=BOUNDARY )
            PilhaPush(&pilha, i, j+1);
        if( j!=0 && M.pontosV[i][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i, j-1);
    }

    return;
}

void AlgoritmoFloodFillCelulaIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    PILHA pilha;

    pilha.topo = NULL;

    //Insere o elemento inicial na pilha
    PilhaPush(&pilha, i, j);

    while( pilha.topo!=NULL ) {
        PilhaPop(&pilha, &i, &j);

        if( i>=M.Nx || i<0 )
            continue;
        if( j>=M.Ny || j<0 )
            continue;

        if( M.celulas[i][j]==Tipo )
            continue;
        if( M.celulas[i][j]==BOUNDARY )
            continue;
        if( M.celulas[i][j]==SURFACE )
            continue;

        M.celulas[i][j] = Tipo;

        //Chamadas recursivas
        if( i!=M.Nx-1 && M.celulas[i+1][j]!=SURFACE )
            PilhaPush(&pilha, i+1, j);
        if( i!=0 && M.celulas[i][j]!=SURFACE )
            PilhaPush(&pilha, i-1, j);
        if( j!=M.Ny-1 && M.celulas[i][j+1]!=SURFACE )
            PilhaPush(&pilha, i, j+1);
        if( j!=0 && M.celulas[i][j]!=SURFACE )
            PilhaPush(&pilha, i, j-1);
    }

    return;
}

void AlgoritmoFloodFillWIterativo(MALHA M, int i, int j, TIPO_CELULA Tipo)
{
    PILHA pilha;

    pilha.topo = NULL;

    //Insere o elemento inicial na pilha
    PilhaPush(&pilha, i, j);

    while( pilha.topo!=NULL ) {
        PilhaPop(&pilha, &i, &j);

        if( i>=M.Nx || i<0 )
            continue;
        if( j>=M.Ny || j<0 )
            continue;

        if( M.pontosW[i][j].tipo==Tipo )
            continue;
        if( M.pontosW[i][j].tipo==BOUNDARY )
            continue;

        M.pontosW[i][j].tipo = Tipo;

        //Chamadas recursivas
        if( i!=M.Nx-1 && M.pontosU[i+1][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i+1, j);
        if( i!=0 && M.pontosU[i][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i-1, j);
        if( j!=M.Ny-1 && M.pontosV[i][j+1].tipo!=BOUNDARY )
            PilhaPush(&pilha, i, j+1);
        if( j!=0 && M.pontosV[i][j].tipo!=BOUNDARY )
            PilhaPush(&pilha, i, j-1);
    }

    return;
}

//VERSAO NOVA!
//double AlgoritmoBisseccao(int Grau, double *Polinomio, double X1, double X2)
//{
//    double meio, pol1, polMeio;
//    static double criterioParada = 1e-13;
//    static int chamadas = 0;
//
//    chamadas++;
//
//    if( fabs(X2-X1)<criterioParada )
//        return 0.5*(X2 + X1);
//
//    if( !(chamadas%2000) )
//        criterioParada *= 10.0;
//
//    meio = 0.5*(X2 + X1);
//    pol1 = ValorPolinomio(X1, Grau, Polinomio);
//    polMeio = ValorPolinomio(meio, Grau, Polinomio);
//
//    if( pol1*polMeio<=0 )
//        return AlgoritmoBisseccao(Grau, Polinomio, X1, meio);
//    else
//        return AlgoritmoBisseccao(Grau, Polinomio, meio, X2);
//}

double AlgoritmoBisseccao(int Grau, double *Polinomio, double X1, double X2)
{
    double meio, pol1, polMeio;

    if( fabs(X2-X1)<1e-14 )
        return 0.5*(X2 + X1);

    meio = 0.5*(X2 + X1);
    pol1 = ValorPolinomio(X1, Grau, Polinomio);
    polMeio = ValorPolinomio(meio, Grau, Polinomio);

    if( pol1*polMeio<=0 )
        return AlgoritmoBisseccao(Grau, Polinomio, X1, meio);
    else
        return AlgoritmoBisseccao(Grau, Polinomio, meio, X2);
}

double ValorPolinomio(double X, int Grau, double *Pol)
{
    double soma;
    int i;

    soma = Pol[Grau];
    for(i=Grau-1; i>=0; i--)
        soma = Pol[i] + X*soma;

    return soma;
}

void EncontraCelulaDaParticula(MALHA *M, double X, double Y, int *iResult, int *jResult)
{
    int i, j;

    for( i=0; i<M->Nx; i++ ) {
        if( X < M->x[i+1] )
            break;
    }

    for( j=0; j<M->Ny; j++ ) {
        if( Y < M->y[j+1] )
            break;
    }


    if( i>=M->Nx || j>=M->Ny ) {
//        ImprimeInterfaceVTK(*M, 1000000);
//        DesenhaMalhaVTK(*M, 0);
        DeuProblema("\nProblema (Malha.c): EncontrarCelulaDaParticula %.15lf %.15lf %.15lf\n", X, Y, M->y[M->Ny]);
    }

    *iResult = i;
    *jResult = j;
    return;
}



// int EmbaracamentoForcado_Raycast(MALHA *M, double Distancia)
// {
//     CURVA *c1, *c2;
//     PONTO *p1, *p2;
//     double xHit, yHit;

//     /// Percorrendo cada ponto p1 e fazendo um raycast a partir dele na direcao normal
//     for( c1=M->interface.curvaInicial; c1; c1=c1->proxCurva ) {
//         LOOP_CURVA(c1) {
//             p1 = loop_stuff.ponto;

//             /// Calculando o vetor normal no ponto p1
//             double nx, ny;
//             VetorNormal(c1, p1, &nx, &ny);

//             /// Percorrendo cada ponto p2 e verificando se o raycast intersecta o segmento de reta a frente de p2
//             for( c2=M->interface.curvaInicial; c2; c2=c2->proxCurva ) {
//                 LOOP_CURVA(c2) {
//                     p2 = loop_stuff.ponto;

//                     /// Faz o raycast e ve se cruza a interface
//                     PONTO pontoRay;
//                     pontoRay.x = p1->x + Distancia*nx;
//                     pontoRay.y = p1->y + Distancia*ny;
//                     char encontrou = PossuiInterseccaoSegmentos(p1, &pontoRay, p2, p2->prox, &xHit, &yHit);



//                     if( encontrou && (c1!=c2) ) {
//                         PrintDebug("MUDANCA TOPOLOGICA: UNIAO DUAS LISTAS\n\n");

//                         EmbaracamentoForcado_DuasListas(M, Distancia, p1, p2, c1, c2);

//                         // Armazena todos os nodes atuais
//                         CURVA *curvaTemp;
//                         int qtdCurvas = 0;
//                         for( curvaTemp=M->interface.curvaInicial; curvaTemp; curvaTemp=curvaTemp->proxCurva )
//                             qtdCurvas++;
//                         PONTO *pontosAnt = (PONTO *)malloc( qtdCurvas*sizeof(PONTO) );

//                         int i = 0;
//                         for( curvaTemp=M->interface.curvaInicial; curvaTemp; curvaTemp=curvaTemp->proxCurva )
//                             pontosAnt[i++] = *(curvaTemp->pontoInicial);


//                         PONTO nodes[2];
//                         RealizaDesembaracamento(&(M->interface), nodes);

//                         // Remove a curva nova que sobrou no meio
//                         int encontrouNode = 0;
//                         for( curvaTemp=M->interface.curvaInicial; curvaTemp; curvaTemp=curvaTemp->proxCurva ) {
//                             PONTO *node = curvaTemp->nodeInicial;
//                             for( i=0; i<qtdCurvas; i++ ) {
//                                 if( pontosAnt[i].x==node->x && pontosAnt[i].y==node->y )
//                                     encontrouNode = 1;
//                             }

//                             // Remove curva
//                             if( encontrouNode==0 ) {
//                                 CURVA *proxCurva = curvaTemp->proxCurva;
//                                 if( curvaTemp->curvaAnt )
//                                     curvaTemp->curvaAnt->proxCurva = proxCurva;
//                                 if( proxCurva )
//                                     proxCurva->curvaAnt = curvaTemp->curvaAnt;
//                                 if( M->interface.curvaInicial==curvaTemp )
//                                     M->interface.curvaInicial = proxCurva;
//                                 LiberaMemoriaCurva(curvaTemp);
//                             }
//                         }

//                         PrintDebug("FINAL: MUDANCA TOPOLOGICA: UNIAO DUAS LISTAS\n\n");
//                         return 1;
//                     }
//                 }
//             }


//         }
//     }


//     return 0;
// }

// void EmbaracamentoForcado_DuasListas(MALHA *M, double Distancia, PONTO *p1C1, PONTO *p1C2, CURVA *c1, CURVA *c2)
// {
//     PONTO *p2C1, *p2C2;
//     double xHit, yHit;
//     int continuar = 1;

//     PONTO *ultP2C1 = NULL, *ultP2C2 = NULL;

//     /// Continuo percorrendo a lista C1 e fazendo raycast ate acabar essa regiao de pontos muito grudadas
//     p2C1 = p1C1;
//     while( continuar ) {
//         p2C1 = p2C1->prox ? p2C1->prox : p2C1->pontoProxCurva;

//         /// Calculando o vetor normal no ponto p1
//         double nx, ny;
//         VetorNormal(c1, p2C1, &nx, &ny);

//         /// Percorrendo cada ponto p2 e verificando se o raycast intersecta o segmento de reta a frente de p2
//         for( p2C2=c2->nodeInicial; p2C2->prox; p2C2=p2C2->prox ) {

//             /// Faz o raycast e ve se cruza a interface
//             PONTO pontoRay;
//             pontoRay.x = p2C1->x + Distancia*nx;
//             pontoRay.y = p2C1->y + Distancia*ny;
//             continuar = PossuiInterseccaoSegmentos(p2C1, &pontoRay, p2C2, p2C2->prox, &xHit, &yHit);

//             if( continuar ) {
//                 ultP2C1 = p2C1;
//                 ultP2C2 = p2C2;
//                 break;
//             }
//         }
//     }

//     p2C1 = ultP2C1;
//     p2C2 = ultP2C2;

//     if( p2C1==NULL || p2C2==NULL )
//         DeuProblema("PROBLEMA: EMBARACAMENTO FORCADO DUAS LISTAS RAYCAST\n\n");

//     /// Agora realiza a troca dos pontos
//     double tempX, tempY;
//     tempX = p1C1->x; tempY = p1C1->y;
//     p1C1->x = p1C2->x; p1C1->y = p1C2->y;
//     p1C2->x = tempX; p1C2->y = tempY;

//     tempX = p2C1->x; tempY = p2C1->y;
//     p2C1->x = p2C2->x; p2C1->y = p2C2->y;
//     p2C2->x = tempX; p2C2->y = tempY;

//     return;
// }

// /// Funcoes para forcar o embaracamento da interface dependendo da classificacao das celulas na malha
// int RealizaEmbaracamentoForcado(MALHA *M)
// {
//     int i, j, listaAtual, qtdCelulasMudanca, embaracou, valorRetorno;
//     int comecouGrupo, adicionouLista, minI, maxI, minJ, maxJ, condicaoI, condicaoJ;
//     PONTO *ponto;
//     CURVA *curva;
//     LISTA_MUDANCA_TOPOL lista[2];
//     PONTO_MUDANCA_TOPOL *pontoMudanca, *proxPontoMudanca;

//     embaracou = 0;
//     qtdCelulasMudanca = 0;

//     //Transformando todas celulas FULL com particulas em MUDANCA_TOPOL
//     for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//         for( ponto=curva->nodeInicial; ponto!=NULL; ponto=ponto->prox ) {
//             EncontraCelulaDaParticula(M, ponto->x, ponto->y, &i, &j);

//             if( M->celulas[i][j]==FULL ) {
//                 M->celulas[i][j] = MUDANCA_TOPOL;
//                 qtdCelulasMudanca++;
//             }
//         }
//     }



//     RemoveCelulasIndesejadasMudancaTopologicaPasso1(M, &qtdCelulasMudanca);
//     RemoveCelulasIndesejadasMudancaTopologicaPasso2(M, &qtdCelulasMudanca);

//     //o valor ZERO indica que nao realizou o embaracamento forcado
//     if( !qtdCelulasMudanca )
//         return 0;


// //    static int a = 0;
// //    a++;
// //    if( a==2 ) {
// //        DesenhaMalhaVTK(*M, 0);
// //        ImprimeInterfaceVTK(*M, 10000);
// //        //DeuProblema("OK!!");
// //    }

//     valorRetorno = 0;
//     do {
//         minI = minJ = 100000;
//         maxI = maxJ =  -100000;

//         //SEPARANDO os pontos em duas listas que vao ser embaracadas uma com a outra
//         lista[0].pontoInicial = NULL;
//         lista[1].pontoInicial = NULL;
//         listaAtual = 0;
//         comecouGrupo = 0;
//         embaracou = 0;
//         for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {

//             adicionouLista = 0;
//             for( ponto=curva->nodeInicial; ponto!=NULL; ponto=ponto->prox ) {
//                 EncontraCelulaDaParticula(M, ponto->x, ponto->y, &i, &j);

//                 if( M->celulas[i][j]!=MUDANCA_TOPOL )
//                     continue;


//                 condicaoI = (i>minI-3) && (i<maxI+3);
//                 condicaoJ = (j>minJ-3) && (j<maxJ+3);

//                 if( listaAtual>1 )
//                     DeuProblema("PROBLEMA EMBARACAMENTO FORCADO!\n");

//                 if( !comecouGrupo ) {
//                     AdicionaNaListaDeMudancaTopol(&lista[listaAtual], ponto);
//                     AtualizaMinimoMaximo(i, j, &minI, &minJ, &maxI, &maxJ);
//                     comecouGrupo = 1;
//                     //PrintDebug("Comecou Grupo: %d %d\n", i, j);

//                     adicionouLista = 1;
//                 }
//                 else if( condicaoI && condicaoJ ) {
//                     AdicionaNaListaDeMudancaTopol(&lista[listaAtual], ponto);
//                     AtualizaMinimoMaximo(i, j, &minI, &minJ, &maxI, &maxJ);
//                     //PrintDebug("Continuou Grupo: %d %d\n", i, j);

//                     adicionouLista = 1;
//                 }


//             }

//             if( adicionouLista )
//                 listaAtual++;

//             if( listaAtual==2 )
//                 break;
//         }

//         if( comecouGrupo && listaAtual==2 ) {
//             embaracou = 1;
//             valorRetorno = 1;
//             ///Faz o esquema de duas listas
//             EmbaracamentoCasoDuasListas(M, lista);
//         }
//         else if( comecouGrupo && ( ((maxI-minI)>1) || ((maxJ-minJ)>1) ) ) {
//             embaracou = 1;
//             valorRetorno = 1;

//             ///Faz o esquema de uma unica lista
//             //EmbaracamentoCasoUmaLista(M, lista[0]);
//             embaracou = 0; // tirar essa linha
//         }


//         ///Liberando memoria das listas
//         for( i=0; i<2; i++ ) {
//             pontoMudanca = lista[i].pontoInicial;
//             while( pontoMudanca!=NULL ) {
//                 proxPontoMudanca = pontoMudanca->prox;
//                 free(pontoMudanca);
//                 pontoMudanca = proxPontoMudanca;
//             }
//         }

//     } while( embaracou );


//     ///Remove todas as flags MUDANCA_TOPOL da malha
//     for( i=0; i<M->Nx; i++ ) {
//         for( j=0; j<M->Ny; j++ ) {
//             if( M->celulas[i][j]==MUDANCA_TOPOL )
//                 M->celulas[i][j] = FULL;
//         }
//     }


//     //o valor UM indica que realizou o embaracamento forcado
//     return valorRetorno;
// }

// int RealizaMudancaTopologica_Uniao(MALHA *M)
// {
//     CURVA *c, *curva1 = NULL, *curva2 = NULL;
//     PONTO *p1_curva1 = NULL, *p2_curva1 = NULL, *p1_curva2 = NULL, *p2_curva2 = NULL, *p;
//     int qtdPontosFull_curva1, qtdPontosFull_curva2;

//     /// Procurando a primeira curva com uma regiao pra ser unida
//     for( c=M->interface.curvaInicial; c; c=c->proxCurva ) {
//         if( curva1 )
//             break;

//         qtdPontosFull_curva1 = 1;

//         for( p=c->nodeInicial; p; p=p->prox ) {
//             int i, j;
//             EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);


//             if( (M->celulas[i][j]==FULL) /*&& (fabs(p->x)>1e-6)*/ && (p->y>1e-6) ) {
//                 if( p1_curva1==NULL )
//                     p1_curva1 = p;
//                 else
//                     qtdPontosFull_curva1++;
//             }
//             else if( p1_curva1 ) {
//                 if( qtdPontosFull_curva1>20 ) {
//                     p2_curva1 = p;
//                     curva1 = c;
//                     // ImprimeInterfaceVTK(*M, 800000);
//                     // PrintDebug("CURVA 1... %d: %lf %lf .... %lf %lf\n", qtdPontosFull_curva1, p1_curva1->x, p1_curva1->y, p2_curva1->x, p2_curva1->y);
//                     break;
//                 }
//                 else {
//                     p1_curva1 = NULL;
//                     qtdPontosFull_curva1 = 1;
//                 }
//             }

            
//         }
//     }

//     /// Se nao tiver encontrado nenhuma curva com regiao FULL, nem precisa continuar
//     if( !curva1 )
//         return 0;

//     PrintDebug("Encontrou curva 1... %lf %lf ... %lf %lf\n", p1_curva1->x, p1_curva1->y, p2_curva1->x, p2_curva1->y);

//     // Pegando o ponto "central" da regiao encontrada na curva 1
//     PONTO *pCentral_curva1 = p1_curva1, *pCentral_curva2;
//     int contador;
//     for( contador=qtdPontosFull_curva1/2; contador>=0; contador--)
//         pCentral_curva1 = pCentral_curva1->prox;

//     // PrintDebug("CENTRAL C1: %lf %lf\n", pCentral_curva1->x, pCentral_curva1->y);

//     /// Procurando a segunda curva com uma regiao pra ser unida com a primeira
//     for( c=M->interface.curvaInicial; c; c=c->proxCurva ) {
//         if( c==curva1 )
//             continue;

//         // PrintDebug("entrou curva\n");

//         qtdPontosFull_curva2 = 1;
//         if( curva2 )
//             break;

//         for( p=c->nodeInicial; p; p=p->prox ) {
//             int i, j;
//             EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);

//             if( (M->celulas[i][j]==FULL) /*&& (fabs(p->x)>1e-6)*/ && (p->y>1e-6) ) {
//                 if( p1_curva2==NULL )
//                     p1_curva2 = p;
//                 else
//                     qtdPontosFull_curva2++;
//             }
//             else if( p1_curva2 ) {
//                 if( (qtdPontosFull_curva2>15) ) {
//                     pCentral_curva2 = p1_curva2;
//                     for( contador=qtdPontosFull_curva2/2; contador>=0; contador--)
//                         pCentral_curva2 = pCentral_curva2->prox;

//                     double dist_x = (pCentral_curva1->x - pCentral_curva2->x);
//                     double dist_y = (pCentral_curva1->y - pCentral_curva2->y);
//                     double distancia_centros = sqrt( dist_x*dist_x + dist_y*dist_y );

//                     // PrintDebug("HAHAHA %lf %lf %lf %lf\n", p1_curva2->x, p1_curva2->y, distancia_centros, 2.0*M->minDxDy);
//                     if( distancia_centros<3.0*M->minDxDy ) {
//                         p2_curva2 = p;
//                         curva2 = c;
//                         // ImprimeInterfaceVTK(*M, 800000);
//                         // PrintDebug("CURVA 2... %d: %lf %lf .... %lf %lf\n", qtdPontosFull_curva2, p1_curva2->x, p1_curva2->y, p2_curva2->x, p2_curva2->y);
//                         // PrintDebug("DISTANCIA CENTROS: %lf\n", distancia_centros);
//                         break;
//                     }
//                     else {
//                         p1_curva2 = NULL;
//                         qtdPontosFull_curva2 = 1;
//                         curva2 = NULL;
//                         p2_curva2 = NULL;
//                     }
//                 }
//                 else {
//                     // printf("QUANTIDADE: %d\n", qtdPontosFull_curva2);
//                     p1_curva2 = NULL;
//                     qtdPontosFull_curva2 = 1;
//                 }
//             }

            
//         }
//     }

//     // PrintDebug("Encontrou curva 2... %lf %lf ... %lf %lf\n", p1_curva2->x, p1_curva2->y, p2_curva2->x, p2_curva2->y);

//     /// Se nao tiver encontrado uma segunda curva tambem, pode parar aqui
//     if( !curva2 )
//         return 0;

//     /// Se chegou ate aqui, eh pq encontrou curva1 e curva2 pra unir. vamos fazer a uniao
//     PONTO *pA = p1_curva1;
//     double distancia1 = (pA->x - p1_curva2->x)*(pA->x - p1_curva2->x) + (pA->y - p1_curva2->y)*(pA->y - p1_curva2->y);
//     double distancia2 = (pA->x - p2_curva2->x)*(pA->x - p2_curva2->x) + (pA->y - p2_curva2->y)*(pA->y - p2_curva2->y);
//     PONTO *pB = (distancia1 < distancia2) ? p1_curva2 : p2_curva2;
//     // PrintDebug("pA pB: %lf %lf ... %lf %lf\n", pA->x, pA->y, pB->x, pB->y);
//     pA->prox = pB;
//     pB->ant = pA;

//     PONTO *ultimo = UltimoPonto(pB);
//     ultimo->prox = ultimo->pontoProxCurva;
//     ultimo->prox->ant = ultimo;
//     ultimo->pontoProxCurva = NULL;

//     pA = p2_curva1;
//     pB = ( pB==p1_curva2 ) ? p2_curva2 : p1_curva2;
//     pB->prox = pA;
//     pA->ant = pB;

//     // Removendo a curva 2 da interface
//     if( M->interface.curvaInicial==curva2 )
//         M->interface.curvaInicial = curva2->proxCurva;
//     else
//         curva2->curvaAnt->proxCurva = curva2->proxCurva;

//     // ImprimeInterfaceVTK(*M, 900000);

//     return 1;
// }

int RealizaEmbaracamentoForcado_Quebra2(MALHA *M, double DistanciaMax, double DistanciaMin, PONTO **pontoRemover, int (*CriterioOpcional)(MALHA*, PONTO*, PONTO*, int) )
{
    PONTO *p1, *p2, *pProx, *p1Escolhido = NULL, *p2Escolhido = NULL;
    CURVA *curva;
    double comprimentoCurva, distancia, distanciaFrente, distanciaTras, distanciaEscolhido;




    ///Cada curva...
    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        distanciaEscolhido = 1e+10;
        p1Escolhido = NULL;
        p2Escolhido = NULL;

        /// Medindo o comprimento total da curva
        comprimentoCurva = 0.0;
        LOOP_CURVA(curva) {
            p1 = loop_stuff.ponto;
            p2 = p1->prox;
            distancia = sqrt( (p2->x - p1->x)*(p2->x - p1->x) + (p2->y - p1->y)*(p2->y - p1->y) );
            comprimentoCurva += distancia;
        }

//        PrintDebug("COMPRIMENTO TOTAL: %lf\n", comprimentoCurva);

        /// Encontrando o primeiro ponto da regiao que sera removida
        LOOP_CURVA(curva) {
            p1 = loop_stuff.ponto;

            distanciaFrente = 0.0;
            distanciaTras = comprimentoCurva;
            for( p2=p1; distanciaTras>=0.0; p2=p2->prox ) {
                pProx = p2->prox;


                /// Medindo a distancia entre p2 e p1 seguindo a curva
                distancia = sqrt( (p2->x - pProx->x)*(p2->x - pProx->x) + (p2->y - pProx->y)*(p2->y - pProx->y) );
                distanciaFrente += distancia;
                distanciaTras -= distancia;

                /// Medindo a distancia entre p2 e p1 em linha reta
                distancia = sqrt( (p2->x - p1->x)*(p2->x - p1->x) + (p2->y - p1->y)*(p2->y - p1->y) );

                // PrintDebug("BB: %lf %lf\n", distancia, DistanciaMax);

                if( (p1!=p2) && (distancia<DistanciaMax) && (distancia<distanciaEscolhido) && (distanciaFrente>DistanciaMin) && (distanciaTras>DistanciaMin) ) {
                    /// Criterios opcionais desta simulacao
                    if( (CriterioOpcional!=NULL) && CriterioOpcional(M, p1, p2, 0) ) {
                        p2 = pProx;
                        continue;
                    }

                    distanciaEscolhido = distancia;
                    p1Escolhido = p1;
                    p2Escolhido = p2;
                }
            }
        }


        /// Se nao tiver encontrado nenhum ponto pra quebrar, apenas vai pra proxima curva
        PrintDebug("AAA: %p\n", p1Escolhido);
        if( p1Escolhido==NULL )
            continue;

//        ImprimeInterfaceVTK(*M, 10000000);
        PrintDebug("p1Escolhido: %lf %lf\n", p1Escolhido->x, p1Escolhido->y);
        PrintDebug("p2Escolhido: %lf %lf\n", p2Escolhido->x, p2Escolhido->y);

        p1 = p1Escolhido;
        p2 = p2Escolhido;
        PONTO p1Copia = *p1, p2Copia = *p2;

        int qtdPontosEmbaraca = 15;
        int i;
        PONTO *A1 = p1, *A2 = p2, *B1 = p1, *B2 = p2;
        PONTO *A1delete, *A2delete, *B1delete, *B2delete;
        int removeuNode = 0;
        for( i=0; i<qtdPontosEmbaraca; i++ ) {
            A1delete = A1;
            A2delete = A2;
            B1delete = B1;
            B2delete = B2;

            A1 = A1->ant;
            A2 = A2->prox;

            B1 = B1->prox;
            B2 = B2->ant;

            if( A1delete==curva->pontoInicial || A2delete==curva->pontoInicial
               || B1delete==curva->pontoInicial || B2delete==curva->pontoInicial )
                removeuNode = 1;

            /// === Precisa liberar essa memoria, seria o correto
            free(A1delete);
            free(A2delete);
            if( B1delete!=A1delete )
                free(B1delete);
            if( B2delete!=A2delete )
                free(B2delete);
        }

        if( A1->prox==NULL )
            DeuProblema("1. Funcao Quebra: A1 era o ponto final");
        if( A2->ant==NULL )
            DeuProblema("2. Funcao Quebra: A2 era o ponto inicial");

        A1->prox = A2;
        A2->ant = A1;

        /// === Se o node foi removido, precisa recolocar o outro node
        if( removeuNode )
            DeuProblema("FALTA IMPLEMENTAR: tentou remover node na funcao de quebra...\n\n");

        CURVA *nova_curva = (CURVA *)malloc( sizeof(CURVA) );
        nova_curva->regiaoEsq = curva->regiaoEsq;
        nova_curva->regiaoDir = curva->regiaoDir;
        nova_curva->proxCurva = NULL;
        nova_curva->curvaAnt = NULL;

        B2->prox = B1;
        B1->ant = B2;
        nova_curva->pontoInicial = B1;
        B1->isNode = 1;
        AdicionaCurvaNaInterface(&(M->interface), nova_curva);
        
        if( CriterioOpcional )
            CriterioOpcional(M, &p1Copia, &p2Copia, 1);
        return 1;
    }


    return 0;
}

// int CriterioOpcional_Quebra(MALHA *M, PONTO *p1, PONTO *p2, int Opcao)
// {
//     static int quebrouBaixo = 0;
//     static int quebrouCima = 0;

//     if( Opcao ) {
//         PrintDebug("quebrou no criterio: %.15lf %.15lf\n", p1->x, p1->y);


//         if(p1->y<0.0)
//             quebrouBaixo = 1;
//         else
//             quebrouCima = 1;
//         return 0;
//     }


//     // No caso da colisao de gotas axissim, vou querer apenas quebras perto do eixo Y
// //    double vec1x = p2->x - p1->x;
// //    double vec1y = p2->y - p1->y;
// //    double vec2x = 1.0;
// //    double vec2y = 0.0;
// //    if( vec1x<0.0 ) {
// //        vec1x *= -1.0;
// //        vec1y *= -1.0;
// //    }

//     /// Criterio da curva front-tracking: Interessado apenas em quebrar "horizontais"
// //    double produto = (vec1x)*(vec2x) + (vec1y)*(vec2y);
// //    double norm1 = sqrt(vec1x*vec1x + vec1y*vec1y);
// //    double norm2 = sqrt(vec2x*vec2x + vec2y*vec2y);
// //    double angulo = (180.0/M_PI)*acos(produto/(norm1*norm2));
// //    if( angulo>15 )
// //        return 1;



//     /// Criterio de malha: interessado em quebrar horizontais
//     int celulasNaoEmpty = 1, k, i, j;
//     EncontraCelulaDaParticula(M, p1->x, p1->y, &i, &j);
//     int iTras, iFrente;
//     for( k=1; k<10; k++ ) {
//         iTras = i-k;
//         iFrente = i+k;
//         if( (iTras>=0) && (iTras<M->Nx) && (M->celulas[iTras][j]!=EMPTY) )
//             celulasNaoEmpty++;
//         if( (iFrente>=0) && (iFrente<M->Nx) && (M->celulas[iFrente][j]!=EMPTY) )
//             celulasNaoEmpty++;
//     }
//     if( celulasNaoEmpty>2 )
//         return 1;



//     if( quebrouBaixo && (p1->y<0.0) )
//         return 1;
//     if( quebrouCima && (p1->y>0.0) )
//         return 1;


// //    PrintDebug("NAO EMPTY: %d %d\n", celulasNaoEmpty, Opcao);
// //    PrintDebug("ANGULO: %lf\n", angulo);
// //    PrintDebug("Ponto: %.15lf %.15lf\n", p1->x, p1->y);



//     return 0;
// }

// void RemoveCurvaDestePonto(MALHA *M, PONTO *pRemover)
// {
//     CURVA *c;
//     PONTO *p;

//     for( c=M->interface.curvaInicial; c; c=c->proxCurva ) {
//         for( p=c->nodeInicial; p; p=p->prox ) {
//             if( p==pRemover ) {

//                 if( c->curvaAnt )
//                     c->curvaAnt->proxCurva = c->proxCurva;
//                 if( c->proxCurva )
//                     c->proxCurva->curvaAnt = c->curvaAnt;
//                 if( M->interface.curvaInicial==c )
//                     M->interface.curvaInicial = c->proxCurva;
//                 LiberaMemoriaCurva(c);
//                 return;
//             }

//         }
//     }

// }

// int RealizaEmbaracamentoForcado_Quebra(MALHA *M)
// {
//     int i, j, jEscolhido, i1, i2, qtdCelulas, precisaEmbaracar = 0, nenhumaSurface = 1;
//     PONTO *p, *p1Cima, *p2Cima, *p1Baixo, *p2Baixo;
//     CURVA *c;
//     double meioY, x, y;



//     ///Detectando se precisa quebrar
//     for( j=M->Ny-1; j>=0; j-- ) {
//         qtdCelulas = 0;
//         nenhumaSurface = 1;
//         for( i=0; i<M->Nx; i++ ) {
//             if( M->celulas[i][j]==SURFACE ) {
//                 qtdCelulas++;
//                 nenhumaSurface = 0;
//             }
//             else if( M->celulas[i][j]==FULL ) {
//                 qtdCelulas = 0;
//                 break;
//             }
//         }

//         if( nenhumaSurface ) {
//             precisaEmbaracar = 0;
//             break;
//         }

//         if( qtdCelulas == 2 ) {
//             precisaEmbaracar = 1;
//             break;
//         }

//     }

//     if( precisaEmbaracar==0 )
//         return 0;


//     jEscolhido = j;

// //    DesenhaMalhaVTK(*M, 0);
// //    ImprimeInterfaceVTK(*M, 1000000);
// //    PrintDebug("\n\n Vai Embaracar %d\n\n", jEscolhido);
// //    getchar();


//     i1 = i2 = -1;
//     for( i=0; i<M->Nx; i++ ) {
//         if( M->celulas[i][jEscolhido+1]==SURFACE && (i1==-1) )
//             i1 = i;
//         else if( M->celulas[i][jEscolhido+1]==SURFACE )
//             i2 = i;
//     }


//     //Encontrando os dois pontos que serao embaracados
//     p1Cima = p2Cima = NULL;
//     meioY = 0.5*( M->y[jEscolhido+1] + M->y[jEscolhido+2] );
//     double melhorY1 = 10.0, melhorY2 = 10.0, dist;
//     for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
//         for( p=c->nodeInicial; p; p=p->prox ) {
//             EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);
//             if( i==i1 ) {
//                 dist = fabs( meioY - p->y );
//                 if( dist<melhorY1 ) {
//                     melhorY1 = dist;
//                     p1Cima = p;
//                 }
//             }
//             if( i==i2 ) {
//                 dist = fabs( meioY - p->y );
//                 if( dist<melhorY2 ) {
//                     melhorY2 = dist;
//                     p2Cima = p;
//                 }
//             }

//         }
//     }

//     if( !p1Cima || !p2Cima )
//         return 0;


//     /// QUEBRANDO DE NOVO EMBAIXO
//     i1 = i2 = -1;
//     for( i=0; i<M->Nx; i++ ) {
//         if( M->celulas[i][jEscolhido-1]==SURFACE && (i1==-1) )
//             i1 = i;
//         else if( M->celulas[i][jEscolhido-1]==SURFACE )
//             i2 = i;
//     }


//     //Encontrando os dois pontos que serao embaracados
//     p1Baixo = p2Baixo = NULL;
//     meioY = 0.5*( M->y[jEscolhido-1] + M->y[jEscolhido] );
//     melhorY1 = 10.0;
//     melhorY2 = 10.0;
//     for( c=M->interface.curvaInicial; c!=NULL; c=c->proxCurva ) {
//         for( p=c->nodeInicial; p; p=p->prox ) {

//             EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);
//             if( i==i1 ) {
//                 dist = fabs( meioY - p->y );
//                 if( dist<melhorY1 ) {
//                     melhorY1 = dist;
//                     p1Baixo = p;
//                 }
//             }
//             if( i==i2 ) {
//                 dist = fabs( meioY - p->y );
//                 if( dist<melhorY2 ) {
//                     melhorY2 = dist;
//                     p2Baixo = p;
//                 }
//             }

//         }
//     }

//     if( !p1Baixo || !p2Baixo )
//         return 0;

//     //Trocando as posicoes pra embaracar
//     x = p1Cima->x;
//     y = p1Cima->y;
//     p1Cima->x = p2Cima->x;
//     p1Cima->y = p2Cima->y;
//     p2Cima->x = x;
//     p2Cima->y = y;


//     x = p1Baixo->x;
//     y = p1Baixo->y;
//     p1Baixo->x = p2Baixo->x;
//     p1Baixo->y = p2Baixo->y;
//     p2Baixo->x = x;
//     p2Baixo->y = y;


// //    DesenhaMalhaVTK(*M, 0);
// //    ImprimeInterfaceVTK(*M, 1000000);
// //    PrintDebug("\n\n Vai desembaracar %d\n\n", jEscolhido);
// //    getchar();

// //    ImprimeInterfaceVTK(*M, 100000);
// //    DeuProblema("AAAA %lf %lf %lf %lf\n", p1->x, p1->y, p2->x, p2->y);

//     return 1;
// }

void AtualizaMinimoMaximo(int i, int j, int *MinI, int *MinJ, int *MaxI, int *MaxJ)
{
    if( i<*MinI )
        *MinI = i;
    if( i>*MaxI )
        *MaxI = i;
    if( j<*MinJ )
        *MinJ = j;
    if( j>*MaxJ )
        *MaxJ = j;

    return;
}

// void AjustaNodesIniciais_CasoDuasListas(MALHA *M, LISTA_MUDANCA_TOPOL L)
// {
//     /// === Verificando se tem algum node na area de mudanca topologica
//     PONTO *pFinal = NULL;
//     CURVA *curvaEscolhida = NULL;
//     PONTO_MUDANCA_TOPOL *p1;
//     int encontrouInicio = 0;

//     for( p1=L.pontoInicial; p1!=NULL; p1=p1->prox ) {
//         CURVA *c;

//         for( c=M->interface.curvaInicial; c; c=c->proxCurva ) {
//             if( p1->ponto==c->nodeInicial ) {
//                 encontrouInicio = 1;
//                 curvaEscolhida = c;
//             }
//         }
//     }

//     if( encontrouInicio ) {
//         // Buscando o ponto final da lista
//         PONTO *p;
//         for( p=curvaEscolhida->nodeInicial; p; p=p->prox ) {
//             int encontrouPonto = 0;
//             for( p1=L.pontoInicial; p1!=NULL; p1=p1->prox ) {
//                 if( p==p1->ponto )
//                     encontrouPonto = 1;
//             }

//             if( !encontrouPonto ) {
//                 pFinal = p->ant;
//                 break;
//             }
//         }

//         if( !pFinal ) {
// //            ImprimeInterfaceVTK(*M, 100000);
// //            DesenhaMalhaVTK(*M, 0);
//             DeuProblema("AH ERA ESSE AQUI\n");
//         }
//         AndaComONodeDaCurva(curvaEscolhida, pFinal, 5);
//     }

//     return;
// }

// void EmbaracamentoCasoDuasListas(MALHA *M, LISTA_MUDANCA_TOPOL *L)
// {
//     PONTO_MUDANCA_TOPOL *p1, *p2, *pEscolhido;
//     double distancia, menorDistancia, vetX, vetY;
//     int i, j;

//     //Primeiramente calculando a direcao dos pontos na lista 0
//     for( p1=L[0].pontoInicial; p1!=NULL; p1=p1->prox ) {

//         //Aproveitando pra voltar essa celula pra FULL
//         EncontraCelulaDaParticula(M, p1->ponto->x, p1->ponto->y, &i, &j);
//         M->celulas[i][j] = FULL;

//         pEscolhido = NULL;
//         menorDistancia = 1e+10;
//         for( p2=L[1].pontoInicial; p2!=NULL; p2=p2->prox ) {
//             vetX = p2->ponto->x - p1->ponto->x;
//             vetY = p2->ponto->y - p1->ponto->y;
//             distancia = sqrt( vetX*vetX + vetY*vetY );
//             if( distancia<menorDistancia ) {
//                 pEscolhido = p2;
//                 menorDistancia = distancia;
//             }
//         }

//         if( pEscolhido!=NULL ) {
//             p1->direcaoMovX = pEscolhido->ponto->x - p1->ponto->x;
//             p1->direcaoMovY = pEscolhido->ponto->y - p1->ponto->y;
//         }
//     }

//     //agora a direcao dos pontos na lista 1
//     for( p1=L[1].pontoInicial; p1!=NULL; p1=p1->prox ) {

//         //Aproveitando pra voltar essa celula pra FULL
//         EncontraCelulaDaParticula(M, p1->ponto->x, p1->ponto->y, &i, &j);
//         M->celulas[i][j] = FULL;

//         pEscolhido = NULL;
//         menorDistancia = 1e+10;
//         for( p2=L[0].pontoInicial; p2!=NULL; p2=p2->prox ) {
//             vetX = p2->ponto->x - p1->ponto->x;
//             vetY = p2->ponto->y - p1->ponto->y;
//             distancia = sqrt( vetX*vetX + vetY*vetY );
//             if( distancia<menorDistancia ) {
//                 pEscolhido = p2;
//                 menorDistancia = distancia;
//             }
//         }

//         if( pEscolhido!=NULL ) {
//             p1->direcaoMovX = pEscolhido->ponto->x - p1->ponto->x;
//             p1->direcaoMovY = pEscolhido->ponto->y - p1->ponto->y;
//         }
//     }

//     //NAO MOVE NADA
//     if( L[0].pontoInicial==NULL || L[1].pontoInicial==NULL )
//         return;


//     /// === Verifica se tem nodeInicial na regiao que vai ser quebrada
//     AjustaNodesIniciais_CasoDuasListas(M, L[0]);
//     AjustaNodesIniciais_CasoDuasListas(M, L[1]);


//     //Movimentando os pontos
//     for( p1=L[0].pontoInicial; p1!=NULL; p1=p1->prox ) {
//         p1->ponto->x += p1->direcaoMovX;
//         p1->ponto->y += p1->direcaoMovY;
//     }
//     for( p1=L[1].pontoInicial; p1!=NULL; p1=p1->prox ) {
//         p1->ponto->x += p1->direcaoMovX;
//         p1->ponto->y += p1->direcaoMovY;
//     }

//     return;
// }

// void EmbaracamentoCasoUmaLista(MALHA *M, LISTA_MUDANCA_TOPOL Lista)
// {
//     //Vou apenas inverter a posicao dos pontos
//     //Tipo assim: o primeiro vira o ultimo. o segundo vira o penultimo. etc. isso vai fazer embaracar
//     PONTO *primPonto, *ultPonto;
//     PONTO_MUDANCA_TOPOL *p;
//     int qtdPontos, i, j;
//     double tempX, tempY;

//     primPonto = Lista.pontoInicial->ponto;
//     ultPonto = Lista.ultimoPonto->ponto;

//     //Contando quantos pontos tem na lista
//     qtdPontos = 0;
//     for( p=Lista.pontoInicial; p!=NULL; p=p->prox ) {
//         //Aproveitando pra voltar essa celula pra FULL
//         EncontraCelulaDaParticula(M, p->ponto->x, p->ponto->y, &i, &j);
//         M->celulas[i][j] = FULL;

//         qtdPontos++;
//     }

//     //Embaracando invertendo os pontos
//     for( ; qtdPontos>1; qtdPontos-=2 ) {
//         tempX = primPonto->x;
//         tempY = primPonto->y;

// //        ImprimeInterfaceVTK(*M, 10000000);
// //        printf("EMBARACOU %lf %lf %lf %lf\n", primPonto->x, primPonto->y, ultPonto->x, ultPonto->y);

//         primPonto->x = ultPonto->x;
//         primPonto->y = ultPonto->y;
//         ultPonto->x = tempX;
//         ultPonto->y = tempY;
//         primPonto = primPonto->prox;
//         ultPonto = ultPonto->ant;

// //        ImprimeInterfaceVTK(*M, 20000000);
// //        DeuProblema("parou\n");
//     }

//     return;
// }

// void RemoveCelulasIndesejadasMudancaTopologicaPasso1(MALHA *M, int *QtdCelulas)
// {
//     int i, j, qtdCelulasInicial;

//     qtdCelulasInicial = *QtdCelulas;

//     //Agora vou remover alguns casos que nao sao de mudanca topologica
//     for( i=0; i<M->Nx; i++ ) {
//         for( j=0; j<M->Ny; j++ ) {




//             if( M->celulas[i][j]!=MUDANCA_TOPOL )
//                 continue;


//             //Vou considerar que nao tem mudanca topologica nas celulas nos extremos do dominio (perto de parede)
//             //vou eliminar celulas que estao a 1 ou 2 celulas de distancia de uma boundary
//             //Depois eh bom verificar se isso nao pode acontecer as vezes...
//             int remove = 0;
//             if( (i<2) || (M->pontosU[i][j].tipo==BOUNDARY) || (M->pontosU[i-1][j].tipo==BOUNDARY) )
//                 remove = 1;
//             else if( (i>=M->Nx-2) || (M->pontosU[i+1][j].tipo==BOUNDARY) || (M->pontosU[i+2][j].tipo==BOUNDARY) )
//                 remove = 1;
//             else if( (j<2) || (M->pontosV[i][j].tipo==BOUNDARY) || (M->pontosV[i][j-1].tipo==BOUNDARY) )
//                 remove = 1;
//             else if( (j>=M->Ny-2) || (M->pontosV[i][j+1].tipo==BOUNDARY) || (M->pontosV[i][j+2].tipo==BOUNDARY) )
//                 remove = 1;
//             if( remove ) {
//                 M->celulas[i][j] = FULL;
//                 (*QtdCelulas)--;
//                 continue;
//             }

//             //Toda celula MUDANCA_TOPOL precisa ter uma outra MUDANCA_TOPOL do seu lado
//             if( M->celulas[i+1][j]!=MUDANCA_TOPOL && M->celulas[i-1][j]!=MUDANCA_TOPOL && M->celulas[i][j+1]!=MUDANCA_TOPOL && M->celulas[i][j-1]!=MUDANCA_TOPOL ) {
//                 M->celulas[i][j] = FULL;
//                 (*QtdCelulas)--;
//                 continue;
//             }


//             //Removendo tambem pedacos que possuem apenas duas celuas M. Nao quero fazer esses
//             if( M->celulas[i+1][j]==MUDANCA_TOPOL && M->celulas[i-1][j]!=MUDANCA_TOPOL && M->celulas[i][j+1]!=MUDANCA_TOPOL && M->celulas[i][j-1]!=MUDANCA_TOPOL ) {
//                 if( i!=M->Nx-2 && M->celulas[i+2][j]!=MUDANCA_TOPOL && M->celulas[i+1][j+1]!=MUDANCA_TOPOL && M->celulas[i+1][j-1]!=MUDANCA_TOPOL  ) {
//                     M->celulas[i][j] = FULL;
//                     M->celulas[i+1][j] = FULL;
//                     (*QtdCelulas) -= 2;
//                 }
//             }
//             if( M->celulas[i+1][j]!=MUDANCA_TOPOL && M->celulas[i-1][j]!=MUDANCA_TOPOL && M->celulas[i][j+1]==MUDANCA_TOPOL && M->celulas[i][j-1]!=MUDANCA_TOPOL ) {
//                 if( j!=M->Ny-2 && M->celulas[i][j+2]!=MUDANCA_TOPOL && M->celulas[i+1][j+1]!=MUDANCA_TOPOL && M->celulas[i-1][j+1]!=MUDANCA_TOPOL  ) {
//                     M->celulas[i][j] = FULL;
//                     M->celulas[i][j+1] = FULL;
//                     (*QtdCelulas) -= 2;
//                 }
//             }

//             //Vou tirar alguns casos de quina que nao sairam no teste acima
// //            if( M->celulas[i+1][j]==SURFACE && M->celulas[i-1][j]!=SURFACE && M->celulas[i][j+1]!=SURFACE && M->celulas[i][j-1]==SURFACE ) {
// //                M->celulas[i][j] = FULL;
// //                qtdCelulasMudanca--;
// //                continue;
// //            }
// //            else if( M->celulas[i+1][j]!=SURFACE && M->celulas[i-1][j]==SURFACE && M->celulas[i][j+1]!=SURFACE && M->celulas[i][j-1]==SURFACE ) {
// //                M->celulas[i][j] = FULL;
// //                qtdCelulasMudanca--;
// //                continue;
// //            }
// //            else if( M->celulas[i+1][j]==SURFACE && M->celulas[i-1][j]!=SURFACE && M->celulas[i][j+1]==SURFACE && M->celulas[i][j-1]!=SURFACE ) {
// //                M->celulas[i][j] = FULL;
// //                qtdCelulasMudanca--;
// //                continue;
// //            }
// //            else if( M->celulas[i+1][j]!=SURFACE && M->celulas[i-1][j]==SURFACE && M->celulas[i][j+1]==SURFACE && M->celulas[i][j-1]!=SURFACE ) {
// //                M->celulas[i][j] = FULL;
// //                qtdCelulasMudanca--;
// //                continue;
// //            }
//         }
//     }

//     if( (*QtdCelulas!=qtdCelulasInicial) && (*QtdCelulas>0) )
//         RemoveCelulasIndesejadasMudancaTopologicaPasso1(M, QtdCelulas);
// }

// void RemoveCelulasIndesejadasMudancaTopologicaPasso2(MALHA *M, int *QtdCelulas)
// {
//     int i, j, qtdCelulasInicial, qtdCorrecao, ultimoI = -1, ultimoJ = -1;
//     CURVA *curva;
//     PONTO *ponto;

//     qtdCelulasInicial = *QtdCelulas;

//     //Corrigindo alguns casos
//     for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//         qtdCorrecao = 0;
//         for( ponto=curva->nodeInicial; ponto!=NULL; ponto=ponto->prox ) {
//             EncontraCelulaDaParticula(M, ponto->x, ponto->y, &i, &j);

//             if( qtdCorrecao!=0 && qtdCorrecao<3 && M->celulas[i][j]!=MUDANCA_TOPOL  ) {
//                 //DeuProblema("ENTROU AQUI!!! %d %d\n\n", ultimoI, ultimoJ);
//                 M->celulas[ultimoI][ultimoJ] = FULL;
//                 qtdCorrecao = 0;
//                 (*QtdCelulas)--;
//             }

//             if( M->celulas[i][j]!=MUDANCA_TOPOL ) {
//                 qtdCorrecao = 0;
//                 continue;
//             }

//             if( qtdCorrecao==0 ) {
//                 qtdCorrecao++;
//                 ultimoI = i;
//                 ultimoJ = j;
//             }
//             else if( i!=ultimoI || j!=ultimoJ ) {
//                 qtdCorrecao++;
//                 ultimoI = i;
//                 ultimoJ = j;
//             }

//         }
//     }

//     if( (*QtdCelulas!=qtdCelulasInicial) && (*QtdCelulas>0) )
//         RemoveCelulasIndesejadasMudancaTopologicaPasso2(M, QtdCelulas);

//     return;
// }

// void AdicionaNaListaDeMudancaTopol(LISTA_MUDANCA_TOPOL *L, PONTO *Ponto)
// {
//     PONTO_MUDANCA_TOPOL *novoItem;

//     novoItem = (PONTO_MUDANCA_TOPOL *)malloc(sizeof(PONTO_MUDANCA_TOPOL));
//     novoItem->ponto = Ponto;
//     novoItem->ant = novoItem->prox = NULL;

//     if( L->pontoInicial==NULL ) {
//         L->pontoInicial = novoItem;
//         L->ultimoPonto = novoItem;
//     }
//     else {
//         L->ultimoPonto->prox = novoItem;
//         novoItem->ant = L->ultimoPonto;
//         L->ultimoPonto = novoItem;
//     }

//     return;
// }

// void RemoveParticulasFULL(MALHA *M)
// {
//     CURVA *curva;
//     PONTO *p, *pInicial, *pFinal;
//     int i, j, iAnt, jAnt, qtdCelulas = 0, qtdPontos;
//     int houveRemocao = 0;

// //    PrintDebug("Comecou remover particulas FULL\n");

//     /// Primeira passagem
//     for( curva=M->interface.curvaInicial; curva; curva=curva->proxCurva ) {

//         pInicial = NULL;
//         pFinal = NULL;
//         iAnt = -1;
//         jAnt = -1;

//         int numPonto = -1;
//         LOOP_CURVA(curva) {
//             p = loop_stuff.ponto;
//             numPonto++;

//             PONTO *p1, *p2, *p3;
//             p1 = p->ant;
//             p2 = p;
//             p3 = p->prox;
//             double vec1x = p2->x - p1->x;
//             double vec1y = p2->y - p1->y;
//             double vec2x = p3->x - p2->x;
//             double vec2y = p3->y - p2->y;
//             double produto = (vec1x)*(vec2x) + (vec1y)*(vec2y);
//             double norm1 = sqrt(vec1x*vec1x + vec1y*vec1y);
//             double norm2 = sqrt(vec2x*vec2x + vec2y*vec2y);
//             double angulo = acos(produto/(norm1*norm2));

//             /// Encontrei uma "quina" na interface. Eh um possivel local de pontos FULL pra remover
//             // if( angulo>M_PI/1.5 ) {
//             if( 1 ) {
//                 EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);

//                 // Soh estamos interessados em regioes FULL
//                 if( M->celulas[i][j]!=FULL )
//                     continue;
// //
//                 qtdPontos = 0;

//                 int iMin = i, iMax = i, jMin = j, jMax = j;
// //
//                 pInicial = p;
//                 while( M->celulas[i][j]==FULL ) {
//                     pInicial = pInicial->ant;
//                     EncontraCelulaDaParticula(M, pInicial->x, pInicial->y, &i, &j);
//                     qtdPontos++;

//                     // Nao quero remover pontos do eixo de simetria no caso axissimetrico
//                     if( (M->tipoCoord == AXI_CILINDRICO) && (pInicial->x<1e-6) )
//                         break;

//                     iMin = (i<iMin) ? i : iMin;
//                     jMin = (j<jMin) ? j : jMin;
//                     iMax = (i>iMax) ? i : iMax;
//                     jMax = (j>jMax) ? j : jMax;
//                 }

//                 pFinal = p;
//                 EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);
//                 while( M->celulas[i][j]==FULL ) {
//                     pFinal = pFinal->prox;
//                     EncontraCelulaDaParticula(M, pFinal->x, pFinal->y, &i, &j);
//                     qtdPontos++;

//                     // Nao quero remover pontos do eixo de simetria no caso axissimetrico
//                     if( (M->tipoCoord == AXI_CILINDRICO) && (pFinal->x<1e-6) )
//                         break;

//                     iMin = (i<iMin) ? i : iMin;
//                     jMin = (j<jMin) ? j : jMin;
//                     iMax = (i>iMax) ? i : iMax;
//                     jMax = (j>jMax) ? j : jMax;
//                 }

//                 /// Encontrou uma regiao com mais de 50 pontos dentro de celulas FULL
//                 /// Vou apagar todos estes pontos!!! (perigoso?)
//                 if( (qtdPontos>50) && (iMax-iMin<5) ) {
//                     PrintDebug("removendo particulas full...\n");
//                     PrintDebug("INICIAL FINAL: %.15lf %.15lf %.15lf %.15lf\n", pInicial->x, pInicial->y, pFinal->x, pFinal->y);
//                     PrintDebug("ponto mudanca: %d %d: %.15lf %.15lf\n", numPonto, qtdPontos, p->x, p->y);
//                     PrintDebug("iMax-iMin: %d ... jMax-jMin: %d\n", iMax-iMin, jMax-jMin);
//                     // ImprimeInterfaceVTK(*M, 400000);

//                     PONTO *proxPonto = pInicial->prox;
//                     pInicial->prox = pFinal;
//                     pFinal->ant = pInicial;

//                     // Precisa tratar este edge case depois... (fiquei com preguica)
// //                    if( pInicial->isNode || pFinal->isNode )
// //                        DeuProblema("1. RemoveParticulasFull: tentou apagar um node..\n");

//                     int apagouNode = 0;
//                     for( p=proxPonto; p!=pFinal; p=proxPonto ) {
//                         // Precisa tratar este edge case depois... (fiquei com preguica)
//                         if( p->isNode ) {
//                             apagouNode = 1;
// //                            DeuProblema("2. RemoveParticulasFull: tentou apagar um node..\n");
//                         }

//                         proxPonto = p->prox;
//                         free(p);
//                     }

//                     if( apagouNode ) {
//                         PrintDebug("Trocou node\n");
//                         curva->pontoInicial = pInicial;
//                         pInicial->isNode = 1;

//                         PONTO *nodeFinal = (PONTO *)malloc(sizeof(PONTO));
//                         nodeFinal->prox = NULL;
//                         nodeFinal->ant = pInicial->ant;
//                         nodeFinal->x = pInicial->x;
//                         nodeFinal->y = pInicial->y;
//                         nodeFinal->fixo = 0;
//                         nodeFinal->outflow = 0;
//                         nodeFinal->isNode = 1;

//                         pInicial->ant->prox = nodeFinal;
//                         pInicial->ant = NULL;
//                     }

//                     // DesenhaMalhaVTK(*M, 0);
//                     // ImprimeInterfaceVTK(*M, 500000);
//                     // DeuProblema("PARTOU\n");

//                     break;
//                 }

//             }
//         }


//     }

// //    PrintDebug("Metade do remover particulas FULL\n");

//     /// Segunda passagem
//     for( curva=M->interface.curvaInicial; curva; curva=curva->proxCurva ) {

//         pInicial = NULL;
//         pFinal = NULL;
//         iAnt = -1;
//         jAnt = -1;

//         for( p=curva->nodeInicial; p; p=p->prox ) {

//             EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);

//             /// Casos para desconsiderar - paredes
//             int paredeCima, paredeBaixo, paredeEsq, paredeDir;
//             paredeEsq = (M->pontosU[i][j].tipo==BOUNDARY) || (i!=0 && M->pontosU[i-1][j].tipo==BOUNDARY); //Estou vendo duas celulas pra esquerda devido a simulacoes de tubo que a interface nao fica na celula da parede
//             paredeDir = (M->pontosU[i+1][j].tipo==BOUNDARY) || ((i!=M->Nx-1) && M->pontosU[i+2][j].tipo==BOUNDARY);
//             paredeBaixo = (M->pontosV[i][j].tipo==BOUNDARY);
//             paredeCima = (M->pontosV[i][j+1].tipo==BOUNDARY);
//             if(  paredeEsq || paredeDir || paredeBaixo || paredeCima  )
//                 pInicial = pFinal = NULL;

//             if( (pInicial==NULL) && (M->celulas[i][j]==FULL) ) {
//                 qtdCelulas = 0;
//                 iAnt = jAnt = -1;
//                 pInicial = p;
//                 pFinal = NULL;
//             }
//             else if( (pInicial!=NULL) && (M->celulas[i][j]==SURFACE) )
//                 pFinal = p;

//             if( (pInicial!=NULL) && (i!=iAnt || j!=jAnt) )
//                 qtdCelulas++;

//             if( pInicial && pFinal ) {
//                 if( qtdCelulas < 4 )
//                     pInicial = pFinal = NULL;
//                 else {
//                     //Remove os pontos entre pInicial e pFinal (incluindo eles)

//                     //Nao quero remover exatamente o primeiro ponto da interface e nem o ultimo
//                     if( pInicial->ant==NULL )
//                         pInicial = pInicial->prox;
//                     if( pFinal->prox==NULL )
//                         pFinal = pFinal->ant;

//                     pInicial->ant->prox = pFinal;
//                     pFinal->ant = pInicial->ant;

//                     PONTO *pRemove, *pProx = NULL;
//                     for( pRemove=pInicial; pRemove!=pFinal; pRemove=pProx ) {
//                         pProx = pRemove->prox;
//                         free(pRemove);
//                     }

//                     houveRemocao = 1;
//                     pInicial = pFinal = NULL;

// //                    DeuProblema("ERMOVE FULLa aaaa\n\n");
//                 }
//             }

//             iAnt = i;
//             jAnt = j;

//         }


//     }

// //    PrintDebug("saiu 2\n");

//     if( houveRemocao ) {
//         AdicionaERemoveParticulas(M, 0.05*M->minDxDy, 0);
// //        ImprimeInterfaceVTK(*M, 200000);
// //        DesenhaMalhaVTK(*M, 0);
// //        DeuProblema("ENCONTROU remocao!\n\n");
//     }

// //    ImprimeInterfaceVTK(*M, 500000);
// //    PrintDebug("Terminou remover particulas FULL\n");
//     return;
// }

// void RemoveParticulasNaSimetriaSURFACE(MALHA *M)
// {
//     CURVA *curva;
//     PONTO *p, *pInicial, *pFinal;
//     int i, j, qtdPontos1, qtdPontos2;
//     int iMin, iMax;
// //    int houveRemocao = 0;

// //    PrintDebug("Comecou remover particulas SURFACE\n");

//     /// Primeira passagem
//     for( curva=M->interface.curvaInicial; curva; curva=curva->proxCurva ) {

//         pInicial = NULL;
//         pFinal = NULL;

//         int numPonto = -1;
//         for( p=curva->nodeInicial; p; p=p->prox ) {
//             numPonto++;

//             if( p->fixo )
//                 continue;

//             PONTO *p1, *p2, *p3;
//             p1 = (p->ant) ? p->ant : p->pontoAntCurva ;
//             p2 = p;
//             p3 = (p->prox) ? p->prox : p->pontoProxCurva;
//             double vec1x = p2->x - p1->x;
//             double vec1y = p2->y - p1->y;
//             double vec2x = p3->x - p2->x;
//             double vec2y = p3->y - p2->y;
//             double produto = (vec1x)*(vec2x) + (vec1y)*(vec2y);
//             double norm1 = sqrt(vec1x*vec1x + vec1y*vec1y);
//             double norm2 = sqrt(vec2x*vec2x + vec2y*vec2y);
//             double angulo = acos(produto/(norm1*norm2));

//             /// Encontrei uma "quina" na interface. Eh um possivel local de pontos SURFACE pra remover
//             if( angulo>M_PI/1.5 ) {
//                 EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);

//                 iMin = i;
//                 iMax = i;

//                 // Soh estamos interessados em regioes SURFACE
//                 if( M->celulas[i][j]!=SURFACE )
//                     continue;
// //
//                 qtdPontos1 = qtdPontos2 = 0;
// //
//                 pInicial = p;
//                 while( M->celulas[i][j]==SURFACE ) {
//                     pInicial = (pInicial->ant) ? pInicial->ant : pInicial->pontoAntCurva;
//                     EncontraCelulaDaParticula(M, pInicial->x, pInicial->y, &i, &j);
//                     qtdPontos1++;

//                     iMin = (i<iMin) ? i : iMin;
//                     iMax = (i>iMax) ? i : iMax;
//                 }

//                 pFinal = p;
//                 EncontraCelulaDaParticula(M, p->x, p->y, &i, &j);
//                 while( M->celulas[i][j]==SURFACE ) {
//                     pFinal = (pFinal->prox) ? pFinal->prox : pFinal->pontoProxCurva;
//                     EncontraCelulaDaParticula(M, pFinal->x, pFinal->y, &i, &j);
//                     qtdPontos2++;

//                     iMin = (i<iMin) ? i : iMin;
//                     iMax = (i>iMax) ? i : iMax;
//                 }

// //                PrintDebug("INICIAL FINAL: %.15lf %.15lf %.15lf %.15lf\n", pInicial->x, pInicial->y, pFinal->x, pFinal->y);
// //                PrintDebug("ponto mudanca: %d %d %d: %.15lf %.15lf\n", numPonto, qtdPontos1, qtdPontos2, p->x, p->y);
// //                ImprimeInterfaceVTK(*M, 400000);
// //                getchar();

//                 /// Encontrou uma regiao com mais de 50 pontos dentro de celulas FULL
//                 /// Vou apagar todos estes pontos!!! (perigoso?)
//                 if( (qtdPontos1+qtdPontos2>40) && (iMax-iMin)<2 ) {
//                     PrintDebug("caso 1 qtdPontos: %d %d %d\n", qtdPontos1, qtdPontos2, iMax-iMin);

//                 //    ImprimeInterfaceVTK(*M, 300000);
//                 //    DesenhaMalhaVTK(*M, 0);
// //                    PrintDebug("MIN MAX: %d %d\n", iMin, iMax);

// //                    PrintDebug("INICIAL FINAL: %.15lf %.15lf %.15lf %.15lf\n", pInicial->x, pInicial->y, pFinal->x, pFinal->y);
// //                    PrintDebug("ponto mudanca: %d %d %d: %.15lf %.15lf\n", numPonto, qtdPontos1, qtdPontos2, p->x, p->y);
// //                    ImprimeInterfaceVTK(*M, 400000);

//                     PONTO *proxPonto = (pInicial->prox) ? pInicial->prox : pInicial->pontoProxCurva;
//                     pInicial->prox = pFinal;
//                     pFinal->ant = pInicial;

//                     // Precisa tratar este edge case depois... (fiquei com preguica)
//                     if( pInicial->isNode || pFinal->isNode )
//                         DeuProblema("1. RemoveParticulasFull: tentou apagar um node..\n");

//                     int apagouNode = 0;
//                     for( p=proxPonto; p!=pFinal; p=proxPonto ) {
//                         // Precisa tratar este edge case depois... (fiquei com preguica)
//                         if( p->isNode ) {
//                             apagouNode = 1;
// //                            DeuProblema("2. RemoveParticulasFull: tentou apagar um node..\n");
//                         }

//                         proxPonto = (p->prox) ? p->prox : p->pontoProxCurva;
//                         free(p);
//                     }

//                     if( apagouNode ) {
//                         PrintDebug("Trocou node\n");
//                         curva->nodeInicial = pInicial;
//                         pInicial->isNode = 1;

//                         PONTO *nodeFinal = (PONTO *)malloc(sizeof(PONTO));
//                         nodeFinal->prox = NULL;
//                         nodeFinal->ant = pInicial->ant;
//                         nodeFinal->pontoProxCurva = pInicial->prox;
//                         nodeFinal->x = pInicial->x;
//                         nodeFinal->y = pInicial->y;
//                         nodeFinal->fixo = 0;
//                         nodeFinal->outflow = 0;
//                         nodeFinal->isNode = 1;

//                         pInicial->ant->prox = nodeFinal;
//                         pInicial->pontoAntCurva = pInicial->ant;
//                         pInicial->ant = NULL;
//                     }

//                     // ImprimeInterfaceVTK(*M, 500000);
//                     // DeuProblema("parou\n");

//                     break;
//                 }

//             }
//         }


//     }



//     /// === Segunda passagem: olhando se tem uma coluna de celulas surface
//     int comecouSurface = 0, terminouSurface = 0;
//     int qtdSurfaces = 0;
//     for( j=0; j<M->Ny; j++ ) {
        
        
//         if( !comecouSurface && (M->celulas[0][j]==SURFACE) && (M->celulas[1][j]==EMPTY) ) {
//             comecouSurface = 1;
//             terminouSurface = 0;
//             qtdSurfaces = 1;
//         }
//         else if( comecouSurface && (M->celulas[0][j]==SURFACE) && (M->celulas[1][j]==EMPTY) )
//             qtdSurfaces++;
//         else if( comecouSurface ) {
//             terminouSurface = 1;
//             comecouSurface = 0;
//         }
        
//         /// Vai remover
//         if( terminouSurface && qtdSurfaces>=3 ) {
//             qtdSurfaces = 0;
//             terminouSurface = 0;

//             /// Encontra um ponto desta celula 
//             for( curva=M->interface.curvaInicial; curva; curva=curva->proxCurva ) {
//                 for( p=curva->nodeInicial; p; p=p->prox ) {
//                     int iMarker, jMarker;
//                     EncontraCelulaDaParticula(M, p->x, p->y, &iMarker, &jMarker);

//                     if( iMarker!=0 || jMarker!=j-1 )
//                         continue;

//                     qtdPontos1 = qtdPontos2 = 0;
//                     pInicial = p;
//                     while( M->celulas[iMarker][jMarker]==SURFACE ) {
//                         pInicial = (pInicial->ant) ? pInicial->ant : pInicial->pontoAntCurva;
//                         EncontraCelulaDaParticula(M, pInicial->x, pInicial->y, &iMarker, &jMarker);
//                         qtdPontos1++;
//                     }

//                     pFinal = p;
//                     EncontraCelulaDaParticula(M, p->x, p->y, &iMarker, &jMarker);
//                     while( M->celulas[iMarker][jMarker]==SURFACE ) {
//                         pFinal = (pFinal->prox) ? pFinal->prox : pFinal->pontoProxCurva;
//                         EncontraCelulaDaParticula(M, pFinal->x, pFinal->y, &iMarker, &jMarker);
//                         qtdPontos2++;
//                     }


//                     /// Encontrou uma regiao com mais de 50 pontos dentro de celulas SURFACE simetria
//                     /// Vou apagar todos estes pontos!!! (perigoso?)
//                     if( (qtdPontos1+qtdPontos2>30) ) {
//                         // ImprimeInterfaceVTK(*M, 500000);

//                         PrintDebug("CASO 2: removeu simetria...\n");

//                         PONTO *proxPonto = (pInicial->prox) ? pInicial->prox : pInicial->pontoProxCurva;
//                         pInicial->prox = pFinal;
//                         pFinal->ant = pInicial;


//                         // Precisa tratar este edge case depois... (fiquei com preguica)
//                         if( pInicial->isNode || pFinal->isNode )
//                             DeuProblema("1. RemoveParticulasFull: tentou apagar um node..\n");

//                         int apagouNode = 0;
//                         for( p=proxPonto; p!=pFinal; p=proxPonto ) {
//                             // Precisa tratar este edge case depois... (fiquei com preguica)
//                             if( p->isNode ) {
//                                 apagouNode = 1;
//     //                            DeuProblema("2. RemoveParticulasFull: tentou apagar um node..\n");
//                             }

//                             proxPonto = (p->prox) ? p->prox : p->pontoProxCurva;
//                             free(p);
//                         }

//                         if( apagouNode ) {
//                             PrintDebug("Trocou node\n");
//                             curva->nodeInicial = pInicial;
//                             pInicial->isNode = 1;

//                             PONTO *nodeFinal = (PONTO *)malloc(sizeof(PONTO));
//                             nodeFinal->prox = NULL;
//                             nodeFinal->ant = pInicial->ant;
//                             nodeFinal->pontoProxCurva = pInicial->prox;
//                             nodeFinal->x = pInicial->x;
//                             nodeFinal->y = pInicial->y;
//                             nodeFinal->fixo = 0;
//                             nodeFinal->outflow = 0;
//                             nodeFinal->isNode = 1;

//                             pInicial->ant->prox = nodeFinal;
//                             pInicial->pontoAntCurva = pInicial->ant;
//                             pInicial->ant = NULL;
//                         }

//                         // DesenhaMalhaVTK(*M, 0);
//                         // ImprimeInterfaceVTK(*M, 600000);
//                         // DeuProblema("PAROU\n");
//                         break;
//                     }

//                 }
//             } 
            
//         }
//         else if( terminouSurface ) {
//             qtdSurfaces = 0;
//             terminouSurface = 0;
//             comecouSurface = 0;
//         }


//     }

// }

// void RemoveCurvasMuitoPequenas(MALHA *M, int MinimoParticulas)
// {
//     CURVA *c, *proxCurva;
//     PONTO *p;//, *pProx;

//     c = M->interface.curvaInicial;
//     while( c ) {
//         proxCurva = c->proxCurva;

//         int qtdParticulas = 0;
//         for( p=c->nodeInicial; p; p=p->prox ) {
//             qtdParticulas++;

//             PONTO *proxPonto = (p->prox) ? p->prox : p->pontoProxCurva;
//             if( (fabs(proxPonto->x - p->x) + fabs(proxPonto->y - p->y))<1e-10 ) {
//                 qtdParticulas = -1; // APenas pra remover a curva
//                 break;
//             }

// //            if( qtdParticulas>=MinimoParticulas )
// //                break;

//         }

//         if( qtdParticulas<MinimoParticulas ) {
//             if( c->curvaAnt )
//                 c->curvaAnt->proxCurva = proxCurva;
//             if( proxCurva )
//                 proxCurva->curvaAnt = c->curvaAnt;
//             if( M->interface.curvaInicial==c )
//                 M->interface.curvaInicial = proxCurva;
//             LiberaMemoriaCurva(c);

//             PrintDebug("RETIROU UMA CURVA MUITO PEQUENA OU ESQUISITA\n");
//         }

//         c = proxCurva;
//     }

//     return;
// }

/// === Realiza um checkup geral na estrutura da interface (pra ver se nao fiz nenhuma coisa errada com os ponteiros...)
// void CheckUpInterface(MALHA *M, const char *NomeTeste)
// {
//     PONTO *p, *pInicial = NULL, *pFinal = NULL;
//     CURVA *c;

// //    PrintDebug("CHECKUP: %s\n", NomeTeste);

//     for( c=M->interface.curvaInicial; c; c=c->proxCurva ) {
//         pInicial = c->pontoInicial;

        
//         for( p=pInicial->prox; p; p=p->prox ) {
//             if( p->ant->prox!=p )
//                 DeuProblema("CheckUpInterface %s: Erro Numero 2\n\n", NomeTeste);

//             if( p->prox==NULL )
//                 pFinal = p;
//             else {
//                 if( p->isNode )
//                     p->isNode = 0;
// //                    DeuProblema("CheckUpInterface %s: Erro Numero 1\n\n", NomeTeste);

//                 if( p->prox->ant!=p )
//                     DeuProblema("CheckUpInterface %s: Erro Numero 3 ... %lf %lf\n\n", NomeTeste, p->x, p->y);
//             }

//             if( p==pInicial )
//                 DeuProblema("CheckUpInterface %s: Erro Numero 8\n\n", NomeTeste);
//         }

//         if( !pInicial->isNode ) {
//             DeuProblema("CheckUpInterface %s: Erro Numero 4  %p\n\n", NomeTeste, pInicial);
//         }
//         if( !pFinal->isNode )
//             DeuProblema("CheckUpInterface %s: Erro Numero 5  %p\n\n", NomeTeste, pInicial);

//         if( pInicial->ant!=NULL )
//             DeuProblema("CheckUpInterface %s: Erro Numero 6\n\n", NomeTeste);

//         if( pInicial->pontoAntCurva!=pFinal->ant )
//             pInicial->pontoAntCurva = pFinal->ant;
//             // DeuProblema("CheckUpInterface %s: Erro Numero 7\n\n", NomeTeste);

//         if( pFinal->pontoProxCurva!=pInicial->prox )
//             pFinal->pontoProxCurva = pInicial->prox;
//             // DeuProblema("CheckUpInterface %s: Erro Numero 7\n\n", NomeTeste);
//     }

// }

// void SuavizaOndulacoes(MALHA *M)
// {
//     int i, iPonto, max, qtdPontos, qtdCurvas, maxPontos, totalPontos;
//     CURVA *curva;
//     PONTO *p, *pontoInicial = NULL;
//     double *beta, *x, *y, deltaT;
//     double t;

//     //DeuProblema("NAO QUERO SUAVIZAR\n");

//     //#define CONDICAO_SUAVIZAR ( p->y<9.85 )
//     //#define CONDICAO_SUAVIZAR ( 1 )

//     PrintDebug("SUAVIZANDO!\n");

//     maxPontos = 5000;
//     beta = (double *)malloc( (maxPontos + 100)*sizeof(double) );
//     x = (double *)malloc( (maxPontos + 100)*sizeof(double) );
//     y = (double *)malloc( (maxPontos + 100)*sizeof(double) );
//     for( curva = M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {

//         pontoInicial = curva->nodeInicial;

//         PrintDebug("PONTO INICIAL: %lf %lf\n", pontoInicial->x, pontoInicial->y);

//         /// Contando quantos pontos serao usados
//         qtdPontos = 0;
//         for( p=pontoInicial; p!=NULL; p=p->prox ) {
// //            if( !(p->y<21.85) )  //Nao vou usar os fixos nesse caso
// //                continue;

//             qtdPontos++;
//         }

//         totalPontos = qtdPontos;

//         /// Contando quantas curvas vamos usar
//         qtdCurvas = 1;
//         while( qtdPontos>maxPontos ) {
//             qtdPontos /= 2;
//             qtdCurvas *= 2;
//         }


//         for( ; qtdCurvas; qtdCurvas-- ) {

//             i = 0;
//             for( p=pontoInicial; p!=NULL; p=p->prox) {
// //                if( !(p->y<21.85) ) { //Nao vou usar os fixos nesse caso
// //                    continue;
// //                }

//                 x[i] = p->x;
//                 y[i] = p->y;

//                 i++;
//                 totalPontos--;
//                 if( (i==qtdPontos) && (qtdCurvas==1) && totalPontos ) { //Completa com os pontos que faltaram na ultima curva
//                     qtdPontos++;
//                 }
//                 else if( i==qtdPontos )
//                     break;
//             }

//             deltaT = 1.0/(qtdPontos-1);


//             t = 0.0;

//             iPonto = 0;
//             for( p=pontoInicial; p!=NULL; p=p->prox) {

//                 //if( !(p->y<20.0) ) {
//                 if(  (p->fixo) ) {
//                     t = t + deltaT;
//                     iPonto++;
//                     continue;
//                 }

//                 //Inicializando beta com os valores de x
//                 for( i=0; i<qtdPontos; i++ )
//                     beta[i] = x[i];

//                 for( max=qtdPontos-1; max>0; max-- ) {
//                     for( i=0; i<max; i++ ) {
//                         beta[i] = beta[i]*(1.0-t) + beta[i+1]*t;
//                     }
//                 }

//                 p->x = beta[0];

//                 //Inicializando beta com os valores de y
//                 for( i=0; i<qtdPontos; i++ )
//                     beta[i] = y[i];

//                 for( max=qtdPontos-1; max>0; max-- ) {
//                     for( i=0; i<max; i++ ) {
//                         beta[i] = beta[i]*(1.0-t) + beta[i+1]*t;
//                     }
//                 }

//                 p->y = beta[0];

//                 if( p->y >= M->yMax )
//                     p->y = M->yMax - 1e-8;


//                 t = t + deltaT;

//                 iPonto++;
//                 if( (iPonto==qtdPontos) && (qtdCurvas!=1) )
//                     break;
//             }


//             pontoInicial = NULL;
//             if( p!=NULL )
//                 pontoInicial = p->prox;
//         }
//     }


//     free(beta);
//     free(x);
//     free(y);
//     return;
// }

// void SuavizaOndulacoesAoRedorDoPonto(MALHA *M, double X, double Y, int qtdPontos)
// {
//     int i, iPonto, max;
//     CURVA *curvaTemp;
//     PONTO *p, *pontoInicial = NULL, *pivo = NULL;
//     double *beta, *x, *y, deltaT;
//     double t;

//     //DeuProblema("NAO QUERO SUAVIZAR\n");

//     //#define CONDICAO_SUAVIZAR ( p->y<9.85 )
//     //#define CONDICAO_SUAVIZAR ( 1 )

//     PrintDebug("SUAVIZANDO!\n");

//     beta = (double *)malloc( (qtdPontos)*sizeof(double) );
//     x = (double *)malloc( (qtdPontos)*sizeof(double) );
//     y = (double *)malloc( (qtdPontos)*sizeof(double) );
//     //for( curva = M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {

//         double minDist = 1.0e+10, distancia;
//         for( curvaTemp = M->interface.curvaInicial; curvaTemp!=NULL; curvaTemp=curvaTemp->proxCurva ) {

//             for( p=curvaTemp->nodeInicial; p!=NULL; p=p->prox ) {
//                 distancia = (p->x - X)*(p->x - X) + (p->y - Y)*(p->y - Y);
//                 if( distancia < minDist ) {
//                     minDist = distancia;
//                     pivo = p;
//                 }
//             }
//         }



// //        if( menorDistancia>M->minDxDy )
// //            continue;


//         pontoInicial = pivo;
//         for( i=0; i<qtdPontos/2; i++ ) {
//             if( pontoInicial->ant == NULL )
//                 break;
//             pontoInicial = pontoInicial->ant;
//         }

//         PrintDebug("PIVO: %lf %lf\nP INICIAL: %lf %lf\n", pivo->x, pivo->y, pontoInicial->x, pontoInicial->y);



//         i = 0;
//         for( p=pontoInicial; p!=NULL; p=p->prox) {

//             x[i] = p->x;
//             y[i] = p->y;

//             i++;
//             if( i==qtdPontos )
//                 break;
//         }
//         qtdPontos = i;

//         deltaT = 1.0/(qtdPontos-1);


//         t = 0.0;

//         iPonto = 0;
//         for( p=pontoInicial; p!=NULL; p=p->prox) {

// //            if(  !(p->y<21.85) || (p->fixo) ) {
// //                t = t + deltaT;
// //                iPonto++;
// //                continue;
// //            }

//             //Inicializando beta com os valores de x
//             for( i=0; i<qtdPontos; i++ )
//                 beta[i] = x[i];

//             for( max=qtdPontos-1; max>0; max-- ) {
//                 for( i=0; i<max; i++ ) {
//                     beta[i] = beta[i]*(1.0-t) + beta[i+1]*t;
//                 }
//             }

//             p->x = beta[0];

//             //Inicializando beta com os valores de y
//             for( i=0; i<qtdPontos; i++ )
//                 beta[i] = y[i];

//             for( max=qtdPontos-1; max>0; max-- ) {
//                 for( i=0; i<max; i++ ) {
//                     beta[i] = beta[i]*(1.0-t) + beta[i+1]*t;
//                 }
//             }

//             p->y = beta[0];

//             if( p->y >= M->yMax )
//                 p->y = M->yMax - 1e-8;


//             t = t + deltaT;

//             iPonto++;
//             if( iPonto==qtdPontos )
//                 break;
//         }

//     //}


//     free(beta);
//     free(x);
//     free(y);
//     return;
// }


///FUNCOES USADAS PRA CRIAR A FUNCAO QUE DA O VALOR PRA UMA CONDICAO DE CONTORNO NOSLIP EM MOVIMENTO
double FuncaoValorNoSlipHorizontal(MALHA M, int i, int j, double X, double T, double **U)
{
    double termo;
    int iV = (i==M.Nx) ? (i-1) : i;

    /// === Vendo o tipo da condicao de contorno nesta parede horizontal
    /// === Precisa olhar dois pontos, pois pode estar bem em uma quina (tipo no caso da Contracao)
    int paredeBaixo = NAO_BOUNDARY, paredeCima = NAO_BOUNDARY;

    // Vendo primeiramente se tem contorno em uma parede de cima
    int tipoDir = M.pontosV[iV][j+1].tipoBoundary;
    int tipoEsq = (iV==0) ? NAO_BOUNDARY : M.pontosV[iV-1][j+1].tipoBoundary;
    if( tipoDir==tipoEsq )
        paredeCima = tipoEsq;
    else if( tipoEsq==NAO_BOUNDARY && tipoDir!=NAO_BOUNDARY )
        paredeCima = tipoDir;
    else if( tipoEsq!=NAO_BOUNDARY && tipoDir==NAO_BOUNDARY )
        paredeCima = tipoEsq;

    // Vendo agora se tem contorno em uma parede de baixo
    tipoDir = M.pontosV[iV][j].tipoBoundary;
    tipoEsq = (iV==0) ? NAO_BOUNDARY :M.pontosV[iV-1][j].tipoBoundary;
    if( tipoDir==tipoEsq )
        paredeBaixo = tipoEsq;
    else if( tipoEsq==NAO_BOUNDARY && tipoDir!=NAO_BOUNDARY )
        paredeBaixo = tipoDir;
    else if( tipoEsq!=NAO_BOUNDARY && tipoDir==NAO_BOUNDARY )
        paredeBaixo = tipoEsq;


    // Se nem a de cima nem a de baixo forem parede, eh pq deu alguma coisa errada nos meus testes acima
    // Ai tem que debugar o que aconteceu
    if( (!paredeBaixo && !paredeCima) || (paredeBaixo==paredeCima) )
        DeuProblema("Problema ValorNoSlipHorizontal: nenhuma parede....");


    /// === Condicao NOSLIP em CIMA
    if( paredeCima ) {
        if( M.pontosV[iV][j+1].movimento==MOVIMENTO_X ) // Se tiver movimento
            // return 8*(1 + tanh(8*T-4))*X*X*(1-X)*(1-X);
            return 1.0;
        else // Se nao tiver, retorna u=0
            return 0.0;
    }

    if( paredeBaixo ) {
        if( M.pontosV[iV][j].movimento==MOVIMENTO_X ) // Se tiver movimento
            // return 8*(1 + tanh(8*T-4))*X*X*(1-X)*(1-X);
            return -1.0;
        else // Se nao tiver, retorna u=0
            return 0.0;
    }

//    if( paredeBaixo==NOSLIP || paredeCima==NOSLIP ) {
////        return M.pontosV[iV][j+1].valorDirichlet;
//        return 8*(1 + tanh(8*T-4))*X*X*(1-X)*(1-X);
//    }

    /// === Condicao SLIP EM UMA PAREDE HORIZONTAL
    //Se estiver em uma parede superior
    if( paredeCima==SLIP ) {
        double k = M.pontosV[iV][j+1].valorDirichlet;
        termo = ( 2.0*k*M.beta )/( M.Re*M.dy[j] );
        return ( termo*U[i][j] )/( 1.0 + termo );
    }

    //Se estiver em uma parede inferior
    if( paredeBaixo==SLIP ) {
        double k = M.pontosV[iV][j].valorDirichlet;
        termo = ( 2.0*k*M.beta )/( M.Re*M.dy[j] );
        return ( termo*U[i][j] )/( 1.0 + termo );
    }


    DeuProblema("VALOR NO SLIP HORIZONTAL: PROBLEMA %d %d %d %d\n", i, j, paredeBaixo, paredeCima);
    return 0.0;

    //return 0.0;
    //return 8*(1 + tanh(8*T-4))*X*X*(1-X)*(1-X);
}

double FuncaoValorNoSlipVertical(MALHA M, int i, int j, double Y, double T, double **V)
{
    double termo;
    int jU = (j==M.Ny) ? (j-1) : j;

    /// === Vendo o tipo da condicao de contorno nesta parede horizontal
    /// === Precisa olhar dois pontos, pois pode estar bem em uma quina (tipo no caso da Contracao)
    int paredeEsq = NAO_BOUNDARY, paredeDir = NAO_BOUNDARY;

    // Vendo primeiramente se tem contorno em uma parede da direita
    int tipoCima = M.pontosU[i+1][jU].tipoBoundary;
    int tipoBaixo = (jU==0) ? NAO_BOUNDARY : M.pontosU[i+1][jU-1].tipoBoundary;
    if( tipoCima==tipoBaixo )
        paredeDir = tipoCima;
    else if( tipoCima==NAO_BOUNDARY && tipoBaixo!=NAO_BOUNDARY )
        paredeDir = tipoBaixo;
    else if( tipoCima!=NAO_BOUNDARY && tipoBaixo==NAO_BOUNDARY )
        paredeDir = tipoCima;

    // Vendo agora se tem contorno em uma parede da esquerda
    tipoCima = M.pontosU[i][jU].tipoBoundary;
    tipoBaixo = (jU==0) ? NAO_BOUNDARY : M.pontosU[i][jU-1].tipoBoundary;
    if( tipoCima==tipoBaixo )
        paredeEsq = tipoCima;
    else if( tipoCima==NAO_BOUNDARY && tipoBaixo!=NAO_BOUNDARY )
        paredeEsq = tipoBaixo;
    else if( tipoCima!=NAO_BOUNDARY && tipoBaixo==NAO_BOUNDARY )
        paredeEsq = tipoCima;

//    PrintDebug("TesteV: %d %d %d %d\n", i, j, paredeEsq, paredeDir);
//    getchar();

    /// === Condicao NOSLIP
    if( paredeDir==NOSLIP || paredeEsq==NOSLIP )
        return 0.0;

    /// === Condicao SLIP EM UMA PAREDE HORIZONTAL
    //Se estiver em uma parede superior
    if( paredeDir==SLIP ) {
        double k = M.pontosU[i+1][jU].valorDirichlet;
        termo = ( 2.0*k*M.beta )/( M.Re*M.dx[i] );
        return ( termo*V[i][j] )/( 1.0 + termo );
    }
    //Se estiver em uma parede inferior
    if( paredeEsq==SLIP ) {
        double k = M.pontosU[i][jU].valorDirichlet;
        termo = ( 2.0*k*M.beta )/( M.Re*M.dx[i] );
        return ( termo*V[i][j] )/( 1.0 + termo );
    }


    DeuProblema("VALOR NO SLIP VERTICAL: PROBLEMA %d %d\n", i, j);
    return 0.0;

    //return 0.0;
    //return 8*(1 + tanh(8*T-4))*X*X*(1-X)*(1-X);
}


/// ADiciona e remove particulas na interface (ela esta aqui, pois depende da malha nao-uniforme)
void AdicionaERemoveParticulas(MALHA *M, double TamanhoMaximo, int PermiteCelulasEmpty)
{
    CURVA *curva;
    PONTO *ponto;//, *proxPonto;//, *novoPonto;
    // double tamanho;
    // int i, j;
    int qtdPontos = 0;

//    return;

    /// ==== Verifica se chegou em um caso degenerado em que ha uma curav que ficou vazia e remove
    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {

        if( curva->pontoInicial==NULL ) {
            if( curva->curvaAnt==NULL ) { // Primeira curva da lista
                M->interface.curvaInicial = curva->proxCurva;
                M->interface.curvaInicial->curvaAnt = NULL;
                free(curva);
                break;
            }
            else {
                curva->curvaAnt->proxCurva = curva->proxCurva;
                curva->proxCurva->curvaAnt = curva->curvaAnt;
                free(curva);
                break;
            }
        }
    }

    /// ==== Removendo particulas
    // for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        
    //     LOOP_CURVA(curva) {
    //         ponto = loop_stuff.ponto;
    //         proxPonto = ponto->prox;

    //         // === nao posso remover nem o ponto inicial nem o final
    //         if( !(ponto->ant) || !(ponto->prox) )
    //             continue;

    //         // Usando a estrategia pra malha nao-uniforme
    //         EncontraCelulaDaParticula(M, ponto->x, ponto->y, &i, &j);
    //         TamanhoMaximo = (M->dx[i] < M->dy[j]) ? M->dx[i] : M->dy[j];
    //         TamanhoMaximo = 0.5 * 0.2*TamanhoMaximo;

    //         // Medindo a distancia entre ponto e pontoAnterior
    //         double distancia, x1, x2, y1, y2;
    //         x1 = ponto->x;
    //         y1 = ponto->y;
    //         x2 = ponto->ant->x;
    //         y2 = ponto->ant->y;
    //         distancia = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    //         if( distancia<TamanhoMaximo ) {
    //             ponto->ant->prox = ponto->prox;
    //             ponto->prox->ant = ponto->ant;
    //             free(ponto);
    //             continue;
    //         }

    //         // Medindo a distancia entre ponto e pontoProximo (apenas no penultimo ponto da curva)
    //         if( proxPonto->prox==NULL ) {
    //             x1 = ponto->x;
    //             y1 = ponto->y;
    //             x2 = ponto->prox->x;
    //             y2 = ponto->prox->y;
    //             distancia = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    //             if( distancia<TamanhoMaximo ) {
    //                 ponto->ant->prox = ponto->prox;
    //                 ponto->prox->ant = ponto->ant;
    //                 free(ponto);
    //                 continue;
    //             }
    //         }

    //     }
    // }

    /// Criando particulas
    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {

        if( curva->pontoInicial==NULL ) {
            DesenhaMalhaVTK(*M, 0);
            DeuProblema("AdicionaERemoveParticulas: CURVA SEM NENHUM PONTO\n\n");
        }

        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;

            qtdPontos++;

            /// === Usando a estrategia pra malha nao-uniforme
            /// === Escolho um espacamento ponderado entre dx e dy dependendo da direcao da interface
//            TamanhoMaximo = (M->dx[i] < M->dy[j]) ? M->dx[i] : M->dy[j];
//            TamanhoMaximo = 0.2*TamanhoMaximo;

            DivideSegmentoSeNecessario(M, ponto, TamanhoMaximo, PermiteCelulasEmpty);
        }
    }

    /// Removendo particulas travadas no outflow
    // for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {

    //     for( ponto=curva->nodeInicial; (ponto) && (ponto->prox) && (ponto->prox->prox); ponto=ponto->prox ) {
    //         PONTO *p1, *p2;
    //         p1 = ponto;
    //         p2 = ponto->prox->prox;

    //         //Quero apenas o caso em que os dois sao outflow
    //         if( !(p1->outflow) || !(p2->outflow)  )
    //             continue;

    //         if( !(p1->ant) || !(p2->prox) )
    //             continue;

    //         //Removo o ponto que esta entre eles
    //         PONTO *meio = ponto->prox;
    //         ponto->prox = meio->prox;
    //         meio->prox->ant = ponto;
    //         free(meio);

    //     }
    // }

//    ///Removendo particulas
//    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
//
//        for( ponto=curva->nodeInicial; (ponto) && (ponto->prox) && (ponto->prox->prox); ponto=ponto->prox ) {
//
////            PONTO *p1, *p2;
////            p1 = ponto;
////            p2 = ponto->prox->prox;
////
////            if( !(p1->ant) || !(p2->prox) )
////                continue;
//
//            // Usando a estrategia pra malha nao-uniforme
//            EncontraCelulaDaParticula(M, ponto->x, ponto->y, &i, &j);
//            TamanhoMaximo = (M->dx[i] < M->dy[j]) ? M->dx[i] : M->dy[j];
//            TamanhoMaximo = 0.2*TamanhoMaximo;
//
//            double distancia, x1, x2, y1, y2;
//            x1 = ponto->x;
//            y1 = ponto->y;
//            x2 = ponto->prox->prox->x;
//            y2 = ponto->prox->prox->y;
//            distancia = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
//            if( distancia<TamanhoMaximo ) {
//                PONTO *meio = ponto->prox;
//                ponto->prox = meio->prox;
//                meio->prox->ant = ponto;
//                free(meio);
//            }
//
//        }
//    }

//    EncontraCelulaDeTodasAsParticulas(M);


    return;
}

void DivideSegmentoSeNecessario(MALHA *M, PONTO *Ponto, double TamanhoMaximo, int PermiteCelulasEmpty)
{
    double tamanho;
    PONTO *proxPonto, *novoPonto;
    static int i, j;


    /// === Testando nova estrategia de determinar o tamanho maximo aqui dentro pra malha nao-uniforme
    EncontraCelulaDaParticula(M, Ponto->x, Ponto->y, &i, &j);
    proxPonto = Ponto->prox;
    static double vetX, vetY, soma, espacamento;
    vetX = fabs(proxPonto->x - Ponto->x);
    vetY = fabs(proxPonto->y - Ponto->y);
    soma = vetX + vetY;

    if( soma<1e-10 )
        return;

    vetX /= soma;
    vetY /= soma;

    // Giving a little extra weight to the smallest direction spacing
    espacamento = (M->dx[i]<M->dy[j]) ? 2.0*vetX*M->dx[i] + vetY*M->dy[j] : vetX*M->dx[i] + 2.0*vetY*M->dy[j];
    espacamento *= 0.5;
    TamanhoMaximo = 0.2*espacamento;

    

    if( isnan(vetX) || isinf(vetY) ) {
        DeuProblema("\n DivideSegmentosSeNecessario: divisao com problema. \n\n");
    }


    proxPonto = Ponto->prox;
    tamanho = sqrt( (proxPonto->x - Ponto->x)*(proxPonto->x - Ponto->x) + (proxPonto->y - Ponto->y)*(proxPonto->y - Ponto->y) );

    //Divide no meio
    if( tamanho>TamanhoMaximo ) {

        // Coordenadas da nova particula
        double x = 0.5*( Ponto->x + proxPonto->x );
        double y = 0.5*( Ponto->y + proxPonto->y );

        // Dividindo usando spline cubica
//        double x, y;
//        DivideSegmentoSpline(M, Ponto, &x, &y);

        EncontraCelulaDaParticula(M, x, y, &i, &j);

        // Se a celula for PAREDE, nao cria esta particula...
        if( M->celulasParede[i][j]==BOUNDARY && !PermiteCelulasEmpty ) {
            return;
//            int criterioDir, criterioEsq, criterioCima, criterioBaixo;
//            criterioEsq = M->pontosU[i][j].tipo==BOUNDARY;
//            criterioDir = M->pontosU[i+1][j].tipo==BOUNDARY;
//            criterioBaixo = M->pontosV[i][j].tipo==BOUNDARY;
//            criterioCima = M->pontosV[i][j+1].tipo==BOUNDARY;
//
//            // Apenas se tiver do lado de uma parede
//            if( criterioBaixo || criterioCima || criterioEsq || criterioDir )
//                return;
        }

	//Tratamento forcado de uma quina:
	//Nesta versao eu classifiquei as celulas BOUNDARY dentro de um bloco paraview.
	if (M->celulas[i][j] == BOUNDARY){
		return;

		//Outra opcao
		//y=M->y1;

		//printf("Dentro de AdicionaERemoveParticulas: Falha no if - criando particulas erroneamente\n");
		//printf("M->pontosU[i][j].tipo=%d, M->pontosV[i][j+1].tipo=%d\n", M->pontosU[i][j].tipo, M->pontosV[i][j+1].tipo);
		//printf("BOUNDARY=%d, M->celulasParede[i][j]=%d\n", BOUNDARY, M->celulasParede[i][j]);
		//exit(1);
	}

        novoPonto = (PONTO *)malloc( sizeof(PONTO) );
        novoPonto->ant = Ponto;
        novoPonto->prox = proxPonto;
        novoPonto->isNH = 0;
        novoPonto->isNode = 0;
        novoPonto->x = x;
        novoPonto->y = y;
        Ponto->prox = novoPonto;
        proxPonto->ant = novoPonto;
        novoPonto->fixo = (Ponto->fixo && proxPonto->fixo) ? 1 : 0;
        novoPonto->outflow = (Ponto->outflow && proxPonto->outflow);

        //Chama recursivamente pra ver se precisa dividir de novo
        DivideSegmentoSeNecessario(M, Ponto, TamanhoMaximo, PermiteCelulasEmpty);
        DivideSegmentoSeNecessario(M, novoPonto, TamanhoMaximo, PermiteCelulasEmpty);
    }

//     return;
}

// void DivideSegmentoSpline(MALHA *M, PONTO *PontoEscolhido, double *NovoX, double *NovoY)
// {
//     int i;
//     int QtdParticulas = 5;

//     // Andando pra tras nas particulas ate chegar na particula inicial
//     PONTO *p = PontoEscolhido;
//     for( i=0; i<QtdParticulas; i++ )
//         p = (p->ant) ? p->ant : p->pontoAntCurva;

//     PONTO *pInicial = p;

//     // Alocando memoria para os vetores
//     double *x, *y, *t_vec;
//     int totalPontos = 2*QtdParticulas + 1;
//     x = (double *)malloc( totalPontos*sizeof(double) );
//     y = (double *)malloc( totalPontos*sizeof(double) );
//     t_vec = (double *)malloc( totalPontos*sizeof(double) );


//     // Fazendo a parametrizacao com o comprimento de cada segmento
//     t_vec[0] = 0.0;
//     for( i=1; i<totalPontos; i++ ) {
//         PONTO *pProx = (p->prox) ? p->prox : p->pontoProxCurva;

//         t_vec[i] = t_vec[i-1] + sqrt( (p->x - pProx->x)*(p->x - pProx->x) + (p->y - pProx->y)*(p->y - pProx->y) );

//         p = pProx;
//     }



//     // Percorrendo os pontos que serao utilizdos na interpolacao splines e jogando em vetores
//     p = pInicial;


//     for( i=0; i<totalPontos; i++ ) {
//         x[i] = p->x;
//         y[i] = p->y;

//         p = (p->prox) ? p->prox : p->pontoProxCurva;
//     }

//     //Calculando os coeficientes das splines
//     POL3 *s_x, *s_y;
//     s_x = (POL3 *)malloc((totalPontos)*sizeof(POL3)); //alocando memoria para os polinomios
//     s_y = (POL3 *)malloc((totalPontos)*sizeof(POL3)); //alocando memoria para os polinomios

//     SplineCubicaNotAKnot(t_vec, x, s_x, totalPontos-1); //Calculando os polinomios
//     SplineCubicaNotAKnot(t_vec, y, s_y, totalPontos-1); //Calculando os polinomios


//     // Finalmente, calculando as coordenadas do novo ponto
//     int tEscolhido = QtdParticulas;
//     double t = 0.5*( t_vec[tEscolhido] + t_vec[tEscolhido+1] );
//     *NovoX = ValorSpline(t, t_vec, s_x, totalPontos-1);
//     *NovoY = ValorSpline(t, t_vec, s_y, totalPontos-1);


//     free(x);
//     free(y);
//     free(t_vec);

//     return;
// }

void EncontraCelulaDeTodasAsParticulas(MALHA *M)
{
    PONTO *ponto;
    CURVA *curva;

    for( curva=M->interface.curvaInicial; curva!=NULL; curva=curva->proxCurva ) {
        LOOP_CURVA(curva) {
            ponto = loop_stuff.ponto;
            EncontraCelulaDaParticula(M, ponto->x, ponto->y, &(ponto->i), &(ponto->j));
        }
    }
    return;
}

void SuavizaInterfaceTSUR(MALHA *M)
{
    CURVA *curva;
    PONTO *p;

    for( curva=M->interface.curvaInicial; curva; curva=curva->proxCurva ) {
        LOOP_CURVA(curva) {
            p = loop_stuff.ponto;

            TSUR(M, p);

            if( p->prox==NULL ) {
                curva->pontoInicial->x = p->x;
                curva->pontoInicial->y = p->y;
            }
        }
    }

}

void TSUR(MALHA *M, PONTO *P1)
{
    double   x0,x1,x2,x3; /* coordenadas x dos pontos p0,p1,p2 e p3  */
    double   y0,y1,y2,y3; /* coordenadas y dos pontos p0,p1,p2 e p3  */
    double   a0;          /* area do triangulo [p0,p1,p2]            */
    double   a1;          /* area do triangulo [p0,p2,p3]            */
    double   a;           /* area total (a0+a1)                      */
    double   l;           /* 1/3 |p0-p3|                             */
    double   h;           /* altura do trapezio                      */
    int i1Antes, j1Antes, i2Antes, j2Antes; /* coordenadas da malha dos pontos p1 e p2 */
    int i1Depois, j1Depois, i2Depois, j2Depois; /* coordenadas da malha dos pontos p1 e p2 */
    int fixo0, fixo1, fixo2, fixo3; //Pra impedir de mover particulas fixas (tipo inflow)

    PONTO *p;
    /* definindo os quatro pontos consecutivos da superficie */
    /* p0, p1, p2 e p3                                       */
    p = P1;
    x0 = p->x;
    y0 = p->y;
    fixo0 = p->fixo;

    p = p->prox;
    x1 = p->x;
    y1 = p->y;
    fixo1 = p->fixo;
    EncontraCelulaDaParticula(M, x1, y1, &i1Antes, &j1Antes);

    p = p->prox;
    x2 = p->x;
    y2 = p->y;
    fixo2 = p->fixo;
    EncontraCelulaDaParticula(M, x2, y2, &i2Antes, &j2Antes);

    p = p->prox;
    x3 = p->x;
    y3 = p->y;
    fixo3 = p->fixo;

    if( fixo0 || fixo1 || fixo2 || fixo3 )
        return;

    /// Impedindo modificao de particulas no canto inicial e final do dominio
    if( (i1Antes==0) || (i2Antes==0) )
        return;
    if( (i1Antes==M->Nx-1) || (i2Antes==M->Nx-1) )
        return;

//    int tsur_extremos = BuscaParametro(M, "tsur_extremos");
    int tsur_extremos = 1;

    if( tsur_extremos==0 ) {
        /// Impedindo modificao de particulas no canto final do dominio
        if( (j1Antes==0) || (j2Antes==0) )
            return;

        /// Impedindo modificao de particulas no canto final do dominio
        if( (j1Antes==M->Ny-1) || (j2Antes==M->Ny-1) )
            return;
    }

    /* coordenadas da malha dos pontos p1 e p2 */

    /* area do triangulo [p0,p1,p2] */
    /// OBS: na verdade falta multiplicar por meio pra virar a area. Depois multiplicamos por meio na definicao de h
    a0 = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);

    /* area do triangulo [p0,p2,p3] */
    /// OBS: na verdade falta multiplicar por meio pra virar a area. Depois multiplicamos por meio na definicao de h
    a1 = (x2-x0)*(y3-y0)-(x3-x0)*(y2-y0);

    /* soma das areas dos trinagulos */
    a  = a0+a1;

    if( fabs(a)<1e-8 )
        return;

    /* calculando 1/3 do comprimento do segmento p0<->p3 */
    l = sqrt((x3-x0)*(x3-x0)+(y3-y0)*(y3-y0))/3.0;

    /* calculando a altura do trapezio com os novos pontos */
    /* p1 e p2 que mantem a area original                  */
    // O calculo da altura depende se eh axissimetrico ou cartesiano 2D
    if( M->tipoCoord==AXI_CILINDRICO ) {
        double centroide_x = (x0 + x1)*(x0*y1 - x1*y0);
        centroide_x += (x1 + x2)*(x1*y2 - x2*y1);
        centroide_x += (x2 + x3)*(x2*y3 - x3*y2);
        centroide_x += (x3 + x0)*(x3*y0 - x0*y3);
        centroide_x /= (3.0*a); // Lembrando que a variavel "a" ta com a area duplicada

        // Teorema do centroide de Pappus
        // Volume = Area * (PerimetroDaRotacaoDoCentroide)
        double volume = (0.5*a)*2.0*M_PI*centroide_x;

        double tau_r = (x3 - x0)/(3.0*l);
        double tau_z = (y3 - y0)/(3.0*l);

        if( fabs(tau_z)>=1e-8 ) {
            double bhaskara_a = 1.0;
            double bhaskara_b = 1.2*( 2.0*(x0/tau_z) + 3.0*l*(tau_r/tau_z) );
            double bhaskara_c = -0.6*volume/(M_PI*l*tau_z);

            double bhaskara_delta = bhaskara_b*bhaskara_b - 4.0*bhaskara_a*bhaskara_c;
            double sol_1 = (- bhaskara_b + sqrt(bhaskara_delta))/(2.0*bhaskara_a);
            double sol_2 = (- bhaskara_b - sqrt(bhaskara_delta))/(2.0*bhaskara_a);

            // Escolhe a solucao que vai dar r1>=0 e r2>=0
            h = sol_1;
            double r1_temp = x0+(x3-x0-h*(y0-y3)/l)/3.0;
            double r2_temp = x0+(2.0*(x3-x0)-h*(y0-y3)/l)/3.0;
            if( r1_temp<0.0 || r2_temp<0.0 )
                h = sol_2;
        }
        else {
            h = (0.5*volume/M_PI)/(2.0*x0*l + 3.0*l*l*tau_r);
        }


    }
    else { // CARTESIANO
        /// Tem uma multiplicacao por meio a mais aqui, pra compensar a area que nao foi multiplicada la em cima
        h = 0.25*a/l;
    }

//    if( isnan(x1) || isnan(y1) || isnan(x2) || isnan(y2) )
//        DeuProblema("1. vish nan\n");

    /* redefinindo o ponto p1 */
    x1 = x0+(x3-x0-h*(y0-y3)/l)/3.0;
    y1 = y0+(y3-y0-h*(x3-x0)/l)/3.0;

    /* redefinindo o ponto p2 */
    x2 = x0+(2.0*(x3-x0)-h*(y0-y3)/l)/3.0;
    y2 = y0+(2.0*(y3-y0)-h*(x3-x0)/l)/3.0;

//    if( isnan(x1) || isnan(y1) || isnan(x2) || isnan(y2) )
//        DeuProblema("2. vish nan: %e\n", l);

    // Denominador igual a zero ali em cima...
    if( fabs(l)<1e-13 )
        return;

    if( x1<=M->xMin ) x1 = M->xMin + 1e-8;
    if( x1>=M->xMax ) x1 = M->xMax - 1e-8;
    if( y1<=M->yMin ) y1 = M->yMin + 1e-8;
    if( y1>=M->yMax ) y1 = M->yMax - 1e-8;

    if( x2<=M->xMin ) x2 = M->xMin + 1e-8;
    if( x2>=M->xMax ) x2 = M->xMax - 1e-8;
    if( y2<=M->yMin ) y2 = M->yMin + 1e-8;
    if( y2>=M->yMax ) y2 = M->yMax - 1e-8;



    /* coordenadas da malha dos pontos p1 e p2 modificados */
    EncontraCelulaDaParticula(M, x1, y1, &i1Depois, &j1Depois);
    EncontraCelulaDaParticula(M, x2, y2, &i2Depois, &j2Depois);

    // Verifando se p1 mudou de tipo de celula ou se p2 mudou de tipo de celula
    TIPO_CELULA celula1Antes, celula1Depois;
    TIPO_CELULA celula2Antes, celula2Depois;

    celula1Antes = M->celulas[i1Antes][j1Antes];
    celula1Depois = M->celulas[i1Depois][j1Depois];

    celula2Antes = M->celulas[i2Antes][j2Antes];
    celula2Depois = M->celulas[i2Depois][j2Depois];

//    if( celula1Antes==EMPTY || celula2Antes==EMPTY ) {
////        DesenhaMalhaVTK(*M, 0);
////        ImprimeInterfaceVTK(*M, 900000000);
//        PrintDebug("AAAAAAAAAAAAAAAA\n\n");
//    }

//    p = P1;
//    p = (p->prox) ? p->prox : p->pontoProxCurva;
//    if( fabs(p->x - 4.001970e-01)<1e-6 ) {
//        PrintDebug("ENTROU AQUI!!!! P1\n");
//    }
//
//    p = (p->prox) ? p->prox : p->pontoProxCurva;
//    if( fabs(p->x - 4.001970e-01)<1e-6 ) {
//        PrintDebug("ENTROU AQUI!!!! P2\n");
//    }

    /// Impedindo que a classificacao de celulas seja alterada
    if( celula1Antes!=celula1Depois )
        return;
    if( celula2Antes!=celula2Depois )
        return;


    /// Impedindo que a particula 1 atravesse uma boundary
    if( (i1Depois>i1Antes) && M->pontosU[i1Depois][j1Antes].tipo==BOUNDARY )
        return;
    else if( (i1Depois<i1Antes) && M->pontosU[i1Antes][j1Antes].tipo==BOUNDARY )
        return;
    if( (j1Depois>j1Antes) && M->pontosV[i1Antes][j1Depois].tipo==BOUNDARY )
        return;
    else if( (j1Depois<j1Antes) && M->pontosV[i1Antes][j1Antes].tipo==BOUNDARY )
        return;

    /// Impedindo que a particula 2 atravesse uma boundary
    if( (i2Depois>i2Antes) && M->pontosU[i2Depois][j2Antes].tipo==BOUNDARY )
        return;
    else if( (i2Depois<i2Antes) && M->pontosU[i2Antes][j2Antes].tipo==BOUNDARY )
        return;
    if( (j2Depois>j2Antes) && M->pontosV[i2Antes][j2Depois].tipo==BOUNDARY )
        return;
    else if( (j2Depois<j2Antes) && M->pontosV[i2Antes][j2Antes].tipo==BOUNDARY )
        return;


    /* atualizando as coordenadas de p1 */
    p = P1;
    p = p->prox;
    p->x = x1;
    p->y = y1;

    /* atualizando as coordenadas de p2 */
    p = p->prox;
    p->x = x2;
    p->y = y2;




    return;
}





///FUNCOES DE VISUALIZACAO
void DesenhaMalhaGnuplot(MALHA Malha)
{
//    FILE *arq;
//    int i, j;
//    double x, y;
//
//    arq = fopen("ArquivosPlot/ScriptMalha.txt", "wt");
//
//    fprintf(arq, "unset key\n");
//    fprintf(arq, "set tic scale 0\n");
//    fprintf(arq, "set xlabel \"x\"\n");
//    fprintf(arq, "set ylabel \"y\"\n");
//    fprintf(arq, "set size square\n");
//    //fprintf(arq, "set ratio -1\n");
//    //fprintf(arq, "#set cbrange [-0.5:1.5]\n");
//    //fprintf(arq, "#set cblabel \"Escala\"\n");
//    //fprintf(arq, "#unset cbtics\n");
//    fprintf(arq, "set xrange [%lf:%lf]\n", Malha.xMin, Malha.xMax);
//    fprintf(arq, "set yrange [-%lf:%lf]\n", Malha.yMin, Malha.yMax);
//    //fprintf(arq, "set view map\n");
//    fprintf(arq, "set term pdf\n");
//    fprintf(arq, "set output \"ArquivosPlot/PlotMalha.pdf\"\n");
//
//    x = Malha.xMin;
//    for( i=0; i<=Malha.Nx; i++ ) {
//        fprintf(arq, "set arrow from %lf, %lf to %lf, %lf nohead\n", x, Malha.yMin, x, Malha.yMax);
//        if( i!=Malha.Nx )
//            x += Malha.dx[i];
//    }
//
//    y = Malha.yMin;
//    for( j=0; j<=Malha.Ny; j++ ) {
//        fprintf(arq, "set arrow from %lf, %lf to %lf, %lf nohead\n", Malha.xMin, y, Malha.xMax, y);
//        if( j!=Malha.Ny )
//            y += Malha.dy[j];
//    }
//
//
//    //fprintf(arq, "set object ellipse center 0.5,0.5 size 0.50,0.50 fillstyle empty lw 2.0 fc rgb \"blue\"\n");
//
//    fprintf(arq, "plot -50\n");
//    fprintf(arq, "unset output\n");
//
//    fclose(arq);
//
//    system("gnuplot ArquivosPlot/ScriptMalha.txt");
//    return;
}



void DesenhaMalhaVTK(MALHA Malha, int n)
{
    FILE *arq;
    char nomeArq[500];
    int i, j;

//    return;

    sprintf(nomeArq, "ArquivosPlot/VTK/%s/MalhaComputacional%d.vtk",  Malha.nomeArquivo, n);


    arq = fopen(nomeArq, "wt");

    fprintf(arq, "# vtk DataFile Version 2.0\n");
    fprintf(arq, "Malha computacional\n");
    fprintf(arq, "ASCII\n");
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



    fprintf(arq, "CELL_DATA %d\n", Malha.Nx*Malha.Ny);
    fprintf(arq, "SCALARS celulas float 3\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for( j=0; j<Malha.Ny; j++ ) {
        for( i=0; i<Malha.Nx; i++ ) {
            if( Malha.celulas[i][j]==FULL /*&& (Malha.celulasParede[i][j]!=SURFACE)*/ )
                fprintf(arq, "%lf %lf %lf\n", 155.0/255.0, 221/255.0, 255/255.0);
            else if( Malha.celulas[i][j]==EMPTY )
                fprintf(arq, "0.9 0.9 0.9\n");
            else if( Malha.celulas[i][j]==SURFACE )
                fprintf(arq, "%lf %lf %lf\n", 255.0/255.0, 0/255.0, 0/255.0);
            else if( Malha.celulas[i][j]==MUDANCA_TOPOL )
                fprintf(arq, "%lf %lf %lf\n", 255.0/255.0, 179.0/255.0, 71.0/255.0);
            else if( Malha.celulas[i][j]==BOUNDARY )
                fprintf(arq, "%lf %lf %lf\n", 176.0/255.0, 95.0/255.0, 60.0/255.0);
            else if( Malha.celulasParede[i][j]==SURFACE )
                fprintf(arq, "%lf %lf %lf\n", 255.0/255.0, 127.0/255.0, 80.0/255.0);
            else
                DeuProblema("1. PROBLEMA NO DESENHO DA MALHA\n");
        }
    }

    fprintf(arq, "\nSCALARS pressao float 3\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for( j=0; j<Malha.Ny; j++ ) {
        for( i=0; i<Malha.Nx; i++ ) {
            if( Malha.pontosP[i][j].tipo==FULL )
                fprintf(arq, "%lf %lf %lf\n", 155.0/255.0, 221/255.0, 255/255.0);
            else if( Malha.pontosP[i][j].tipo==EMPTY )
                fprintf(arq, "0.9 0.9 0.9\n");
            else if( Malha.pontosP[i][j].tipo==SURFACE )
                fprintf(arq, "%lf %lf %lf\n", 255.0/255.0, 0/255.0, 0/255.0);
            else if( Malha.pontosP[i][j].tipo==BOUNDARY )
                fprintf(arq, "%lf %lf %lf\n", 176.0/255.0, 95.0/255.0, 60.0/255.0);
            else
                DeuProblema("2. PROBLEMA NO DESENHO DA MALHA %d\n");
        }
    }

    fprintf(arq, "\nSCALARS celulasParede float 3\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    for( j=0; j<Malha.Ny; j++ ) {
        for( i=0; i<Malha.Nx; i++ ) {
            if( Malha.celulasParede[i][j]==FULL )
                fprintf(arq, "%lf %lf %lf\n", 155.0/255.0, 221/255.0, 255/255.0);
            else if( Malha.celulasParede[i][j]==EMPTY )
                fprintf(arq, "0.9 0.9 0.9\n");
            else if( Malha.celulasParede[i][j]==SURFACE )
                fprintf(arq, "%lf %lf %lf\n", 255.0/255.0, 0/255.0, 0/255.0);
            else if( Malha.celulasParede[i][j]==BOUNDARY )
                fprintf(arq, "%lf %lf %lf\n", 176.0/255.0, 95.0/255.0, 60.0/255.0);
            else
                DeuProblema("3. PROBLEMA NO DESENHO DA MALHA %d\n");
        }
    }

    fclose(arq);

    PrintDebug("\nDesenhou a Malha!\n");
    return;
}



///Funcoes auxiliares pra usar a pilha do flood fill
void PilhaPush(PILHA *P, int i, int j)
{
    ITEM_PILHA *novoItem;

    //PROBLEMA
    if( P==NULL )
        return;

    novoItem = (ITEM_PILHA*)malloc(sizeof(ITEM_PILHA));
    novoItem->i = i;
    novoItem->j = j;
    novoItem->prox = NULL;

    if( P->topo==NULL )
        P->topo = novoItem;
    else {
        novoItem->prox = P->topo;
        P->topo = novoItem;
    }

    return;
}

void PilhaPop(PILHA *P, int *i, int *j)
{
    ITEM_PILHA *temp;

    if( P==NULL || P->topo==NULL)
        return;

    *i = P->topo->i;
    *j = P->topo->j;

    //Removendo
    temp = P->topo;

    P->topo = P->topo->prox;
    free(temp);

    return;
}


///Usada pra criar a pasta dos arquivos vtk
int SubstituiString(char *Str, const char *Orig, char *Sub)
{
    static char temp[4096];
    static char buffer[4096];
    char *p;

    strcpy(temp, Str);

    if(!(p = strstr(temp, Orig)))  // Is 'orig' even in 'temp'?
        return 0;

    strncpy(buffer, temp, p-temp); // Copy characters from 'temp' start to 'orig' str
    buffer[p-temp] = '\0';

    sprintf(buffer + (p - temp), "%s%s", Sub, p + strlen(Orig));
    sprintf(Str, "%s", buffer);
    return 1;
}

void StringMaiuscula(char *String)
{
    int i;

    if( String==NULL )
        return;

    for( i=0; String[i]!='\0'; i++ )
        String[i] = toupper(String[i]);

    return;
}

void ArquivoMaiusculoTemporario(FILE *ArqOriginal, FILE *ArqNovo)
{
    int c;
    int continuar = 1;

    while ((c=fgetc(ArqOriginal))!=EOF) {
        // Se encontrar o caractere / (inicio de um path), interrompe a transformação toupper
        if (c == '/') {
            continuar = 0;
        // Se encontra o caractere "espaço" (fim de um path), reativa a transformação toupper
        } else if (c== ' ') {
            continuar = 1;
        }

        if (continuar == 1) {
            fprintf(ArqNovo, "%c", toupper(c));
        } else {
            fprintf(ArqNovo, "%c", c);
        }
    }
}

//Funcao criado por mim, Irineu, afim de classificar as celulas BOUNDARY
void ArmazenaQuinasBlocoParaview(MALHA *M, BLOCO_PARAVIEW *PrimBloco, int NumFaces){

	//double **mat;	
        int i;

	/*
	//printf("vamos entrar na allocation\n");
	mat = (double **)malloc(sizeof(double *)*(NumFaces+1));
	for (i=0; i<NumFaces+1; i++){
		mat[i] = (double *)malloc(sizeof(double)*3);
	}

	//printf("ja alocamos\n");

	M->PontosBlocoParaview = mat;
	*/

	//printf("NumFaces=%d\n", NumFaces);getchar();
	//BLOCO_PARAVIEW novoBloco=NULL;
	//BLOCO_PARAVIEW *bloco;
    	//bloco = (BLOCO_PARAVIEW *)malloc( sizeof(BLOCO_PARAVIEW) );
    	//novoBloco->matrizPontos = (double **)AlocaMatriz(NumFaces + 1, 2, sizeof(double), 0);

    	//Escrevendo os pontos
    	//for( bloco=PrimBloco; bloco!=NULL; bloco=bloco->prox ) {
        //for( i=0; i<NumFaces+1; i++ ){
		//mat[i][0] = PrimBloco->matrizPontos[i][0];
		//mat[i][1] = PrimBloco->matrizPontos[i][1];
		//printf("M->x[i]=%f, mat[%d][%d]=%lf, mat[%d][%d]=%lf \n", M->x[i], i, 0, M->PontosBlocoParaview[i][0], i, 0, M->PontosBlocoParaview[i][1]);
	//}
    	//}
	M->x1 = PrimBloco->matrizPontos[0][0];
	M->y1 = PrimBloco->matrizPontos[0][1];

	M->x2 = PrimBloco->matrizPontos[1][0];
	M->y2 = PrimBloco->matrizPontos[1][1];

	M->x3 = PrimBloco->matrizPontos[2][0];
	M->y3 = PrimBloco->matrizPontos[2][1];

	M->x4 = PrimBloco->matrizPontos[3][0];
	M->y4 = PrimBloco->matrizPontos[3][1];

	/*
	for (i=0; i<NumFaces+1; i++){
		free(mat[i]);
	}
	free(mat);
	//free(bloco);
	*/	
}


void AdicionaBlocoParaview(MALHA *M, BLOCO_PARAVIEW **PrimBloco, int NumFaces, int *Faces)
{
    BLOCO_PARAVIEW *novoBloco;
    double valorInicial = 0.0;

    novoBloco = (BLOCO_PARAVIEW *)malloc( sizeof(BLOCO_PARAVIEW) );
    novoBloco->matrizPontos = (double **)AlocaMatriz(NumFaces + 1, 2, sizeof(double), &valorInicial);
    novoBloco->numFaces = NumFaces;
    novoBloco->prox = *PrimBloco;
    *PrimBloco = novoBloco;



    //Colocando os valores dos pontos
    int i;
    int numPontos = 0;
    double ultimoX, ultimoY;
    for( i=0; i<NumFaces; i++ ) {

        PrintDebug("PONTO %d\n", i);


        // Encontrando o bloco com o id desta face
        BLOCO *blocoBoundary;
        for( blocoBoundary=M->primBloco; blocoBoundary!=NULL; blocoBoundary=blocoBoundary->prox ) {

            if( blocoBoundary->id==Faces[i] )
                break;
        }

        if( blocoBoundary==NULL )
            DeuProblema("AdicionaBlocoParaview: BLOCO NAO ENCONTRADO\n");

        if( i==0 ) {
            novoBloco->matrizPontos[numPontos][0] = blocoBoundary->x1;
            novoBloco->matrizPontos[numPontos++][1] = blocoBoundary->y1;

            novoBloco->matrizPontos[numPontos][0] = blocoBoundary->x2;
            novoBloco->matrizPontos[numPontos++][1] = blocoBoundary->y2;

            ultimoX = blocoBoundary->x2;
            ultimoY = blocoBoundary->y2;
        }
        else {

            PrintDebug("ULTIMO: %lf %lf\n", ultimoX, ultimoY);
            PrintDebug("p1: %lf %lf\n", blocoBoundary->x1, blocoBoundary->y1);
            PrintDebug("p2: %lf %lf\n", blocoBoundary->x2, blocoBoundary->y2);

            if( ZERO(blocoBoundary->x1 - ultimoX) && ZERO(blocoBoundary->y1 - ultimoY) ) {
                novoBloco->matrizPontos[numPontos][0] = blocoBoundary->x2;
                novoBloco->matrizPontos[numPontos++][1] = blocoBoundary->y2;
                ultimoX = blocoBoundary->x2;
                ultimoY = blocoBoundary->y2;
            }
            else if( ZERO(blocoBoundary->x2 - ultimoX) && ZERO(blocoBoundary->y2 - ultimoY) ) {
                novoBloco->matrizPontos[numPontos][0] = blocoBoundary->x1;
                novoBloco->matrizPontos[numPontos++][1] = blocoBoundary->y1;
                ultimoX = blocoBoundary->x1;
                ultimoY = blocoBoundary->y1;
            }
            else
                DeuProblema("AdicionaBlocoParaview: nao encontrou o proximo ponto.\n");


        }

    }

    PrintDebug("\n\n\n");
    for( i=0; i<=NumFaces; i++ )
        PrintDebug("%lf %lf\n", novoBloco->matrizPontos[i][0], novoBloco->matrizPontos[i][1]);




    return;
}

void EscreveArquivoBlocosParaview(BLOCO_PARAVIEW *PrimBloco, MALHA M)
{
    FILE *arqBlocosParaview;
    BLOCO_PARAVIEW *bloco;

    char nomeArq[500];
    sprintf(nomeArq, "ArquivosPlot/VTK/%s/blocos_paraview.vtk", M.nomeArquivo);

    arqBlocosParaview = fopen(nomeArq, "wt");
    fprintf(arqBlocosParaview, "# vtk DataFile Version 2.0\n");
    fprintf(arqBlocosParaview, "Structured Grid Example\n");
    fprintf(arqBlocosParaview, "ASCII\n");
    fprintf(arqBlocosParaview, "DATASET POLYDATA\n");


    // Contando a quantidade de pontos e blocos
    int qtdPontos = 0, qtdBlocos = 0;
    for( bloco=PrimBloco; bloco!=NULL; bloco=bloco->prox ) {
        qtdPontos += bloco->numFaces;
        qtdBlocos++;
    }
    fprintf(arqBlocosParaview, "POINTS %d float\n", qtdPontos);

    //Escrevendo os pontos
    for( bloco=PrimBloco; bloco!=NULL; bloco=bloco->prox ) {
        int i;
        for( i=0; i<bloco->numFaces; i++ )
            fprintf(arqBlocosParaview, "%lf %lf 0.0\n", bloco->matrizPontos[i][0], bloco->matrizPontos[i][1]);
    }

    //Escrevendo os poligonos
    fprintf(arqBlocosParaview, "POLYGONS %d %d\n", qtdBlocos, qtdPontos + qtdBlocos);
    int cont = 0;
    for( bloco=PrimBloco; bloco!=NULL; bloco=bloco->prox ) {
        fprintf(arqBlocosParaview, "%d", bloco->numFaces);
        int i;
        for( i=0; i<bloco->numFaces; i++ )
            fprintf(arqBlocosParaview, " %d", cont++);
        fprintf(arqBlocosParaview, "\n");
    }

    fclose(arqBlocosParaview);
}

// Funcao usada em simulacoes de queda de gota pra calcular a velocidade dela
// A velocidade calculada eh na ponta inferior da gota
double CalculaVelocidadeDaGota(MALHA *M, double **U, double **V, int Passo)
{
    int i, j, qtdPontos;
    PONTO *p, *pEscolhido = NULL;
    double menorDistancia, distancia, mediaY;

    if( M->interface.curvaInicial==NULL )
        return 0.0;

    if( (Passo-1)%100!=0 )
        return 0.0;

    // Abrindo o arquivo de saida
    static FILE *arqVeloc = NULL;

    if( arqVeloc==NULL ) {
        char nomeArq[500];
        sprintf(nomeArq, "ArquivosPlot/VTK/%s/velocidadeGota.txt", M->nomeArquivo);
        arqVeloc = fopen(nomeArq, "wt");
    }


    mediaY = 0.0;
    qtdPontos = 0;
    LOOP_CURVA(M->interface.curvaInicial) {
        p = loop_stuff.ponto;

        mediaY += p->y;
        qtdPontos++;
    }
    mediaY /= qtdPontos;

    // Encontrando a particula que esta mais no centro do dominio do lado de baixo
//    menorDistancia = 1e+10;
//    for( p=M->interface.curvaInicial->nodeInicial; p; p=p->prox ) {
//        if( p->y >= mediaY )
//            continue;
//
//        if( fabs(p->x - centroDominio)<menorDistancia ) {
//            menorDistancia = fabs(p->x - centroDominio);
//            pEscolhido = p;
//        }
//    }

    //Encontrando a particula que esta mais embaixo da gota
    double menorY = 1000.0;
    LOOP_CURVA(M->interface.curvaInicial) {
        p = loop_stuff.ponto;

        if( p->y < menorY ) {
            menorY = p->y;
            pEscolhido = p;
        }
    }

    // Encontrando o ponto da malha mais prox dessa particula
    menorDistancia = 1e+10;
    double x, y;
    int iEscolhido = 0, jEscolhido = 0;
    for( i=0; i<M->Nx; i++ ) {
        for( j=0; j<=M->Ny; j++ ) {

            // Pegando apenas pontos FULL
            if( M->pontosV[i][j].tipo!=FULL )
                continue;

            x = 0.5*(M->x[i] + M->x[i+1]);
            y = M->y[j];
            distancia = sqrt( (pEscolhido->x-x)*(pEscolhido->x-x) + (pEscolhido->y-y)*(pEscolhido->y-y) );

            if( distancia < menorDistancia ) {
                menorDistancia = distancia;
                iEscolhido = i;
                jEscolhido = j;
            }
        }
    }

    if( (Passo-1)%100==0 ) {
        fprintf(arqVeloc, "%lf %lf %lf %lf\n", M->dt*Passo, pEscolhido->x, pEscolhido->y, V[iEscolhido][jEscolhido]);
        fflush(arqVeloc);
    }

    return 0.0;
}






/// ============= FUNCOES QUAD TREE
// void CriaMalhaQuadTree(MALHA *M)
// {
//     CELULA_QUAD *celula;
//     CURVA *curva;
//     PONTO *ponto_interface;
//     PILHA_GERAL pilha;

//     pilha.prim = NULL;

//     celula = NovaCelulaQuad(NULL, M->xMin, M->xMax, M->yMin, M->yMax);
//     celula->qtdPontos = 0;
//     celula->ponto_celula = NULL;
//     M->qtdCelulas = 1;
//     M->quad_celula = celula;

//     /// === Colocando todos os pontos na lista da primeira celula
//     /// === Vou ignorar o primeiro ponto, pois nas curvas ciclicas eh igual ao primeiro
//     for( curva=M->interface.curvaInicial; curva; curva=curva->proxCurva ) {
//         for( ponto_interface=curva->nodeInicial->prox; ponto_interface; ponto_interface=ponto_interface->prox ) {
//             AdicionaPontoQuadTree(M, celula, ponto_interface, 2);
//         }
//     }

//     /// === Adicionando a celula inicial na pilha
//     PushPilha(&pilha, celula);

//     /// === Enquanto tiver celulas na pilha, continua...
// //    int i=0;
//     while( pilha.prim ) {
//         celula = PopPilha(&pilha);


//         /// === Se a celula tiver no maximo Max pontos, nao precisa fazer nada
//         if( celula->qtdPontos<=1 )
//             continue;

//         /// === Vamos dividir a celula em quatro
//         DivideCelulaQuad(M, celula);

//         /// === Adiciona as 4 novas celulas na pilha
//         PushPilha(&pilha, celula->cima_esq);
//         PushPilha(&pilha, celula->cima_dir);
//         PushPilha(&pilha, celula->baixo_esq);
//         PushPilha(&pilha, celula->baixo_dir);
//     }

//     return;
// }

// void DivideCelulaQuad(MALHA *M, CELULA_QUAD *C)
// {
//     static double meioX, meioY;
//     static PONTO_CELULA_QUAD *ponto, *proxPonto;

//     meioX = 0.5*(C->xMin + C->xMax);
//     meioY = 0.5*(C->yMin + C->yMax);
//     C->cima_esq = NovaCelulaQuad(C, C->xMin, meioX, meioY, C->yMax);
//     C->cima_dir = NovaCelulaQuad(C, meioX, C->xMax, meioY, C->yMax);
//     C->baixo_esq = NovaCelulaQuad(C, C->xMin, meioX, C->yMin, meioY);
//     C->baixo_dir = NovaCelulaQuad(C, meioX, C->xMax, C->yMin, meioY);

//     M->qtdCelulas += 3;

//     // Distribuindo os pontos entre as quatro celulas
//     for( ponto=C->ponto_celula; ponto; ponto=proxPonto ) {
//         proxPonto = ponto->prox;

//         // Verificando em qual das quatro celulas este ponto esta
//         if( C->cima_esq ) {
//             if( ponto->ponto->x<meioX ) {
//                 if( ponto->ponto->y<meioY )
//                     MovePontoQuadTree(C->baixo_esq, ponto);
//                 else
//                     MovePontoQuadTree(C->cima_esq, ponto);
//             }
//             else {
//                 if( ponto->ponto->y<meioY )
//                     MovePontoQuadTree(C->baixo_dir, ponto);
//                 else
//                     MovePontoQuadTree(C->cima_dir, ponto);
//             }
//         }
//     }

//     return;
// }

// CELULA_QUAD *NovaCelulaQuad(CELULA_QUAD *Pai, double xMin, double xMax, double yMin, double yMax)
// {
//     CELULA_QUAD *celula;

//     celula = (CELULA_QUAD*)malloc(sizeof(CELULA_QUAD));
//     celula->xMin = xMin;
//     celula->xMax = xMax;
//     celula->yMin = yMin;
//     celula->yMax = yMax;

//     celula->pai = Pai;

//     celula->cima_esq = NULL;
//     celula->cima_dir = NULL;
//     celula->baixo_esq = NULL;
//     celula->baixo_dir = NULL;

//     celula->ponto_celula = NULL;
//     celula->qtdPontos = 0;

//     return celula;
// }

// void DesenhaMalhaQuadVTK(MALHA M, int n)
// {
//     FILE *arq;
//     char nomeArq[500];
//     int i;

//     sprintf(nomeArq, "ArquivosPlot/VTK/%s/MalhaQuad%d.vtk",  M.nomeArquivo, n);


//     arq = fopen(nomeArq, "wt");

//     fprintf(arq, "# vtk DataFile Version 2.0\n");
//     fprintf(arq, "Malha computacional\n");
//     fprintf(arq, "ASCII\n");
//     fprintf(arq, "DATASET POLYDATA\n");

//     fprintf(arq, "POINTS %d float\n", 4*(M.qtdCelulas));
//     ImprimePontosQuadTree(arq, M.quad_celula);

//     fprintf(arq, "POLYGONS %d %d\n", M.qtdCelulas, 5*(M.qtdCelulas));
//     for( i=0; i<M.qtdCelulas; i++ )
//         fprintf(arq, "4 %d %d %d %d\n", 4*i, 4*i + 1, 4*i + 2, 4*i + 3);

//     fclose(arq);
//     return;
// }

// void ImprimePontosQuadTree(FILE *Arq, CELULA_QUAD *C)
// {

//     if( C->cima_esq ) {
//         ImprimePontosQuadTree(Arq, C->cima_esq);
//         ImprimePontosQuadTree(Arq, C->cima_dir);
//         ImprimePontosQuadTree(Arq, C->baixo_esq);
//         ImprimePontosQuadTree(Arq, C->baixo_dir);
//         return;
//     }

//     fprintf(Arq, "%lf %lf %lf\n", C->xMin, C->yMin, 0.0);
//     fprintf(Arq, "%lf %lf %lf\n", C->xMax, C->yMin, 0.0);
//     fprintf(Arq, "%lf %lf %lf\n", C->xMax, C->yMax, 0.0);
//     fprintf(Arq, "%lf %lf %lf\n", C->xMin, C->yMax, 0.0);
//     return;
// }

// void MovePontoQuadTree(CELULA_QUAD *C, PONTO_CELULA_QUAD *P)
// {
//     P->prox = C->ponto_celula;
//     C->ponto_celula = P;
//     (C->qtdPontos)++;
//     return;
// }

// void AdicionaPontoQuadTree(MALHA *M, CELULA_QUAD *C, PONTO *PontoInt, int MaxPontos)
// {
//     PONTO_CELULA_QUAD *novoP;
// //    static double meioX, meioY;

//     novoP = (PONTO_CELULA_QUAD*)malloc(sizeof(PONTO_CELULA_QUAD));
//     novoP->ponto = PontoInt;
//     novoP->prox = C->ponto_celula;
//     C->ponto_celula = novoP;
//     (C->qtdPontos)++;
//     return;

//     /// === Se esta celula for folha, mas ja tiver o maximo de pontos, vou quebrar ela
// //    if( !(C->cima_esq) && (C->qtdPontos==MaxPontos) ) {
// //        DivideCelulaQuad(M, C, MaxPontos);
// //        AdicionaPontoQuadTree(M, C, PontoInt, MaxPontos);
// //        return;
// //    }
// //
// //    / === Se essa celula nao for uma folha, chama recursivamente. Soh vou adicionar nas folhas
// //
// }



double BuscaParametro(MALHA *M, const char *NomeParametro)
{
    ITEM_PILHA_GERAL *item;
    char temp[300];

    // Jogando pra maiusculo
    int i;
    for( i=0; NomeParametro[i]!='\0'; i++ )
        temp[i] = toupper( NomeParametro[i] );
    temp[i] = '\0';

    for( item=M->parametrosOpcionais.prim; item!=NULL; item=item->prox ) {
        PARAMETRO_OPCIONAL *param = (PARAMETRO_OPCIONAL *)item->valor;
        if( !strcmp(temp, param->nomeParametro) )
            return param->valorParametro;
    }

    DeuProblema("NAO FOI ENCONTRADO O PARAMETRO OPCIONAL %s\n", temp);
    return 0.0;
}

// Calcula a área do fluido através da fórmula de Gauss
float calcAreaFluido(MALHA malha) {
    CURVA *curva;
    PONTO *ponto;
    float area, areaTotal = 0;

    for (curva = malha.interface.curvaInicial; curva != NULL; curva = curva->proxCurva) {
        area = 0;
	ponto = curva->pontoInicial;
        do {
            if (ponto->prox != NULL) {
                area += ponto->x * ponto->prox->y - ponto->prox->x * ponto->y;
            } else {
                area += ponto->x * curva->pontoInicial->y - curva->pontoInicial->x * ponto->y;
            }
            ponto = ponto->prox;
        } while(ponto != curva->pontoInicial);
        area = fabs(area * 1/2);
        areaTotal += area;
    }

    return 2 * areaTotal;
}

// Calcula o spread máximo do fluido
float calcSpread(MALHA malha) {
    CURVA *curva;
    PONTO *ponto;
    float xmin = malha.xMax, xmax = malha.xMin;

    for (curva = malha.interface.curvaInicial; curva != NULL; curva = curva->proxCurva) {
        ponto = curva->pontoInicial;
        do {
            if (ponto->x < xmin) {
                xmin = ponto->x;
            } else if (ponto->x > xmax) {
                xmax = ponto->x;
            }
            ponto = ponto->prox;
        } while(ponto != curva->pontoInicial);
    }

    return 2 * (xmax - xmin);
}

// Calcula a altura máxima do fluido que passa pelo buraco
float calcAltura(MALHA malha) {
    CURVA *curva;
    PONTO *ponto;
    float ymin = malha.yMax, ymax = malha.yMin;

    for (curva = malha.interface.curvaInicial; curva != NULL; curva = curva->proxCurva) {
        ponto = curva->pontoInicial;
        do {
            if (ponto->y < ymin) {
                ymin = ponto->y;
            } else if (ponto->y > ymax) {
                ymax = ponto->y;
            }
            ponto = ponto->prox;
        } while(ponto != curva->pontoInicial);
    }

    return (ymax - ymin);
}

void fprintfVetorFormatado(FILE *arq, float *vet, int tam) {
    fprintf(arq, "[");
    for (int i = 0; i < tam; i++) {
        if (i != tam-1) {
            fprintf(arq, "%f, ", vet[i]);
        } else {
            fprintf(arq, "%f", vet[i]);
        }
    }
    fprintf(arq, "]");
}

void pontoEsquerda(MALHA malha, int *i, int *j) {
    CURVA *curva;
    PONTO *ponto, *ponto_xmin;

    ponto_xmin = malha.interface.curvaInicial->pontoInicial;
    for (curva = malha.interface.curvaInicial; curva != NULL; curva = curva->proxCurva) {
        ponto = curva->pontoInicial;
        do {
            if (ponto->x < ponto_xmin->x) {
                ponto_xmin = ponto;
            } 
            ponto = ponto->prox;
        } while(ponto != curva->pontoInicial);
    }

    EncontraCelulaDaParticula(&malha, ponto_xmin->x, ponto_xmin->y, i, j);

    return;
}

void pontoDireita(MALHA malha, int *i, int *j) {
    CURVA *curva;
    PONTO *ponto, *ponto_xmax;

    ponto_xmax = malha.interface.curvaInicial->pontoInicial;
    for (curva = malha.interface.curvaInicial; curva != NULL; curva = curva->proxCurva) {
        ponto = curva->pontoInicial;
        do {
            if (ponto->x > ponto_xmax->x) {
                ponto_xmax = ponto;
            } 
            ponto = ponto->prox;
        } while(ponto != curva->pontoInicial);
    }

    EncontraCelulaDaParticula(&malha, ponto_xmax->x, ponto_xmax->y, i, j);

    return;
}

void pontoMeio(MALHA malha, int *i, int *j) {
    CURVA *curva;
    PONTO *ponto;
    float xmin = malha.xMax, xmax = malha.xMin;
    float ymin = malha.yMax, ymax = malha.yMin;

    for (curva = malha.interface.curvaInicial; curva != NULL; curva = curva->proxCurva) {
        ponto = curva->pontoInicial;
        do {
            if (ponto->x < xmin) {
                xmin = ponto->x;
            } else if (ponto->x > xmax) {
                xmax = ponto->x;
            }
            if (ponto->y < ymin) {
                ymin = ponto->y;
            } else if (ponto->y > ymax) {
                ymax = ponto->y;
            }
            ponto = ponto->prox;
        } while(ponto != curva->pontoInicial);
    }

    EncontraCelulaDaParticula(&malha, (xmax - xmin) / 2, (ymax - ymin) / 2, i, j);

    return;
}
