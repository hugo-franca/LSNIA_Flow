static char help[] = "Simulador Navier-Stokes. Projeto de mestrado, Hugo Franca\n\n";

#include <petscksp.h>

#include "FuncoesAuxiliares.h"
#include "Malha.h"
#include "InterfaceFrontTracking.h"
#include "MetodoProjecao.h"
#include "Visualizacao.h"

//Funcoes criadas para o estudo do problema da gota. 12/05/2025
#include "DropFunctions.h"

char nomePastaArquivos[300] = "";

const int velocidadeNewtoniana = 1;

double **debugsCell = NULL;


void LogOutputs(MALHA *M)
{
    static FILE *arq = NULL;

    if( M->interface.curvaInicial==NULL )
        return;

    if( arq==NULL ) {
        char nomeArq[500];
        sprintf(nomeArq, "ArquivosPlot/VTK/%s/logfile.txt", M->nomeArquivo);
        arq = fopen(nomeArq, "wt");
    }

    if( M->passoTemporal%M->intervaloPrint )
        return;

    // Assuming just one curve
    double droplet_radius = 0.0;
    LOOP_CURVA(M->interface.curvaInicial) {
        PONTO *p = loop_stuff.ponto;

        if( p->x>droplet_radius )
            droplet_radius = p->x;
    }


//    fprintf(arq, "%lf %.10lf %.10lf %.10lf\n", M->passoTemporal*M->dt, x_max - x_min, y_min, y_max);
    fprintf(arq, "%lf %.10lf\n", M->passoTemporal*M->dt, droplet_radius);
    fflush(arq);
    return;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
    PetscErrorCode ierr;
    SOLVER_PROJECAO solver;
    MALHA malha;
    LISTA_CELULAS_SURFACE listaCelulasSurf;
    double residuo = 0.0;
    int n, rank, tam;
    float *areaFluido, *spread, *altura, energiaAnt = 0.0;
    float *energiaCinetica, *energiaGravitacional, *energiaSuperficial, *energiaDissipativa;
    
    PetscInitialize(&argc, &args, (char*)0, help);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    listaCelulasSurf.prim = NULL;
    listaCelulasSurf.ult = NULL;


    /// ==== Pegando o nome do arquivo na linha de comando e lendo o arquivo
    char arquivoSim[200];
    PetscBool flagArquivo;
    ierr = PetscOptionsGetString(NULL, NULL, "-arquivo", arquivoSim, 200, &flagArquivo); CHKERRQ(ierr);
    if( !flagArquivo )
        DeuProblema("\n\nERRO: Forneca o nome do arquivo de entrada.\n\n");

    //printf("Vamos ler o modelo do arquivo:\n");
    malha = LerModeloDoArquivo(arquivoSim);
    //printf("Arquivo já foi lido:\n");


    // Cria um arquivo para armazenar algumas informações (spread, área, energias)
    FILE *arqInfo;
    char nomeArqInfo[250];

    nomeArqInfo[0] = '\0';
    sprintf(nomeArqInfo, "%s_INFO", arquivoSim);

    arqInfo = fopen(nomeArqInfo, "w");


    PrintDebug("Inicializando estruturas\n");
    strcpy(nomePastaArquivos, malha.nomeArquivo);

    /// ==== Aloca memoria para todas matrizes de solucoes
    ierr = InicializaEstruturasProjecao(&solver, malha); CHKERRQ(ierr);

    n = -1;
    // CarregarEstadoBinario(&malha, solver.UVelho, solver.VVelho, solver.PVelho, solver.TxxVelho, solver.TxyVelho, solver.TyyVelho, solver.lambdaVelho, solver.muVelho, solver.nuVelho, &listaCelulasSurf, &n);
    //CarregarEstado(&malha, solver.UVelho, solver.VVelho, solver.PVelho, solver.TxxVelho, solver.TxyVelho, solver.TyyVelho, &listaCelulasSurf, &n);

    /// ==== Imprimeces a solucao inicial e a malha computacional inicial nos arquivos vtk
//    ImprimeArquivoVTK(solver.UVelho, solver.VVelho, solver.WVelho, solver.PVelho,
//                              solver.TxxVelho, solver.TxyVelho, solver.TyyVelho, solver.TxtVelho, solver.TytVelho, solver.TttVelho,
//                              solver.lambdaVelho, solver.muVelho, solver.nuVelho, malha, n+1);;
    ImprimeInterfaceVTK(malha, n+1);
    DesenhaMalhaVTK(malha, 0);
//    CriaMalhaQuadTree(&malha);
//    DesenhaMalhaQuadVTK(malha, 0);
//    DeuProblema("desenhou\n");

    // ImprimeArquivoBinario(solver.UVelho, solver.VVelho, malha, n+1);

    // Estrutura para armazenar a área ocupada pelo fluido em alguns passos de tempo
    tam = (int) (malha.Nt / malha.intervaloPrint);
    areaFluido = (float *) calloc(tam, sizeof(float));
    spread = (float *) calloc(tam, sizeof(float));
    altura = (float *) calloc(tam, sizeof(float));
    energiaCinetica = (float *) calloc(tam, sizeof(float));
    energiaGravitacional = (float *) calloc(tam, sizeof(float));
    energiaSuperficial = (float *) calloc(tam, sizeof(float));
    energiaDissipativa = (float *) calloc(tam, sizeof(float));

    double t = 0.0;

    PrintDebug("Inicializando passos temporais \n");
    for( n++; n<malha.Nt; n++ ) {
        malha.passoTemporal = n;

        clock_t start, end;
        double cpu_time_used;
        start = clock();

        LogOutputs(&malha);
        // CalculaVorticeContracao(&malha, solver.UVelho);
        

        ierr = ExecutaPassoTemporalProjecao(&solver, &malha, n, rank, &listaCelulasSurf, &residuo); CHKERRQ(ierr);

        CalculaTodasEnergias(&malha, solver.UNovo, solver.VNovo, solver.PNovo, n);
        // CalculatePenetrationSpreading(&malha, solver.VVelho);
        // CalculateBridgeHeight(&malha);

//        CalculaTodasEnergias(&malha, solver.UNovo, solver.VNovo, solver.PNovo, n);
    //    DeuProblema("para\n");

        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;	

	    energiaAnt = CalculaEnergiaDissipativaVisc(&malha, solver.UNovo, solver.VNovo, energiaAnt);

        if( n%malha.intervaloPrint==0 && rank==0 ) {
            PrintDebug("Passo %d finalizado. residuo=%.15e CPU_TIME: %g\n\n", n+1, residuo, cpu_time_used);
            areaFluido[n / malha.intervaloPrint] = calcAreaFluido(malha);
            spread[n / malha.intervaloPrint] = calcSpread(malha);
            altura[n / malha.intervaloPrint] = calcAltura(malha);
            energiaCinetica[n / malha.intervaloPrint] = CalculaEnergiaCinetica(&malha, solver.UNovo, solver.VNovo);
            energiaGravitacional[n / malha.intervaloPrint] = CalculaEnergiaGrav(&malha);
            energiaSuperficial[n / malha.intervaloPrint] = CalculaEnergiaSuperficie(&malha);
            energiaDissipativa[n / malha.intervaloPrint] = energiaAnt;

	}

        if( isnan(residuo) || isinf(residuo) )
            DeuProblema("\n\n Valor estranho no residuo: %lf\n\n", residuo);


        /// ==== Imprimindo o arquivo VTK que contem as propriedades no dominio (U, V, P, T, etc...)
        if( (malha.intervaloVTK_prop) && ((n+1)%(malha.intervaloVTK_prop)==0) && rank==0 ) {

             ImprimeArquivoVTK(solver.UVelho, solver.VVelho, solver.WVelho, solver.PVelho,
                              solver.TxxVelho, solver.TxyVelho, solver.TyyVelho, solver.TxtVelho, solver.TytVelho, solver.TttVelho,
                              solver.lambdaVelho, solver.muVelho, solver.nuVelho,
                              solver.evptStructureVelho, solver.norm_tau_dev,
                              solver.MuNewtGenVelho,
                              solver.termoDominanteXX, solver.termoDominanteXY, solver.termoDominanteYY,
                              solver.residualFaceU, solver.residualFaceV, solver.residualCell,
                              malha, n+1);
            DesenhaMalhaVTK(malha, 0);
	  
	    //Diametro da gota: essa funcao deve ser comentada para outras geometrias. 12/05/2025.
	    //DropDiameterAxisymmetric(malha, t);
	    //EstudoCodigo(malha, listaCelulasSurf, n);
            // ImprimeArquivoBinario(solver.UVelho, solver.VVelho, malha, n+1);
        }


        /// ==== Imprimindo o arquivo VTK que contem a superficie livre (se houver)
        if( (malha.intervaloVTK_surf) && ((n+1)%(malha.intervaloVTK_surf)==0) && rank==0 ) {
            ImprimeInterfaceVTK(malha, n+1);
            DesenhaMalhaVTK(malha, 0);
//            getchar();
        }

        /// ==== Imprimindo arquivo de Backup usado pra debugar problemas...
        if( (malha.intervaloEstado) &&((n+1)%(malha.intervaloEstado)==0) && rank==0 )
            SalvarEstadoBinario(malha, solver.UVelho, solver.VVelho, solver.PVelho, solver.TxxVelho, solver.TxyVelho, solver.TyyVelho, solver.lambdaVelho, solver.muVelho, solver.nuVelho, listaCelulasSurf, n);


//        if( n==54500 ) {
//            SalvarEstadoBinario(malha, solver.UVelho, solver.VVelho, solver.PVelho, solver.TxxVelho, solver.TxyVelho, solver.TyyVelho, solver.lambdaVelho, solver.muVelho, solver.nuVelho, listaCelulasSurf, n);
//            DeuProblema("parou\n\n");
//        }

        t = t + malha.dt;

        // DeuProblema("parou\n");
    }

    fprintf(arqInfo, "\n\nÁrea do fluido calculada em alguns passos de tempo:\n");
    fprintfVetorFormatado(arqInfo, areaFluido, tam);

    fprintf(arqInfo, "\n\nSpread máximo calculado em alguns passos de tempo:\n");
    fprintfVetorFormatado(arqInfo, spread, tam);

    fprintf(arqInfo, "\n\nAltura máxima calculada em alguns passos de tempo:\n");
    fprintfVetorFormatado(arqInfo, altura, tam);

    fprintf(arqInfo, "\n\nEnergia cinética calculada em alguns passos de tempo:\n");
    fprintfVetorFormatado(arqInfo, energiaCinetica, tam);

    fprintf(arqInfo, "\n\nEnergia gravitacional calculada em alguns passos de tempo:\n");
    fprintfVetorFormatado(arqInfo, energiaGravitacional, tam);

    fprintf(arqInfo, "\n\nEnergia superficial calculada em alguns passos de tempo:\n");
    fprintfVetorFormatado(arqInfo, energiaSuperficial, tam);

    fprintf(arqInfo, "\n\nEnergia dissipativa calculada em alguns passos de tempo:\n");
    fprintfVetorFormatado(arqInfo, energiaDissipativa, tam);

    /// === Imprimindo o arquivo do ultimo passo temporal
    if( rank==0 ) {
        //ImprimeArquivoVTK(solver.UVelho, solver.VVelho, solver.WVelho, solver.PVelho, solver.TxxVelho, solver.TxyVelho, solver.TyyVelho, malha, n);
        //ImprimeArquivoVTK(solver.UVelho, solver.VVelho, solver.WVelho, solver.PVelho, solver.TxxVelho, solver.TxyVelho, solver.TyyVelho, solver.lambdaVelho, solver.muVelho, solver.nuVelho, malha, n);
        //ImprimeInterfaceVTK(malha, n);
    }


    FinalizaPrograma(&malha, &listaCelulasSurf, &solver);
    ierr = PetscFinalize();
    FinalizouNormalmente("\n\nAcabou Normalmente!\n\n");
    return 0;
}
