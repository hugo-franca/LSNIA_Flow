#include "FuncoesAuxiliares.h"

extern char nomePastaArquivos[];

clock_t start_timer, end_timer;



void PushPilha(PILHA_GERAL *P, void *Valor)
{
    ITEM_PILHA_GERAL *novoItem;

    novoItem = (ITEM_PILHA_GERAL *)malloc(sizeof(ITEM_PILHA_GERAL));
    novoItem->valor = Valor;
    novoItem->prox = P->prim;
    P->prim = novoItem;
    return;
}

void *PopPilha(PILHA_GERAL *P)
{
    ITEM_PILHA_GERAL *item;
    void *valor;

    item = P->prim;
    P->prim = item->prox;
    valor = item->valor;
    free(item);
    return valor;
}

void InsereParametroOpcional(PILHA_GERAL *P, char *StringParametro)
{
    PARAMETRO_OPCIONAL *novo;
    int i, j = 0;

    novo = (PARAMETRO_OPCIONAL *)malloc( sizeof(PARAMETRO_OPCIONAL) );

    // Separando o nome do parametro
    for( i=0; StringParametro[i]!='\0' && StringParametro[i]!='='; i++ )
        novo->nomeParametro[i] = StringParametro[i];
    novo->nomeParametro[i] = '\0';
    if( StringParametro[i]!='=' )
        DeuProblema("PROBLEMA NA DEFINICAO EM ALGUM DOS PARAMETROS OPCIONAIS\n");

    // Separando o valor do parametro
    char stringValor[300];
    for( i++; StringParametro[i]!='\0'; i++ ) {
        stringValor[j++] = StringParametro[i];
    }
    stringValor[j] = '\0';
    sscanf(stringValor, "%lf", &(novo->valorParametro));

    PushPilha(P, novo);
    return;
}

void **AlocaMatriz(int Linhas, int Colunas, size_t TamanhoElemento, void *ValorInicial)
{
    void **matriz;
    int i, j, byteAtualValor, byteAtualMatriz;
    unsigned char *charValorInicial, *charMatriz;

    charValorInicial = (unsigned char *)ValorInicial;

    matriz = (void **)malloc( Linhas*sizeof(void *) );
    for( i=0; i<Linhas; i++ )
        matriz[i] = malloc( Colunas*TamanhoElemento );

    //Inicializando todos os elementos com o valor inicial (soh funciona se for um tipo int)
    if( ValorInicial!=NULL ) {
        for(i=0; i<Linhas; i++) {
            charMatriz = (unsigned char *)matriz[i];
            byteAtualMatriz = 0;
            for(j=0; j<Colunas; j++) {
                for( byteAtualValor = 0; byteAtualValor<TamanhoElemento; byteAtualValor++ )
                    charMatriz[byteAtualMatriz++] = charValorInicial[byteAtualValor];
            }
        }
    }


    return matriz;
}

void DesalocaMatriz(void **Matriz, int Linhas, int Colunas)
{
    int i;

    if( Matriz==NULL )
        return;

    for( i=0; i<Linhas; i++ )
        free(Matriz[i]);
    free(Matriz);
}

void TrocaPonteirosMatrizes(double ***Matriz1, double ***Matriz2)
{
    double **Temp;

    Temp = *Matriz1;
    *Matriz1 = *Matriz2;
    *Matriz2 = Temp;

    return;
}

void InicializaSolverKSP(KSP *Ksp, Mat A)
{
    //PC pc;

    KSPCreate(PETSC_COMM_WORLD, Ksp);
    KSPSetOperators(*Ksp, A, A);
    //KSPGetPC(*Ksp, &pc);
    //PCSetType(pc, PCNONE);
    KSPSetTolerances(*Ksp, 1.e-9, 1.e-9, PETSC_DEFAULT, 50000);
    KSPSetInitialGuessNonzero(*Ksp, PETSC_TRUE); //Nao vou usar ZERO como chute inicial. Vou usar a solucao do passo anterior
    KSPSetFromOptions(*Ksp);
    KSPSetUp(*Ksp);
    return;
}

void ImprimeMatrizPraMatlab(Mat M, const char *NomeArquivo)
{
    PetscViewer viewer;
    MPI_Comm comm;

    PetscObjectSetName((PetscObject)M, "matriz");
    PetscObjectGetComm((PetscObject)M, &comm);
    PetscViewerCreate(comm, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERASCII);
    PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    PetscViewerFileSetName(viewer, NomeArquivo);
    MatView(M, viewer);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);

    printf("Imprimiu matriz no arquivo\n");
    return;
}

void ImprimeVetorPraMatlab(Vec V, const char *NomeArquivo)
{
    PetscViewer viewer;
    MPI_Comm comm;

    PetscObjectSetName((PetscObject)V, "vetor");
    PetscObjectGetComm((PetscObject)V, &comm);
    PetscViewerCreate(comm, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERASCII);
    PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    PetscViewerFileSetName(viewer, NomeArquivo);
    VecView(V, viewer);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);

    printf("Imprimiu vetor no arquivo\n");
    return;
}

void BroadcastMatrizMPI(double **Matriz, int Linhas, int Colunas)
{
//    int i;
//
//    ///Nao vou usar essa funcao no momento. Por enquanto vou rodar apenas sequencial
//    ///Essa funcao eh apenas pra paralelizacao
//    return;
//
//    for( i=0; i<Linhas; i++ ) {
//        MPI_Bcast(Matriz[i], Colunas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    }
//
//    return;
}

double DeuProblema(const char* format, ...)
{
    FILE *arq;
    char nomeArq[300];
    va_list argptr;

    va_start(argptr, format);
    vprintf(format, argptr);
    va_end(argptr);

    sprintf(nomeArq, "ArquivosPlot/VTK/%s/DeuProblema.txt", nomePastaArquivos);
    arq = fopen(nomeArq, "wt");
    va_start(argptr, format);
    vfprintf(arq, format, argptr);
    va_end(argptr);
    fclose(arq);

    exit(0);
    return -100000;
}

double FinalizouNormalmente(const char* format, ...)
{
    FILE *arq;
    char nomeArq[300];
    va_list argptr;

    va_start(argptr, format);
    vprintf(format, argptr);
    va_end(argptr);

    sprintf(nomeArq, "ArquivosPlot/VTK/%s/FinalizouOK.txt", nomePastaArquivos);
    arq = fopen(nomeArq, "wt");
    va_start(argptr, format);
    vfprintf(arq, format, argptr);
    va_end(argptr);
    fclose(arq);

    exit(0);
    return -100000;
}

void PrintDebug(const char* format, ...)
{
    static FILE *arq = NULL;

    if( (arq==NULL) && (strlen(nomePastaArquivos)>0) ) {
        char nomeArq[500];
        sprintf(nomeArq, "ArquivosPlot/VTK/%s/print_output.txt", nomePastaArquivos);
        arq = fopen(nomeArq, "wt");
    }

    //Imprimindo no prompt
    va_list argptr;
    va_start(argptr, format);
    vprintf(format, argptr);
    va_end(argptr);
    fflush(stdout);

    //Imprimindo no arquivo
    if( arq ) {
        va_start(argptr, format);
        vfprintf(arq, format, argptr);
        va_end(argptr);
        fflush(arq);
    }

    return;
}


void IniciaTimer()
{
    start_timer = clock();
}

double FinalizaTimer()
{
    end_timer = clock();
    return ((double) (end_timer - start_timer)) / CLOCKS_PER_SEC;
}

void Check_Read(int QtdLida, int QtdEsperada)
{
    if( QtdLida==QtdEsperada )
        return;

    DeuProblema("Problema em algum fread ou fscanf.\n\n");
    return;
}


