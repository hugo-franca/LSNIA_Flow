# CFLAGS	  = -I"/usr/include/python3.8/" -lpython3.8 -v
CFLAGS	  = -I"/usr/include/python3.10/" -lpython3.10 -v
FFLAGS    =
CPPFLAGS  =
FPPFLAGS  =
LOCDIR    = src/ksp/ksp/examples/tests/

# Se o usuario nao tiver fornecido o nome do executavel, vamos usar o nome padrao "main"
ifndef NOME_EXEC
	NOME_EXEC = main
endif

#PETSC_DIR=/home/diego/petsc-3.18.5
#PETSC_DIR=/home/irineu/Documents/FCT/Pesquisa/Codes/petsc-3.16.2
# PETSC_DIR=/home/irineu/OpenFOAM/irineu-9/ThirdParty/petsc-3.16.5
# PETSC_DIR=/home/irineu/OpenFOAM/irineu-9/ThirdParty/petsc-3.16.5
PETSC_DIR=/mnt/d/share/basilisk_codes/petsc/petsc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

main: ContornoSuperficie.o Derivadas.o EquacaoConstitutiva.o EVPT_Equacoes.o EquacaoMovimento.o EquacaoPoisson.o FuncoesAuxiliares.o InterfaceFrontTracking.o main.o Malha.o MetodoProjecao.o Visualizacao.o DropFunctions.o
	-${CLINKER} $^ -o ${NOME_EXEC} ${PETSC_KSP_LIB} ${CFLAGS}
