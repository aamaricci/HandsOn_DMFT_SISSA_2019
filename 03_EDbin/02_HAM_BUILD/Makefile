FC=gfortran

EXE1 = 01_build_ham
EXE2 = 02_diag_test

DIREXE = ./test

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

OBJS     = M1_COMMON_VARS.o M2_SPARSE_MATRIX.o M3_SETUP.o M4_BUILD_HAM.o  

FFLAG = -O3 -ffree-line-length-none
#FFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none


ARGS=-llapack -lblas


all: $(OBJS)
	${FC} ${FFLAG} ${OBJS} ${EXE1}.f90 -o ${DIREXE}/${EXE1} ${ARGS}
	${FC} ${FFLAG} ${OBJS} ${EXE2}.f90 -o ${DIREXE}/${EXE2} ${ARGS}
	@echo ""


objs: ${OBJS}

clean: 
	@echo "Cleaning:"
	@rm -fv *.mod
	@rm -fv *.o
	@rm -fv *~
	@rm -fv ${DIREXE}/${EXE1}
	@rm -fv ${DIREXE}/${EXE2} 
	@echo ""

.f90.o:	
	${FC} ${FFLAG} -c $< 


