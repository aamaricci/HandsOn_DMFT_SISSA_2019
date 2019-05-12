FC=gfortran

EXE1 = 01_build_ham
EXE2 = 02_lanc_test
EXE3 = 03_arpack_test

DIREXE = ./test

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

OBJS     = M1_COMMON_VARS.o M2_SPARSE_MATRIX.o M3_LANCZOS.o M4_SETUP.o M5_MATVEC_PRODUCT.o  

FFLAG = -O3 -ffree-line-length-none
#FFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none


ARGS=-llapack -lblas
AR_INCS=$(shell pkg-config --cflags arpack)
AR_ARGS=$(shell pkg-config --libs arpack)


all: $(OBJS)
	${FC} ${FFLAG} ${OBJS} ${EXE1}.f90 -o ${DIREXE}/${EXE1} ${ARGS}
	${FC} ${FFLAG} ${OBJS} ${EXE2}.f90 -o ${DIREXE}/${EXE2} ${ARGS}
	${FC} ${FFLAG} ${OBJS} ${EXE3}.f90 -o ${DIREXE}/${EXE3} ${AR_INCS} ${AR_ARGS} ${ARGS} 
	@echo ""


objs: ${OBJS}

clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@rm -f ${DIREXE}/${EXE1}
	@rm -f ${DIREXE}/${EXE2} 
	@echo ""

.f90.o:	
	${FC} ${FFLAG} -c $< 


