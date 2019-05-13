FC=gfortran

EXE1 = 01_lanc_gf
EXE2 = 02_arpack_gf

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

OBJS     = M1_COMMON_VARS.o M2_SPARSE_MATRIX.o M3_LANCZOS.o M4_SETUP.o M5_MATVEC_PRODUCT.o M6_GF_NORMAL.o

FFLAG = -O3 -ffast-math -march=native -funroll-all-loops
#FFLAG = -O3 -ffree-line-length-none
#FFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none


ARGS=-llapack -lblas
AR_INCS=$(shell pkg-config --cflags arpack)
AR_ARGS=$(shell pkg-config --libs arpack)

all: $(OBJS)
	${FC} ${FFLAG} ${OBJS} ${EXE1}.f90 -o ${EXE1} ${ARGS}
	${FC} ${FFLAG} ${OBJS} ${EXE2}.f90 -o ${EXE2} ${AR_INCS} ${AR_ARGS} ${ARGS} 
	@echo ""


objs: ${OBJS}

clean: 
	@echo "Cleaning:"
	@rm -fv *.mod
	@rm -fv *.o
	@rm -fv *~
	@rm -fv ${EXE1}
	@rm -fv ${EXE2} 
	@echo ""

.f90.o:	
	${FC} ${FFLAG} -c $< 


