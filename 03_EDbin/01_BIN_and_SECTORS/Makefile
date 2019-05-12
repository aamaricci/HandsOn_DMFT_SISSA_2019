FC=gfortran

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

OBJS     = COMMON_VARS.o AUX_FUNX.o


all: $(OBJS)
	${FC} ${OBJS} 01_bit_operations.f90 -o 01_bit_operations ${INCS} ${ARGS}
	${FC} ${OBJS} 02_build_sectors.f90 -o 02_build_sectors ${INCS} ${ARGS}
	@echo ""


clean: 
	@echo "Cleaning:"
	@rm -fv *.mod *.o *~ 01_bit_operations 02_build_sectors
	@echo ""

.f90.o:	
	$(FC) $(FFLAG) -c $< $(INCS)
