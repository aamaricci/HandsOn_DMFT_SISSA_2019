#CODE TO BE COMPILED:
EXE  = ipt_hm_matsubara

FC   = gfortran
FFLAG = -O3 -funroll-all-loops -Wno-argument-mismatch
HERE = `pwd`
DIREXE = .


####### TEXT BELOW IS HARDLY TO BE CHANGED ########
OBJS   = IPT_FFTGF.o

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90
.PHONY: code 

#Compile the FFTpack library and then the code
all: code

#this compile the code, specified above in the EXE variable
code: 	$(OBJS)
	${FC} ${FFLAG} ${OBJS} ${EXE}.f90 -o ${DIREXE}/${EXE}
	@echo ""

.f90.o:	
	$(FC) $(FFLAG) -c $<

clean: 
	@echo "Cleaning:"
	@rm -fv *.mod
	@rm -fv *.o
	@rm -fv *~
	@rm -fv ${DIREXE}/${EXE}
	@echo ""
