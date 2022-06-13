.SUFFIXES:
.SUFFIXES: .o .f90 .f

FC := gfortran
FFLAGS := -O3 -g

ALGLIB := $(shell if [ -e $(CURDIR)/libalgencan.a ]; then echo true; fi)

ifneq ($(ALGLIB),true)

all: algencan MOPsolver

algencan: 
	$(MAKE) -C $(CURDIR)/algencan-3.1.1 
	mv -f $(CURDIR)/algencan-3.1.1/lib/libalgencan.a $(CURDIR)
endif

OBJECTS = globals.o lapack.o myproblem.o scalefactor.o evalfuns.o checkdF.o innersolver.o armijo.o MOPsolverNewton.o MOPsolverNewtonWoSafeg.o MOPsolverNG.o MOPsolverSD.o MOPsolverNewtonScalar.o  main.o

MOPsolver: $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) -L$(CURDIR) -lalgencan

globals.o: globals.f90
	$(FC) -c $(FFLAGS) globals.f90

globals.mod: globals.o

lapack.o: lapack.f
	$(FC) -c $(FFLAGS) lapack.f

myproblem.o: globals.mod myproblem.f90
	$(FC) -c $(FFLAGS) myproblem.f90

myproblem.mod: myproblem.o

scalefactor.o: globals.mod myproblem.mod scalefactor.f90
	$(FC) -c $(FFLAGS) scalefactor.f90

evalfuns.o: globals.mod myproblem.mod evalfuns.f90
	$(FC) -c $(FFLAGS) evalfuns.f90

checkdF.o: myproblem.mod checkdF.f90
	$(FC) -c $(FFLAGS) checkdF.f90
	
innersolver.o: globals.mod myproblem.mod innersolver.f90
	$(FC) -c $(FFLAGS) innersolver.f90

armijo.o: armijo.f90
	$(FC) -c $(FFLAGS) armijo.f90

MOPsolverNewton.o: globals.mod MOPsolverNewton.f90
	$(FC) -c $(FFLAGS) MOPsolverNewton.f90

MOPsolverNewtonWoSafeg.o: globals.mod MOPsolverNewtonWoSafeg.f90
	$(FC) -c $(FFLAGS) MOPsolverNewtonWoSafeg.f90

MOPsolverNG.o: globals.mod MOPsolverNG.f90
	$(FC) -c $(FFLAGS) MOPsolverNG.f90

MOPsolverSD.o: globals.mod MOPsolverSD.f90
	$(FC) -c $(FFLAGS) MOPsolverSD.f90

MOPsolverNewtonScalar.o: globals.mod MOPsolverNewtonScalar.f90
	$(FC) -c $(FFLAGS) MOPsolverNewtonScalar.f90

main.o: globals.mod myproblem.mod main.f90
	$(FC) -c $(FFLAGS) main.f90

CLEAN:
	rm -f *.mod *.o MOPsolver

.PHONY: CLEAN
