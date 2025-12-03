#FC= xlf90
FC=mpif90 
FLAGS= -c   -O4 
SUFFIX= -qsuffix=f=f90
LF= 
LIB=modules_MPI_90.a

all: 
	(cp ../Modules/*.mod .) ;\
        (cp ../Modules/maxent_stoch.G90 . ) ;\
	(make -f PRE1 ) ;\
	(make -f Compile FC="$(FC)" FLAGS="$(FLAGS)"  LIB="$(LIB)" )

clean:
	(make -f PRE1 clean ) ;\
	(rm maxent_stoch.G90 ) ;\
	rm *.mod ;\
	(make -f Compile clean )
