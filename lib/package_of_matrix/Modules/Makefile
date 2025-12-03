FC= gfortran
#FC= mpxlf90
#FLAGS= -c  -q64  -O4 
FLAGS= -c -O1 -ffree-line-length-none
SUFFIX=
LF= 
LIB=modules_90.a

all: 
	(make -f Compile FC="$(FC)" FLAGS="$(FLAGS)"  LIB="$(LIB)" )

clean:
	(make -f Compile  clean ) ;\
 	rm *.mod
