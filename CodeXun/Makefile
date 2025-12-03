FC=mpiifort -fc=ifort 
FLAGS= -c -O3
# -c -fast -xavx -O3 -opt-report
# -openmp -openmp-report -xavx -vec-report2 -c -O3
SUFFIX= 
# -check all -traceback
LF = -mkl
# -openmp 
#LF =  -q64  
HOME = ../../../Lib_90_new
LIBS= $(HOME)/Modules/modules_90.a \
      $(HOME)/MyEis/libeis.a \
      $(HOME)/MyNag/libnag.a \
      $(HOME)/MyLin/liblin.a \
      $(HOME)/Ran/libran.a \
#      $(HOME)/LaPack/lapack.a \
#      $(HOME)/Blas/libblas.a
all:
	cp $(HOME)/Modules/*.mod . ;\
	(make -f Compile  FC="$(FC)" `mpif90 -showme:compile`  LF="$(LF)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" SUFFIX="$(SUFFIX)") 
clean: 	
	(make -f Compile  clean );\
	rm *.mod
