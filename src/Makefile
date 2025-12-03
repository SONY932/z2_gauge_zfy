FC=mpiifort -fc=ifort 

# 默认编译选项（优化模式）
FLAGS_OPT = -c -O3
# 调试编译选项（debug模式）
FLAGS_DEBUG = -c -O0 -g -traceback -check all -fpe0 -warn all -debug extended

# 根据 MODE 变量选择编译模式（默认为优化模式）
# 使用方法：make           （优化模式）
#          make MODE=debug （调试模式）
ifeq ($(MODE),debug)
    FLAGS = $(FLAGS_DEBUG)
else
    FLAGS = $(FLAGS_OPT)
endif

LF = -mkl
HOME = /home/lijiaqiang/lib/package_of_matrix
LIBS= $(HOME)/Modules/modules_90.a \
      $(HOME)/MyEis/libeis.a \
      $(HOME)/MyNag/libnag.a \
      $(HOME)/MyLin/liblin.a \
      $(HOME)/Ran/libran.a \
#      $(HOME)/LaPack/lapack.a \
#      $(HOME)/Blas/libblas.a

all:
	@echo "编译模式: $(if $(filter debug,$(MODE)),调试模式,优化模式)"
	cp $(HOME)/Modules/*.mod . ;\
	(make -f Compile  FC="$(FC)" LF="$(LF)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)") 

debug:
	@echo "强制使用调试模式编译..."
	$(MAKE) MODE=debug

clean: 	
	(make -f Compile  clean );\
	rm -f *.mod

.PHONY: all debug clean

