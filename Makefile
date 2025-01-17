#================================================================
#   Copyright (C) 2021 Sangfor Ltd. All rights reserved.
#   
#   File Name：Makefile
#   Author: Wenqiang Wang
#   Created Time:2021-10-30
#   Discription:
#
#================================================================

GPU_CUDA := #ON
FLOAT16 := ON
SVE := ON

XFAST := ON
ZFAST := #ON

FREE_SURFACE := ON
PML := ON
SOLVE_DISPLACEMENT := ON
Terrain_Smooth := #ON 
DealWithFirstLayer := ON

LayeredStructureTerrain := ON
StructureTerrain := ON


SET_BASIN := #ON

SRCDIR := ./SRC-fp16-sve


CCHOME   := /opt/HPCKit/24.6.30/compiler/bisheng
CUDAHOME := /workspace/public/software/tool/cuda/cuda-11.6
MPIHOME  := /opt/HPCKit/24.6.30/hmpi/bisheng/hmpi
PROJHOME := /home/Shenchao/wangwq/software/proj-8.1.0


CC := clang

#General Compiler
ifdef GPU_CUDA
GC := $(CUDAHOME)/bin/nvcc -rdc=true -maxrregcount=64 -arch=sm_70 #-Xptxas=-v 
LIBS := -L$(CUDAHOME)/lib64 -lcudart -lcublas
INCS := -I$(CUDAHOME)/include 
else
GC := clang++
endif



LIBS += -L$(MPIHOME)/lib -lmpi
INCS += -I$(MPIHOME)/include 


LIBS += -L$(PROJHOME)/lib -lproj
INCS += -I$(PROJHOME)/include  

LIBS += -lunwind

LIBS += --rtlib=compiler-rt
#LDLIBS = --rtlib=compiler-rt


OBJDIR := ./obj
BINDIR := ./bin


CFLAGS := -c -O2 -Wno-all

CFLAGS += -mcpu=native -march=armv8-a+sve
wave_deriv_file := fp16_wave_deriv_sve.o
wave_rk_file := wave_rk.o

#LFLAGS := -fopenmp -O2
LFLAGS := -O2

GCFLAGS := 

ifdef GPU_CUDA
#LFLAGS += -Xptxas=-v 

#LFLAGS += -arch=sm_70 -rdc=true -Xptxas=-v 
#GCFLAGS += --fmad=false 
GCFLAGS += -x cu
endif

vpath

vpath % $(SRCDIR)
vpath % $(OBJDIR)
vpath % $(BINDIR)


DFLAGS_LIST := XFAST ZFAST GPU_CUDA FLOAT16 FREE_SURFACE PML SOLVE_DISPLACEMENT \
			   Terrain_Smooth DealWithFirstLayer SET_BASIN LayeredStructureTerrain StructureTerrain

DFLAGS := $(foreach flag,$(DFLAGS_LIST),$(if $($(flag)),-D$(flag)))


OBJS := cjson.o init_gpu.o init_grid.o init_MPI.o main.o getParams.o create_dir.o \
		run.o printInfo.o modelChecking.o  cpu_Malloc.o MPI_send_recv.o data_io.o \
		coord.o terrain.o  medium.o dealMedium.o crustMedium.o calc_CFL.o CJM.o \
		contravariant.o MPI_send_recv.o MPI_send_recv_FLOAT.o multiSource.o $(wave_deriv_file) $(wave_rk_file) \
		propagate.o freeSurface.o singleSource.o station.o PGV.o addMoment.o \
		init_pml_para.o pml_deriv.o pml_rk.o pml_freeSurface.o

OBJS := $(addprefix $(OBJDIR)/,$(OBJS))


$(BINDIR)/main: $(OBJS)
	$(GC) $(LFLAGS) $(LIBS) $^ -o $@ 


$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(GC) $(CFLAGS) $(DFLAGS) $(GCFLAGS) $(INCS)  $^ -o $@


$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	-rm $(OBJDIR)/* -rf
	-rm $(BINDIR)/* -rf
	-rm output -rf

