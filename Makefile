# Warnings
WFLAGS	:= -Wall -Wextra

# Optimization and architecture
OPT		:= -O3

# Language standard
CCSTD	:= -std=c99
CXXSTD	:= -std=c++11


# Linker options
LDOPT 	:= $(OPT)
LDFLAGS := 
BIN = "/usr/local/gcc/6.4.0/bin/gcc"
.DEFAULT_GOAL := all

.PHONY: debug
debug : OPT  := -O0 -g -fno-omit-frame-pointer -fsanitize=address
debug : LDFLAGS := -fsanitize=address
debug : ARCH :=
debug : $(EXEC)

all : main mainJobA helperJobB

main: main.cu
	module load cuda;nvcc -o main $(OPT) main.cu -fopenmp -ccbin $(BIN)
	
mainJobA: mainJobA_godunov.cu
	module load cuda;nvcc -o mainJobA_godunov $(OPT) mainJobA_godunov.cu -fopenmp -ccbin $(BIN)

helperJobB: helperJobB_kernel.cu
	module load cuda;nvcc -o helperJobB_kernel $(OPT) helperJobB_kernel.cu -fopenmp -ccbin $(BIN)

.PHONY: clean
clean:
	@ rm -f $(EXEC) *.o
