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

main: main.cpp
	mpicxx $(CXXSTD) $(WFLAGS) $(OPT) -fopenmp -o $@ $<
	
mainJobA: mainJobA_godunov.cu
	$(CXXSTD) $(WFLAGS) $(OPT) -fopenmp -o $@ $<

helperJobB: helperJobB_kernel.cu
	$(CXXSTD) $(WFLAGS) $(OPT) -fopenmp -o $@ $<

.PHONY: clean
clean:
	@ rm -f $(EXEC) *.o
