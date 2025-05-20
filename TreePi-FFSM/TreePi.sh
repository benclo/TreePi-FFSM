#!/bin/bash
#SBATCH -p main      # Queue name
#SBATCH --mem=200000

# Compiler and flags
GCC="g++ -O3"

# List all source files
SRC_FILES="adj_matrix.cpp embedding.cpp TreePi.cpp BTreePlus.cpp Graph.cpp graphFFSM.cpp p1.cpp p1_helper.cpp runner.cpp"

# List all object files
OBJ_FILES="adj_matrix.o embedding.o TreePi.o BTreePlus.o Graph.o graphFFSM.o p1.o p1_helper.o runner.o"

# Compile individual object files
for src in $SRC_FILES; do
    obj="${src%.cpp}.o"
    $GCC -c $src -o $obj
done

$GCC -o TreePi $OBJ_FILES

./TreePi
