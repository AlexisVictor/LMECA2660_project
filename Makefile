INC_DIR := -I/opt/homebrew/Cellar/petsc/3.19.0/include#/usr/local/opt/petsc/include
LIB_DIR := -L/opt/homebrew/Cellar/petsc/3.19.0/lib#/usr/local/opt/petsc/lib

LIB := -lpetsc.3.19.0

LIB := -lpetsc

CXX_FLAGS := -O0 -Wall -g -Werror 

#Compilation
all :
	gcc -o project project.c poisson.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)
	./project -pc_type lu -ksp_type fgmres
#Delete of executable files
clean :
	rm projet

#Delete of results
clean_txt :
	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*