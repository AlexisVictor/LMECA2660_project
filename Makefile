INC_DIR := -I/home/thoussiau/Projet_2660/lib_petsc/include
LIB_DIR := -L/home/thoussiau/Projet_2660/lib_petsc/lib -Wl,-rpath=/home/thoussiau/Projet_2660/lib_petsc/lib


LIB := -lpetsc

CXX_FLAGS := -O0 -Wall -g -Werror 

#Compilation
all :
	gcc -o project project.c poisson.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)
	./project -pc_type lu -ksp_type fgmres
#Delete of executable files
clean :
	rm project

#Delete of results
clean_txt :
	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*
