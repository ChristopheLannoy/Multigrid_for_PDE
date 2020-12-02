# librairies de PRIMME
LIBP = -L./primme/ -lprimme

# librairies de JADAMILU_2.0
#LIBJ = -L./JADAMILU_2.0/lib/INT32GNU/ -ljadamilu

# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/ 

# librairies de SuiteSparse
L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
L3 = SuiteSparse/AMD/Lib/libamd.a 
L4 = SuiteSparse/CAMD/Lib/libcamd.a  
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a
# toutes les librairies
LIB = $(LIBP) $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack

COPT = -O3 -Wall

default: main

clean: 
	rm *.o 
	rm main
	rm *.txt

main: main.c prob.o time.o plot.o plot_vector.o plot_eigen.o umfpack.o  mat_vec.o diagonal.o norme.c smoothing.o norme_residu.o solve_L.o solve_U.o multi_grid.o elem_number.o coarse_to_fine.o fine_to_coarse.o primme.o
	gcc $(COPT) $^ -o $@ $(LIB)

umfpack.o: umfpack.c
	gcc $(COPT) -c $< -o $@ -ISuiteSparse/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include

primme.o: primme.c
	gcc $(COPT) -c $< -o $@ $(INCP)

%.o: %.c
	gcc $(COPT) -c $< -o $@ 


