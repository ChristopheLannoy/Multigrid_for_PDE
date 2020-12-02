#include <stdio.h>
#include <stdlib.h>

int plot_vector(int n, double *x)
/*
 * Plot le vecteur x (axe y) 
 * n : dimension du vecteur x
 * (axe x : 0, 1, 2, 3, 4, ... n)
 * 
 */
 
{
	FILE *data = fopen("data_vector.txt", "w"); //écriture de fichier
	
	int i;
	for(i = 0; i <= n; i++){
		fprintf(data, "%d %e\n", i, x[i]);
	}
	fclose(data);
	
	// fichier commande
	FILE *cmd = fopen("cmd_vect.txt", "w");
	fprintf(cmd,
		"unset key\n"
		//"set pm3d at b \n"
		"set logscale y\n"
		"set title 'Evolution de la norme residuelle pour %d iterations'\n" 
		"plot 'data_vector.txt' using 1:2 with lines\n",n
		 );
	fclose(cmd);

	//exécution
	system("gnuplot -persistent cmd_vect.txt");
	
	return 0;
} 
