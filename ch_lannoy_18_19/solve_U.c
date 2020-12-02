void solve_U(int n,int *ia, int *ja, double *a, double *r, double *temp)

/*
 * Solveur triangulaire supérieur
 * temp = U^(-1) * r
 * U est la pratie triangulaire supérieur de A 
 * 
 * 
 * 
 */ 

{
	int i, j;
	for(i = n-1; i>= 0; i--){
		temp[i] = r[i];
		for(j = ia[i+1] - 1; j >= ia[i]; j--){
			if(ja[j]>i){
				temp[i] -= a[j]*temp[ja[j]];
			} else if(ja[j]==i){
				temp[i] /= a[j]; 
			}	
		}	
	}
}
