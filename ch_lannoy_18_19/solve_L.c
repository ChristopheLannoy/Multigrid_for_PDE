void solve_L(int n,int *ia, int *ja, double *a, double *r, double *temp)

/*
 * Solveur triangulaire inférieur
 * temp = L^(-1) * r
 * L est la pratie triangulaire inférieur de A 
 * 
 * 
 * 
 */ 

{
	int i, j;
	for(i = 0; i < n; i++){
		temp[i] = r[i];
		for(j = ia[i]; j < ia[i+1]; j++){
			if(ja[j]<i){
				temp[i] -= a[j]*temp[ja[j]];
			} else if(ja[j]==i){
				temp[i] /= a[j]; 
			}	
		}	
	}
}
