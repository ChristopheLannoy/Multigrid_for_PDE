void fine_to_coarse(int m_fine, double *r_fine, double *r_coarse);

void coarse_to_fine(int m_coarse, double *x_fine, double *x_coarse);

int elem_number(int m, int x, int y);

int smoothing(int *ia, int *ja, double *a,double *b, double *u, int n_iter, int n); 
