#ifndef my_header_stuff     //La modification de ces données implique une recompilation complete
#define my_header_stuff     //    !!!! make clean ==>  make !!!!!
    #define x_ratio 5.0/8   // point that define the problem shape
	#define y_ratio 2.0/8   //  (remove the upper right rectangle)
	#define L 2.0           // problem dimension (square's sides lenght [m])
#endif

/* prototypes utilisés dans main.h */
double mytimer();
int prob(int m, int *n, int **ia, int **ja, double **a, double **b);

int plot(int m, double *evecs);
int plot_vector(int n, double *x);
int plot_eigen(int m, double *evecs);
                 
void mat_vec(double *vx, double *vy, int n, int *ia, int *ja, double *a);

int smoothing(int *ia, int *ja, double *a,double *b, double *u, int n_iter, int n); 

int diagonal(int *ia, int *ja, double *a , double **d,int n);

double norm_residu( int n, double *b, double *x, int *ia,
                  int *ja, double *a);
double norme(int n, double *x);
                  
void solve_L(int n,int *ia, int *ja, double *a, double *r, double *temp);
void solve_U(int n,int *ia, int *ja, double *a, double *r, double *temp);

int multi_grid(int m_fine, int n_fine, int *ia, int *ja, double *a, 
				double *b, double *u, double m_coarsest_level);

int primme(int primme_n, int* primme_ia, int* primme_ja, double* primme_a, 
           int nev, double *evals, double *evecs, int max);


