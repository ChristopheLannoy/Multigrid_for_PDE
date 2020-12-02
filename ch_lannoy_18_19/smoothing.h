/* prototype */
void mat_vec(double *vx, double *vy, int n, int *ia, int *ja, double *a);

int diagonal(int *ia, int *ja, double *a , double **d,int n);
          
void solve_L(int n,int *ia, int *ja, double *a, double *r, double *temp);
void solve_U(int n,int *ia, int *ja, double *a, double *r, double *temp);
