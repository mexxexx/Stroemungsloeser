#ifndef UVP_H
#define UVP_H

#ifdef DEBUG
#define PRINT_RES_DEBUG 1
#else
#define PRINT_RES_DEBUG 0
#endif

/* Löst das LGS Laplace(p)=rhs mit einem effizienten SOR Verfahren, 
 * wobei p eine (imax+2) x (jmax+2) Matrix ist die imax bzw. jmax 
 * innere Gitterpunkte und eine Ghost-Schicht besitzt. */
void solvePoisson(double *p, double *rhs, double omega, double epsilon, int itermax, 
					double deltaX, double deltaY, int imax, int jmax);

void computeDelt(double *delt, int imax, int jmax, double delx, double dely, 
					double umax, double vmax, double Re, double tau);

void computeRHS(double *F, double *G, double *rhs, int imax, int jmax, double delt, double delx, double dely);
					
void computeFG(double *U, double *V, double *F, double *G, char *FLAG, int imax, 
					int jmax, double delt, double delx, double dely, double GX, double GY, double alpha, double Re);
					
void adapUV(double *U, double *V, double *F, double *G, double *P, int imax, 
					int jmax, double delt, double delx, double dely, double *umax, double *vmax);
#endif
