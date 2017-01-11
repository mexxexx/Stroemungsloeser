#ifndef VISUAL_H
#define VISUAL_H

void printScalarField(double *field, int imax, int jmax, double xlength, double ylength, char *filename);
void printVectorField(double *xvalues, double *yvalues, int imax, int jmax, double xlength, double ylength, char *filename);
void printObstacles(char *field, int imax, int jmax, double xlength, double ylength, char *filename);
void COMP_PSI_ZETA(double *U, double *V, int imax, int jmax, double xlength, double ylength, double *PSI, double *ZETA, char *FLAG);
void COMPT_HEAT(double *U, double*V, double *TEMP, double *HEAT, double *FLAG, double Re, double Pr, int imax, int jmax, double delx, double dely);

#endif
