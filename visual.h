#ifndef VISUAL_H
#define VISUAL_H

void printScalarField(double *field, int imax, int jmax, double xlength, double ylength, char *filename);
void printVectorField(double *xvalues, double *yvalues, int imax, int jmax, double xlength, double ylength, char *filename);

#endif