#ifndef BOUNDARY_H
#define BOUNDARY_H

/* Wendet die homogene Neumann-Randbedingung Laplace(p|Gamma)=0 
 * auf eine Poissonmatrix p mit (imax x jmax) inneren Gitterpunkten mit Ghost-Schicht an. */
void applyHomogenousNeumannBC(double *p, int imax, int jmax);

void setBoundaryCont(double *U, double *V, int imax, int jmax);

#endif
