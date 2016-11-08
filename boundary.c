#include "init.h"
#include "boundary.h"

/* p ist Gitter mit (imax x jmax) inneren Zellen */
void applyHomogenousNeumannBC(double *p, int imax, int jmax) {
	const int imaxPlus2 = imax+2;
	const int imaxPlus1 = imax+1;
	const int jmaxPlus1 = jmax+1;
	const int lowerBound = (imaxPlus1 < jmaxPlus1) ? imaxPlus1 : jmaxPlus1;
	int k;
	
	for (k = 1; k < lowerBound; k++) {
		p[POS2D(k, 0, imaxPlus2)] = p[POS2D(k, 1, imaxPlus2)];
		p[POS2D(k, jmaxPlus1, imaxPlus2)] = p[POS2D(k, jmax, imaxPlus2)];
		p[POS2D(0, k, imaxPlus2)] = p[POS2D(1, k, imaxPlus2)];
		p[POS2D(imaxPlus1, k, imaxPlus2)] = p[POS2D(imax, k, imaxPlus2)];
	}
	
	for (k = lowerBound; k < imaxPlus1; k++) {
		p[POS2D(k, 0, imaxPlus2)] = p[POS2D(k, 1, imaxPlus2)];
		p[POS2D(k, jmaxPlus1, imaxPlus2)] = p[POS2D(k, jmax, imaxPlus2)];
	}
	
	for (k = lowerBound; k < jmaxPlus1; k++) {
		p[POS2D(0, k, imaxPlus2)] = p[POS2D(1, k, imaxPlus2)];
		p[POS2D(imaxPlus1, k, imaxPlus2)] = p[POS2D(imax, k, imaxPlus2)];
	}
}

void setBoundaryCont(double *U, double *V, int imax, int jmax) {
	for (int i = 1; i <= imax; i++) {
		V[POS2D(i, 0, imax+2)] = 0;
		V[POS2D(i, jmax, imax+2)] = 0;
		U[POS2D(i, 0, imax+2)] = -U[POS2D(i, 1, imax+2)];
		U[POS2D(i, jmax+1, imax+2)] = -U[POS2D(i, jmax, imax+2)];
	}
	for (int j = 1; j <= jmax; j++) {
		U[POS2D(0, j, imax+2)] = 0;
		U[POS2D(imax, j, imax+2)] = 0;
		V[POS2D(0, j, imax+2)] = -V[POS2D(1, j, imax+2)];
		V[POS2D(imax+1, j, imax+2)] = -V[POS2D(imax, j, imax+2)];
	}
}
