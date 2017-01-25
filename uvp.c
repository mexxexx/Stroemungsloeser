#include <stdio.h>
#include <math.h>
#include "init.h"
#include "boundary.h"
#include "uvp.h"

void solvePoisson(double *p, double *rhs, char *FLAG, double omega, double epsilon, int itermax, 
						double deltaX, double deltaY, int imax, int jmax, int numFluidCells, int* iter, double *res) {
	const double oneMinusOmega = 1 - omega;
	const double oneOverDeltaXSquared = 1 / (deltaX * deltaX);
	const double oneOverDeltaYSquared = 1 / (deltaY * deltaY);
	const double omegaRelaxation = omega / (2 * (oneOverDeltaXSquared + oneOverDeltaYSquared));
	const int imaxPlus2 = imax+2;
	
	int ijlocation = 0;
	
	*res = 1;
	double sum = 0;
	double twoPij = 0;
	*iter = 0;
	int i,j;
	while (*res >= epsilon && *iter < itermax) {
		*res = 0;
		applyHomogenousNeumannBC(p, imax, jmax);
		obstacleBC(p, FLAG, imax, jmax, deltaX, deltaY);
			
		for (j = 1; j <= jmax; j++) {
			for (i = 1; i <= imax; i++) {				
				ijlocation = POS2D(i, j, imaxPlus2);
				if (!FLAG[ijlocation]) {
					p[ijlocation] = oneMinusOmega * p[ijlocation] + omegaRelaxation * (
							oneOverDeltaXSquared * (p[ijlocation-1] + p[ijlocation+1]) + 
							oneOverDeltaYSquared * (p[ijlocation-imaxPlus2] + p[ijlocation+imaxPlus2]) - 
							rhs[ijlocation]);
				}				
			}
		}
		
		for (i = 1; i <= imax; i++) {
			for (j = 1; j <= jmax; j++) {
				ijlocation = POS2D(i, j, imaxPlus2);
				if (!FLAG[ijlocation]) {
					twoPij = 2 * p[ijlocation];
					sum = oneOverDeltaXSquared * (p[ijlocation+1] + p[ijlocation-1] - twoPij) +
						oneOverDeltaYSquared * (p[ijlocation+imaxPlus2] + p[ijlocation-imaxPlus2] - twoPij) - 
						rhs[ijlocation];
					*res += sum * sum;
				}
			}
		}
		*res /= numFluidCells;
		*res = sqrt(*res);
		
		(*iter)++;
	}
	
	if (*iter == itermax) 
		printf("Abgebrochen nach %i iterationen\n", itermax);
}

double min (double a, double b) {
	double min = a;
	if (b < min) min = b;
	return min;
}

void computeDelt(double *delt, double delt_min, int imax, int jmax, double delx, double dely, double umax, double vmax, double Re, double Pr, double tau) {
	if (tau >= 0) {
		(*delt) = min(0.5 * Re * 1/(1/(delx*delx) + 1/(dely*dely)), min(delx/umax, min(dely/vmax, 0.5 * Re * Pr * 1/(1/(delx*delx) + 1/(dely*dely)))));
		(*delt) = min(*delt, delt_min);
		(*delt) *= tau;
	}
}

void computeFG(double *U, double *V, double *F, double *G, double *TEMP, char *FLAG, int imax,  int jmax, double delt, 
				double delx, double dely, double GX, double GY, double alpha, double Re, double beta) {
	double oneOverDelX = 1 / delx;	
	double oneOverDelY = 1 / dely;
	
	double uij, uiPj, uiMj, uijP, uijM, uiMjP;
	double vij, viPj, viMj, vijP, vijM, viPjM;
	double duux, duvy, ddux, dduy;
	double duvx, dvvy, ddvx, ddvy;
	
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			switch(FLAG[POS2D(i,j,imax+2)]){	
				case 29:		//29 entspricht B_N Nord-Kantenzelle
					G[POS2D(i,j,imax+2)]=V[POS2D(i,j,imax+2)];
					break;
				case 27:		//27 entspricht B_S Sued-Kantenzelle
					G[POS2D(i,j-1,imax+2)]=V[POS2D(i,j-1,imax+2)];
					break;	
				case 23:		//23 entspricht B_W West-Kantenzelle
					F[POS2D(i-1,j,imax+2)]=U[POS2D(i-1,j,imax+2)];
					break;		
				case 15:		//15 entspricht B_O Ost-Kantenzelle
					F[POS2D(i,j,imax+2)]=U[POS2D(i,j,imax+2)];
					break;		
				case 13:		//13 entspricht B_NO Nord-Ost-Kantenzelle
					F[POS2D(i,j,imax+2)]=U[POS2D(i,j,imax+2)];
					G[POS2D(i,j,imax+2)]=V[POS2D(i,j,imax+2)];
					break;			
				case 11:		//11 entspricht B_SO Sued-Ost-Kantenzelle
					F[POS2D(i,j,imax+2)]=U[POS2D(i,j,imax+2)];
					G[POS2D(i,j-1,imax+2)]=V[POS2D(i,j-1,imax+2)];
					break;	
				case 21:		//21 entspricht B_NW Nord-West-Kantenzelle
					F[POS2D(i-1,j,imax+2)]=U[POS2D(i-1,j,imax+2)];
					G[POS2D(i,j,imax+2)]=V[POS2D(i,j,imax+2)];
					break;
				case 19:		//19 entspricht B_SW Sued-West-Kantenzelle
					F[POS2D(i-1,j,imax+2)]=U[POS2D(i-1,j,imax+2)];
					G[POS2D(i,j-1,imax+2)]=V[POS2D(i,j-1,imax+2)];
					break;
				default:
					uij = U[POS2D(i,j,imax+2)];
					uiPj = U[POS2D(i+1,j,imax+2)];
					uiMj = U[POS2D(i-1,j,imax+2)];
					uijP = U[POS2D(i,j+1,imax+2)];
					uijM = U[POS2D(i,j-1,imax+2)];
					uiMjP = U[POS2D(i-1,j+1,imax+2)];
					
					vij = V[POS2D(i,j,imax+2)];
					viPj = V[POS2D(i+1,j,imax+2)];
					viMj = V[POS2D(i-1,j,imax+2)];
					vijP = V[POS2D(i,j+1,imax+2)];
					vijM = V[POS2D(i,j-1,imax+2)];
					viPjM = V[POS2D(i+1,j-1,imax+2)];
					
					if (i < imax && (!FLAG[POS2D(i, j, imax+2)] && !FLAG[POS2D(i+1, j, imax+2)])) {
						ddux = oneOverDelX * oneOverDelX * (uiPj - uij - uij + uiMj);			
						dduy = oneOverDelY * oneOverDelY * (uijP - uij - uij + uijM);
						duux = oneOverDelX * 0.25 * (((uij+uiPj)*(uij+uiPj) - (uiMj+uij)*(uiMj+uij)) +
													alpha * (fabs(uij+uiPj)*(uij-uiPj) - fabs(uiMj+uij)*(uiMj-uij)));
						duvy = oneOverDelY * 0.25 * (((vij+viPj)*(uij+uijP) - (vijM+viPjM)*(uijM+uij)) +
													alpha * (fabs(vij+viPj)*(uij-uijP) - fabs(vijM+viPjM)*(uijM-uij)));
								
						F[POS2D(i,j,imax+2)] = uij + delt * (1/Re * (ddux+dduy) - duux - duvy + GX);
					}
					
					if (j < jmax && (!FLAG[POS2D(i, j, imax+2)] && !FLAG[POS2D(i, j+1, imax+2)])) {				
						ddvx = oneOverDelX * oneOverDelX * (viPj - vij - vij + viMj);			
						ddvy = oneOverDelY * oneOverDelY * (vijP - vij - vij + vijM);
						dvvy = oneOverDelY * 0.25 * (((vij+vijP)*(vij+vijP)-(vijM+vij)*(vijM+vij)) +
													alpha * (fabs(vij+vijP)*(vij-vijP)-fabs(vijM+vij)*(vijM-vij)));
						duvx = oneOverDelX * 0.25 * (((uij+uijP)*(vij+viPj) - (uiMj+uiMjP)*(viMj+vij)) + 
													alpha * (fabs(uij+uijP)*(vij-viPj) - fabs(uiMj+uiMjP)*(viMj-vij)));
								
						G[POS2D(i,j,imax+2)] = vij + delt * (1/Re * (ddvx+ddvy) - duvx - dvvy + GY);
					}
					break;
			}
		}
	}

	for (int j = 1; j <= jmax; j++) {
		F[POS2D(0,j,imax+2)] = U[POS2D(0,j,imax+2)];
		F[POS2D(imax,j,imax+2)] = U[POS2D(imax,j,imax+2)];
	}
	
	for (int i = 1; i <= imax; i++) {
		G[POS2D(i,0,imax+2)] = V[POS2D(i, 0,imax+2)];
		G[POS2D(i,jmax,imax+2)] = V[POS2D(i,jmax,imax+2)];
	}

/*Erweiterung um Boussinesq-Term nach Energiegleichung: */

	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {	
			F[POS2D(i,j,imax+2)] -= beta*0.5*delt*(TEMP[POS2D(i,j,imax+2)]+TEMP[POS2D(i+1,j,imax+2)])*GX;
			G[POS2D(i,j,imax+2)] -= beta*0.5*delt*(TEMP[POS2D(i,j,imax+2)]+TEMP[POS2D(i,j+1,imax+2)])*GY;
		}
	}
}

void computeRHS(double *F, double *G, double *rhs, char *FLAG, int imax, int jmax, double delt, double delx, double dely) {
	double oneOverDelX = 1 / delx;	
	double oneOverDelY = 1 / dely;
	
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			int ijlocation = POS2D(i, j, imax+2);
			if (!FLAG[ijlocation]) {
				rhs[ijlocation] = 1/delt * 
					(oneOverDelX * (F[ijlocation]-F[ijlocation-1]) +
					oneOverDelY * (G[ijlocation]-G[ijlocation-imax-2]));
			}
		}
	}
}

void adapUV(double *U, double *V, double *F, double *G, double *P, char *FLAG, int imax, 
					int jmax, double delt, double delx, double dely, double *umax, double *vmax) {	
	double oneOverDelX = 1 / delx;	
	double oneOverDelY = 1 / dely;
	double u, v;
	*umax = *vmax = -1;
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			if (i < imax && (!FLAG[POS2D(i, j, imax+2)] && !FLAG[POS2D(i+1, j, imax+2)])) {
				u = F[POS2D(i, j, imax+2)] - delt * oneOverDelX * (P[POS2D(i+1, j, imax+2)] - P[POS2D(i, j, imax+2)]);
				U[POS2D(i, j, imax+2)] = u;
				u = fabs(u);
				if (u > *umax)
					*umax = u;
			}
			if (j < jmax && (!FLAG[POS2D(i, j, imax+2)] && !FLAG[POS2D(i, j+1, imax+2)])) {
				v = G[POS2D(i, j, imax+2)] - delt * oneOverDelY * (P[POS2D(i, j+1, imax+2)] - P[POS2D(i, j, imax+2)]);
				V[POS2D(i, j, imax+2)] = v; 
				v = fabs(v);
				if (v > *vmax)
					*vmax = v;
			}
		}
	}
}

void computeTEMP(double *U, double *V, double *TEMP, char *FLAG,  int imax, 
					int jmax, double delt, double delx, double dely, double alpha, double Re, double Pr){

	double duTx, dvTy, ddTxx, ddTyy, dTt;
	double oneOverDelX = 1 / delx;	
	double oneOverDelY = 1 / dely;
	double uij, uiMj, vij, vijM, Tij, TiPj, TiMj, TijP, TijM;

	for(int i=1; i<=imax; i++){
		for(int j=1; j<=jmax; j++){
			uij = U[POS2D(i, j, imax+2)];
			uiMj = U[POS2D(i-1, j, imax+2)];
			vij = V[POS2D(i, j, imax+2)];
			vijM = V[POS2D(i, j-1, imax+2)];
			
			Tij = TEMP[POS2D(i, j, imax+2)];
			TiPj = TEMP[POS2D(i+1, j, imax+2)];
			TiMj = TEMP[POS2D(i-1, j, imax+2)];
			TijP = TEMP[POS2D(i, j+1, imax+2)];
			TijM = TEMP[POS2D(i, j-1, imax+2)];

			//duTx=(oneOverDelX*(U[POS2D(i, j, imax+2)]*((TEMP[POS2D(i, j, imax+2)]+TEMP[POS2D(i+1, j, imax+2)])*0.5)-U[POS2D(i-1, j, imax+2)]*((TEMP[POS2D(i-1, j, imax+2)]+TEMP[POS2D(i, j, imax+2)])*0.5)))+alpha*oneOverDelX*(fabs(U[POS2D(i, j, imax+2)])*((TEMP[POS2D(i, j, imax+2)]-TEMP[POS2D(i+1, j, imax+2)])*0.5)-fabs(U[POS2D(i-1, j, imax+2)])*((TEMP[POS2D(i-1, j, imax+2)]-TEMP[POS2D(i, j, imax+2)])*0.5));
			duTx = oneOverDelX * 0.5 * (uij*(Tij+TiPj)-uiMj*(TiMj+Tij)+alpha*(fabs(uij)*(Tij-TiPj)-fabs(uiMj)*(TiMj-Tij)));
			dvTy = oneOverDelY * 0.5 * (vij*(Tij+TijP)-vijM*(TijM+Tij)+alpha*(fabs(vij)*(Tij-TijP)-fabs(vijM)*(TijM-Tij)));
			
			//dvTy=(oneOverDelY*(V[POS2D(i, j, imax+2)]*((TEMP[POS2D(i, j, imax+2)]+TEMP[POS2D(i, j+1, imax+2)])*0.5)-V[POS2D(i, j-1, imax+2)]*((TEMP[POS2D(i, j-1, imax+2)]+TEMP[POS2D(i, j, imax+2)])*0.5)))+alpha*oneOverDelY*(fabs(V[POS2D(i, j, imax+2)])*((TEMP[POS2D(i, j, imax+2)]-TEMP[POS2D(i, j+1, imax+2)])*0.5)-fabs(V[POS2D(i, j-1, imax+2)])*((TEMP[POS2D(i, j-1, imax+2)]-TEMP[POS2D(i, j, imax+2)])*0.5));

			//ddTxx=oneOverDelX*oneOverDelX*(TEMP[POS2D(i+1, j, imax+2)]-2*TEMP[POS2D(i, j, imax+2)]+TEMP[POS2D(i-1, j, imax+2)]);
			ddTxx = oneOverDelX*oneOverDelX*(TiPj-2*Tij+TiMj);
			ddTyy = oneOverDelY*oneOverDelY*(TijP-2*Tij+TijM);
			
			//ddTyy=oneOverDelY*oneOverDelY*(TEMP[POS2D(i, j+1, imax+2)]-2*TEMP[POS2D(i, j, imax+2)]+TEMP[POS2D(i, j-1, imax+2)]);

			dTt=(1./Re)*(1./Pr)*(ddTxx+ddTyy)-duTx-dvTy;
	
			//TEMP[POS2D(i+1, j, imax+2)]=delt*dTt+TEMP[POS2D(i+1, j, imax+2)];
			TEMP[POS2D(i, j, imax+2)] += delt*dTt;
		}
	}

}
