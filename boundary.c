#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "init.h"
#include "boundary.h"

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

void setBoundaryCond(double *U, double *V, double *TEMP, char *FLAG, double delx, double dely, int imax, int jmax,int wl, int wr, int wt, int wb,
					int tl, double tl_value, int tr, double tr_value, int tt, double tt_value, int tb, double tb_value) {
	double oneOverdelxSquaredPlusdelySquared=(1/(delx*delx+dely*dely));
	
	switch(wl){
		case 1:
			for (int j = 1; j <= jmax; j++) {
				U[POS2D(0, j, imax+2)] = 0;
				V[POS2D(0, j, imax+2)] = -V[POS2D(1, j, imax+2)];
			}
			break;
		case 2:
			for (int j = 1; j <= jmax; j++) {
				U[POS2D(0,j,imax+2)]=0;
				V[POS2D(0,j,imax+2)]=V[POS2D(1,j,imax+2)];				
			}
			break;
		case 3:
			for (int j = 1; j <= jmax; j++) {
				U[POS2D(0,j,imax+2)]=U[POS2D(1,j,imax+2)];
				V[POS2D(0,j,imax+2)]=V[POS2D(1,j,imax+2)];				
			}
			break;
		default: 
			printf("Randbedingung nicht zulaessig.\n"); 
			exit(EXIT_FAILURE);
			break;
	}
	
	switch(tl){
		case 1:
			for (int j = 1; j <= jmax; j++)
				TEMP[POS2D(0, j, imax+2)] = 2*tl_value*((j-0.5)*dely) - TEMP[POS2D(1, j, imax+2)];
			break;
		case 2:
			for (int j = 1; j <= jmax; j++) 
				TEMP[POS2D(0, j, imax+2)] = TEMP[POS2D(1, j, imax+2)] + delx*tl_value*((j-0.5)*dely);
			break;
		default: 
			printf("Randbedingung nicht zulaessig.\n"); 
			exit(EXIT_FAILURE);
			break;
	}

	switch(wr){
		case 1:
			for (int j = 1; j <= jmax; j++) {
				U[POS2D(imax, j, imax+2)] = 0;
				V[POS2D(imax+1, j, imax+2)] = -V[POS2D(imax, j, imax+2)];
			}
			break;
		case 2:		
			for (int j = 1; j <= jmax; j++) {
				U[POS2D(imax, j, imax+2)] = 0;
				V[POS2D(imax+1, j, imax+2)] = V[POS2D(imax, j, imax+2)];
			}
			break;
		case 3:
			for (int j = 1; j <= jmax; j++) {
				U[POS2D(imax, j, imax+2)] = U[POS2D(imax-1, j, imax+2)];
				V[POS2D(imax+1, j, imax+2)] = V[POS2D(imax, j, imax+2)];
			}
			break;
		default: 
			printf("Randbedingung nicht zulaessig.\n"); 
			exit(EXIT_FAILURE);
			break;
	}

	switch(tr){
		case 1:
			for (int j = 1; j <= jmax; j++)
				TEMP[POS2D(imax+1, j, imax+2)] = 2*tl_value*((j-0.5)*dely) - TEMP[POS2D(imax, j, imax+2)];
			break;
		case 2:
			for (int j = 1; j <= jmax; j++) 
				TEMP[POS2D(imax+1, j, imax+2)] = TEMP[POS2D(imax, j, imax+2)] + delx*tl_value*((j-0.5)*dely);
			break;
		default: 
			printf("Randbedingung nicht zulaessig.\n"); 
			exit(EXIT_FAILURE);
			break;
	}

	switch(wt){
		case 1:
			for (int i = 1; i <= imax; i++) {
				V[POS2D(i, jmax, imax+2)] = 0;
				U[POS2D(i, jmax+1, imax+2)] = -U[POS2D(i, jmax, imax+2)];
			}
			break;
		case 2:
			for (int i = 1; i <= imax; i++) {
				V[POS2D(i, jmax, imax+2)] = 0;
				U[POS2D(i, jmax+1, imax+2)] = U[POS2D(i, jmax, imax+2)];
			}
			break;
		case 3:
			for (int i = 1; i <= imax; i++) {
				V[POS2D(i, jmax, imax+2)] = V[POS2D(i, jmax-1, imax+2)];
				U[POS2D(i, jmax+1, imax+2)] = U[POS2D(i, jmax, imax+2)];
			}
			break;		
		default: 
			printf("Randbedingung nicht zulaessig.\n"); 
			exit(EXIT_FAILURE);
			break;	
	}

	switch(tt){
		case 1:
			for (int i = 1; i <= imax; i++)
				TEMP[POS2D(i, 0, imax+2)] = 2*tl_value*((i-0.5)*delx) - TEMP[POS2D(i, 1, imax+2)];
			break;
		case 2:
			for (int i = 1; i <= imax; i++)
				TEMP[POS2D(i, 0, imax+2)] = TEMP[POS2D(i, 1, imax+2)] + dely*tl_value*((i-0.5)*delx);
			break;
		default: 
			printf("Randbedingung nicht zulaessig.\n"); 
			exit(EXIT_FAILURE);
			break;
	}
	
	switch(wb){
		case 1:
			for (int i = 1; i <= imax; i++) {
				V[POS2D(i, 0, imax+2)] = 0;
				U[POS2D(i, 0, imax+2)] = -U[POS2D(i, 1, imax+2)];
			}
			break;
		case 2:
			for (int i = 1; i <= imax; i++) {
				V[POS2D(i, 0, imax+2)] = 0;
				U[POS2D(i, 0, imax+2)] = U[POS2D(i, 1, imax+2)];
			}
			break;
		case 3:
			for (int i = 1; i <= imax; i++) {
				V[POS2D(i, 0, imax+2)] = V[POS2D(i, 1, imax+2)];
				U[POS2D(i, 0, imax+2)] = U[POS2D(i, 1, imax+2)];
			}
			break;
		default: 
			printf("Randbedingung nicht zulaessig.\n"); 
			exit(EXIT_FAILURE);
			break;
	}

	switch(tb){
		case 1:
			for (int i = 1; i <= imax; i++)
				TEMP[POS2D(i, jmax+1, imax+2)] = 2*tl_value*((i-0.5)*delx) - TEMP[POS2D(i, jmax, imax+2)];
			break;
		case 2:
			for (int i = 1; i <= imax; i++)
				TEMP[POS2D(i, jmax+1, imax+2)] = TEMP[POS2D(i, jmax, imax+2)] + dely*tl_value*((i-0.5)*delx);
			break;
		default: 
			printf("Randbedingung nicht zulaessig.\n"); 
			exit(EXIT_FAILURE);
			break;
	}
	
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			switch(FLAG[POS2D(i,j,imax+2)]){
				case 29:		//29 entspricht B_N Nord-Kantenzelle
					U[POS2D(i,j,imax+2)]=-U[POS2D(i,j+1,imax+2)];
					U[POS2D(i-1,j,imax+2)]=-U[POS2D(i-1,j+1,imax+2)];
					V[POS2D(i,j,imax+2)]=0;
					TEMP[POS2D(i,j,imax+2)] = TEMP[POS2D(i,j+1,imax+2)];
					break;
				case 27:		//27 entspricht B_S Sued-Kantenzelle
					U[POS2D(i,j,imax+2)]=-U[POS2D(i,j-1,imax+2)];
					U[POS2D(i-1,j,imax+2)]=-U[POS2D(i-1,j-1,imax+2)];
					V[POS2D(i,j-1,imax+2)]=0;
					TEMP[POS2D(i,j,imax+2)] = TEMP[POS2D(i,j-1,imax+2)];
					break;	
				case 23:		//23 entspricht B_W West-Kantenzelle
					V[POS2D(i,j,imax+2)]=-V[POS2D(i-1,j,imax+2)];
					V[POS2D(i,j-1,imax+2)]=-V[POS2D(i-1,j-1,imax+2)];
					U[POS2D(i-1,j,imax+2)]=0;
					TEMP[POS2D(i,j,imax+2)] = TEMP[POS2D(i-1,j,imax+2)];
					break;		
				case 15:		//15 entspricht B_O Ost-Kantenzelle
					V[POS2D(i,j,imax+2)]=-V[POS2D(i+1,j,imax+2)];
					V[POS2D(i,j-1,imax+2)]=-V[POS2D(i+1,j-1,imax+2)];
					U[POS2D(i,j,imax+2)]=0;
					TEMP[POS2D(i,j,imax+2)] = TEMP[POS2D(i+1,j,imax+2)];
					break;		
				case 13:		//13 entspricht B_NO Nord-Ost-Kantenzelle
					U[POS2D(i-1,j,imax+2)]=-U[POS2D(i-1,j+1,imax+2)];
					V[POS2D(i,j-1,imax+2)]=-V[POS2D(i+1,j-1,imax+2)];
					U[POS2D(i,j,imax+2)]=0;
					V[POS2D(i,j,imax+2)]=0;
					TEMP[POS2D(i,j,imax+2)]=oneOverdelxSquaredPlusdelySquared*(delx*delx*TEMP[POS2D(i,j+1,imax+2)]+dely*dely*TEMP[POS2D(i+1,j,imax+2)]);
					break;			
				case 11:		//11 entspricht B_SO Sued-Ost-Kantenzelle
					U[POS2D(i-1,j,imax+2)]=-U[POS2D(i-1,j-1,imax+2)];
					V[POS2D(i,j,imax+2)]=-V[POS2D(i+1,j,imax+2)];
					U[POS2D(i,j,imax+2)]=0;
					V[POS2D(i,j-1,imax+2)]=0;
					TEMP[POS2D(i,j,imax+2)]=oneOverdelxSquaredPlusdelySquared*(delx*delx*TEMP[POS2D(i,j-1,imax+2)]+dely*dely*TEMP[POS2D(i+1,j,imax+2)]);
					break;	
				case 21:		//21 entspricht B_NW Nord-West-Kantenzelle
					U[POS2D(i,j,imax+2)]=-U[POS2D(i,j+1,imax+2)];
					V[POS2D(i,j-1,imax+2)]=-V[POS2D(i-1,j-1,imax+2)];
					U[POS2D(i-1,j,imax+2)]=0;
					V[POS2D(i,j,imax+2)]=0;
					TEMP[POS2D(i,j,imax+2)]=oneOverdelxSquaredPlusdelySquared*(delx*delx*TEMP[POS2D(i,j+1,imax+2)]+dely*dely*TEMP[POS2D(i-1,j,imax+2)]);
					break;
				case 19:		//19 entspricht B_SW Sued-West-Kantenzelle
					U[POS2D(i,j,imax+2)]=-U[POS2D(i,j-1,imax+2)];
					V[POS2D(i,j,imax+2)]=-V[POS2D(i-1,j,imax+2)];
					U[POS2D(i-1,j,imax+2)]=0;
					V[POS2D(i,j-1,imax+2)]=0;
					TEMP[POS2D(i,j,imax+2)]=oneOverdelxSquaredPlusdelySquared*(delx*delx*TEMP[POS2D(i,j-1,imax+2)]+dely*dely*TEMP[POS2D(i-1,j,imax+2)]);
					break;
			}
		}
	}		
}


void obstacleBC(double *p, char *FLAG, int imax, int jmax, double deltaX, double deltaY){
	double oneOverDeltaXSquaredPlusDeltaYSquared=(1/(deltaX*deltaX+deltaY*deltaY));

	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			switch(FLAG[POS2D(i,j,imax+2)]){
				case 29:		/*29 entspricht B_N Nord-Kantenzelle*/
					p[POS2D(i,j,imax+2)]=p[POS2D(i,j+1,imax+2)];
					break;
				case 27:		/*27 entspricht B_S Sued-Kantenzelle*/
					p[POS2D(i,j,imax+2)]=p[POS2D(i,j-1,imax+2)];
					break;	
				case 23:		/*23 entspricht B_W West-Kantenzelle*/
					p[POS2D(i,j,imax+2)]=p[POS2D(i-1,j,imax+2)];
					break;		
				case 15:		/*15 entspricht B_O Ost-Kantenzelle*/
					p[POS2D(i,j,imax+2)]=p[POS2D(i+1,j,imax+2)];
					break;		
				case 13:		/*13 entspricht B_NO Nord-Ost-Kantenzelle*/
					p[POS2D(i,j,imax+2)]=oneOverDeltaXSquaredPlusDeltaYSquared*(deltaX*deltaX*p[POS2D(i,j+1,imax+2)]+deltaY*deltaY*p[POS2D(i+1,j,imax+2)]);
					break;			
				case 11:		/*11 entspricht B_SO Sued-Ost-Kantenzelle*/
					p[POS2D(i,j,imax+2)]=oneOverDeltaXSquaredPlusDeltaYSquared*(deltaX*deltaX*p[POS2D(i,j-1,imax+2)]+deltaY*deltaY*p[POS2D(i+1,j,imax+2)]);
					break;	
				case 21:		/*21 entspricht B_NW Nord-West-Kantenzelle*/
					p[POS2D(i,j,imax+2)]=oneOverDeltaXSquaredPlusDeltaYSquared*(deltaX*deltaX*p[POS2D(i,j+1,imax+2)]+deltaY*deltaY*p[POS2D(i-1,j,imax+2)]);
					break;
				case 19:		/*19 entspricht B_SW Sued-West-Kantenzelle*/
					p[POS2D(i,j,imax+2)]=oneOverDeltaXSquaredPlusDeltaYSquared*(deltaX*deltaX*p[POS2D(i,j-1,imax+2)]+deltaY*deltaY*p[POS2D(i-1,j,imax+2)]);
					break;
			}
		}
	}		
}


void initDrivenCavity(double *U, double *V, int imax, int jmax) {
	for (int i = 1; i <= imax; i++) {
		U[POS2D(i, jmax+1, imax+2)] = (2.0 - U[POS2D(i, jmax, imax+2)]);
		//U[POS2D(i, 0, imax+2)] = -(2.0 - U[POS2D(i, 1, imax+2)]);
	}
	/*for (int j = 1; j <= jmax; j++) {
		V[POS2D(0, j, imax+2)] = (2.0 - V[POS2D(1, j, imax+2)]);
		V[POS2D(imax+1, j, imax+2)] = -(2.0 - V[POS2D(imax, j, imax+2)]);
	}*/
}

void initWest(double *U, double *V, int imax, int jmax) {
	for (int j = 1; j <= jmax; j++) {
		U[POS2D(0,j,imax+2)]=1.0; 
		V[POS2D(0,j,imax+2)]=-V[POS2D(1,j,imax+2)];	//Links einströmen
	}
}



void initEast(double *U, double *V, int imax, int jmax) {
	for (int j = 1; j <= jmax; j++) {
		U[POS2D(imax,j,imax+2)]=-1.0; 
		V[POS2D(imax+1,j,imax+2)]= -V[POS2D(imax,j,imax+2)];	//Rechts einströmen	
	}
}

void initNorth(double *U, double *V, int imax, int jmax) {
	for (int i = 1; i <= imax; i++) {
		V[POS2D(i,jmax,imax+2)] = -1.0; 
		U[POS2D(i,jmax+1,imax+2)] = -U[POS2D(i,jmax,imax+2)];	//Oben einströmen	
	}
}

void initSouth(double *U, double *V, int imax, int jmax) {
	for (int i = 1; i <= imax; i++) {
		V[POS2D(i,0,imax+2)] = 1.0; 
		U[POS2D(i,0,imax+2)] = -U[POS2D(i,1,imax+2)];	//Unten einströmen	
	}
}

void setSpecialBoundaryCond (double *U, double*V, double *TEMP, int imax, int jmax, char *problem){
	if((strcmp(problem,"Driven cavity") == 0) || (strcmp(problem,"DC") == 0) || (strcmp(problem,"driven cavity") == 0)){
		initDrivenCavity(U, V, imax, jmax);
	}
	else if(strcmp(problem, "West") == 0 || (strcmp(problem, "west") == 0)){
		initWest(U, V, imax, jmax);
	}
	else if(strcmp(problem, "East") == 0 || (strcmp(problem, "east") == 0)){
		initEast(U, V, imax, jmax);
	}
	else if(strcmp(problem, "North") == 0 || (strcmp(problem, "north") == 0)){
		initNorth(U, V, imax, jmax);
	}
	else if(strcmp(problem, "South") == 0 || (strcmp(problem, "south") == 0)){
		initSouth(U, V, imax, jmax);
	}
	else if(strcmp(problem, "Stufe") == 0){
		initWest(U, V, imax, jmax);
	}
	else
		printf("Unbekanntes Problem");	
}
