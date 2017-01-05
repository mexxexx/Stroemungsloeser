#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "particles.h"
#include "init.h"

void particleVelocity(double *U, double *V, double delx, double dely, int imax, int jmax, Particle *particles, int partCount) {
	double oneOverDelxDely = 1.0 / (delx * dely);
	for (int k = 0; k < partCount; k++) {
		if (!particles[k].isActive) continue;
		
		double x = particles[k].x;
		double y = particles[k].y;
		
		int i = (int)(x/delx)+1;
		int j = (int)((y+0.5*dely)/dely)+1;
		
		double x1 = (i-1)*delx;
		double x2 = i*delx;
		double y1 = (j-1.5)*delx;
		double y2 = (j-0.5)*delx;
		
		double u1 = U[POS2D(i-1, j-1, imax+2)];
		double u2 = U[POS2D(i, j-1, imax+2)];
		double u3 = U[POS2D(i-1, j, imax+2)];
		double u4 = U[POS2D(i, j, imax+2)];
		
		particles[k].u = oneOverDelxDely * ((x2-x)*(y2-y)*u1+(x-x1)*(y2-y)*u2+(x2-x)*(y-y1)*u3+(x-x1)*(y-y1)*u4);
		
		i=(int)((x+0.5*delx)/delx)+1;
		j=(int)(y/dely)+1;
		
		x1 = (i-1.5)*delx;
		x2 = (i-0.5)*delx;
		y1 = (j-1)*delx;
		y2 = j*delx;
		
		double v1 = V[POS2D(i-1, j-1, imax+2)];
		double v2 = V[POS2D(i, j-1, imax+2)];
		double v3 = V[POS2D(i-1, j, imax+2)];
		double v4 = V[POS2D(i, j, imax+2)];
		
		particles[k].v = oneOverDelxDely * ((x2-x)*(y2-y)*v1+(x-x1)*(y2-y)*v2+(x2-x)*(y-y1)*v3+(x-x1)*(y-y1)*v4);
	}
}

void particleTransport(Particle *particles, double delt, int partCount, double xlength, double ylength) {
	for (int k = 0; k < partCount; k++) {
		Particle *p = &particles[k];
		if (!p->isActive) continue;
		
		p->x += delt * p->u;
		p->y += delt * p->v;
		
		if (p->x < 0 || p->x >= xlength || p->y < 0 || p->y >= ylength) {
			p->isActive = 0;
			p->x = -10;
			p->y = -10;
		}
	}
}


void particleSeed(Particle *particles, double posx1, double posy1, double posx2, double posy2, int partCount, int anzahl){
	int count=0;
	double horizontal=(1./anzahl)*fabs(posx2-posx1);
	double vertical=(1./anzahl)*fabs(posy1-posy2);
	for (int i=0;i<partCount;i++) {
		if(count>anzahl) break;

		if(particles[i].isActive) continue;

		else{
			particles[i].x=posx1+count*horizontal;
			particles[i].y=posy1+count*vertical;
			particles[i].isActive=1;	
			count++;	
		}
	}
}


void particleInit(Particle *particles, int partCount){
	for(int i=0;i<partCount;i++){	
		particles[i].x=-10;
		particles[i].y=-10;
		particles[i].u=0;
		particles[i].v=0;
		particles[i].isActive=0;
	}
}

void printParticles(Particle *particles, int partCount, char *filename) {
	FILE *f = fopen(filename, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		return;
	}
	
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Particles\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET POLYDATA\n");
	fprintf(f, "POINTS %i double\n", partCount);
	for (int i = 0; i < partCount; i++)
		fprintf(f, "%f %f 0.0\n", particles[i].x, particles[i].y);
	fprintf(f, "-10.0 -10.0 0.0\n");
	
	fclose(f);
}


void bmp_verify(char* FLAG, int imax, int jmax){
	int sides=0;

	for (int i = 1; i <= imax; i++) {  
		for (int j = 1; j <= jmax; j++) {
			//int posij = POS2D(i, j, imax+2);
			sides=0;
			if (FLAG[POS2D(i, j, imax+2)]!=0) {
				if (FLAG[POS2D(i+1, j, imax+2)]!=0) sides++;
				if (FLAG[POS2D(i-1, j, imax+2)]!=0) sides++;
				if (FLAG[POS2D(i, j+1, imax+2)]!=0) sides++;
				if (FLAG[POS2D(i, j-1, imax+2)]!=0) sides++;
			}
			if(sides>4){
				printf("Oops. There appears to be a forbidden cell in your obstacle. It's at position i=%d, j=%d\n",i,j);
				exit(-1);
			}
		}
	}
}
